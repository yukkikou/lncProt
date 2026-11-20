#!/usr/bin/env python
"""47_strand_controls.py — strict strand-parity controls by lncRNA source.

Modes:
  S sources   — plus/minus balance by class, LTR status, age, activity,
                with PCG reference, Wilson CI and per-chromosome tests.
  T te        — TSS-donating TE orientation by order/family/class, with a
                chromosome-matched gene-strand permutation null and logistic
                model.
  C cage      — CAGE-peak plus bias controls: score/length/GC/chromosome
                strata, PCG-TSS proximity, and logistic model.
Outputs -> results/strand_controls/*.tsv + step summary.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest, chi2_contingency, fisher_exact
from statsmodels.formula.api import logit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "strand_controls")
DATA = os.path.join(WORKDIR, "data")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def chrom_het(pos, neg, chrom):
    tab = pd.crosstab(chrom, np.where(pos, "plus", "minus"))
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return np.nan
    try:
        return chi2_contingency(tab.values)[1]
    except Exception:
        return np.nan


def main():
    mode = (sys.argv[1] if len(sys.argv) > 1 else "S").upper()
    meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                       sep="\t", index_col=0, low_memory=False)
    te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                     index_col=0, low_memory=False)[["tss_in_ltr"]]
    act = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                      index_col=0, low_memory=False)[["active"]]
    age = pd.read_csv(os.path.join(RES, "comparative",
                                   "E_transcript_age.tsv"), sep="\t",
                      index_col=0)[["age_class"]]
    m = meta[["chrom", "strand", "cls", "group"]].join(te).join(act) \
        .join(age)

    with step("47_strand_" + mode) as st:
        os.makedirs(OUT, exist_ok=True)
        if mode == "S":
            groups = []
            for c in CLASSES:
                groups.append((c, m.cls == c))
            groups += [("PCG", m.group == "pcg"),
                       ("pseudogene", m.cls == "pseudogene"),
                       ("lnc_LTR", (m.group == "lnc") & (m.tss_in_ltr == 1)),
                       ("lnc_nonLTR",
                        (m.group == "lnc") & (m.tss_in_ltr != 1)),
                       ("intergenic_LTR",
                        (m.cls == "intergenic") & (m.tss_in_ltr == 1)),
                       ("intergenic_nonLTR",
                        (m.cls == "intergenic") & (m.tss_in_ltr != 1)),
                       ("lnc_active", (m.group == "lnc") & (m.active == 1)),
                       ("lnc_inactive",
                        (m.group == "lnc") & (m.active == 0))]
            for a in ["human_specific", "primate", "mammalian", "eutherian"]:
                groups.append(("age_" + a, (m.group == "lnc") &
                               (m.age_class == a)))
            pcg = m[m.group == "pcg"]
            ppos = int((pcg.strand == "+").sum())
            pneg = int((pcg.strand == "-").sum())
            rows = []
            for lab, mask in groups:
                sub = m[mask]
                pos = int((sub.strand == "+").sum())
                neg = int((sub.strand == "-").sum())
                n = pos + neg
                if n < 20:
                    continue
                lo, hi = wilson(pos, n)
                p0 = binomtest(pos, n, 0.5).pvalue
                pf = fisher_exact([[pos, neg], [ppos, pneg]])[1]
                rows.append(dict(source=lab, n=n, plus=pos, minus=neg,
                                 plus_frac=pos / n, ci_low=lo, ci_high=hi,
                                 p_vs_half=p0, p_vs_PCG=pf,
                                 chrom_het_p=chrom_het(
                                     sub.strand == "+", sub.strand == "-",
                                     sub.chrom)))
            T = pd.DataFrame(rows)
            T["fdr_vs_half"] = bh_fdr(T.p_vs_half.values)
            T["fdr_vs_PCG"] = bh_fdr(T.p_vs_PCG.values)
            T.to_csv(safe_write_path(os.path.join(OUT, "S_sources.tsv")),
                     sep="\t", index=False)
            st.note("sources: " + "; ".join(
                "{} {:.3f}".format(r.source, r.plus_frac)
                for r in T.itertuples() if r.source in
                ("divergent", "intergenic", "intergenic_LTR",
                 "intergenic_nonLTR")))
            st.summary_md("results/strand_controls/step_summaries/47_s.md",
                          "47S strand balance by source", df_to_md(T, ".3f"))

        elif mode == "T":
            tt = pd.read_csv(os.path.join(DATA, "tss_in_te.tsv"), sep="\t",
                             header=None,
                             names=["c", "s", "e", "gene_id", "tc", "ts",
                                    "te", "tid", "t_sub", "tstrand", "fam",
                                    "order"])
            tt = tt.merge(meta[["chrom", "strand", "cls", "group"]],
                          left_on="gene_id", right_index=True, how="inner")
            tt["same"] = (tt.tstrand == tt.strand)
            tt["len"] = tt.te - tt.ts
            rows = []
            for lab, sub in [("all", tt), ("LTR", tt[tt.order == "LTR"]),
                             ("nonLTR", tt[tt.order != "LTR"])]:
                k, n = int(sub.same.sum()), len(sub)
                rows.append(dict(group=lab, n=n, same=k,
                                 same_frac=k / n,
                                 p_vs_half=binomtest(k, n, 0.5).pvalue))
            # family-level (LTR, n>=50)
            fams = tt[tt.order == "LTR"].t_sub.value_counts()
            for f in fams[fams >= 50].index[:25]:
                sub = tt[(tt.order == "LTR") & (tt.t_sub == f)]
                k, n = int(sub.same.sum()), len(sub)
                rows.append(dict(group="LTR:" + f, n=n, same=k,
                                 same_frac=k / n,
                                 p_vs_half=binomtest(k, n, 0.5).pvalue))
            for c in CLASSES:
                sub = tt[tt.cls == c]
                if len(sub) >= 50:
                    k, n = int(sub.same.sum()), len(sub)
                    rows.append(dict(group="class:" + c, n=n, same=k,
                                     same_frac=k / n,
                                     p_vs_half=binomtest(k, n, 0.5).pvalue))
            T = pd.DataFrame(rows)
            T["fdr"] = bh_fdr(T.p_vs_half.values)
            T.to_csv(safe_write_path(os.path.join(
                OUT, "T_te_orientation.tsv")), sep="\t", index=False)

            # chromosome-matched gene-strand permutation null
            rng = np.random.default_rng(42)
            obs_all = float(tt.same.mean())
            null = []
            by_chrom = {c: sub.strand.values
                        for c, sub in m.groupby("chrom")}
            for _ in range(1000):
                sim = []
                for c, sub in tt.groupby("chrom"):
                    arr = by_chrom.get(c)
                    if arr is None or len(arr) == 0:
                        continue
                    picks = rng.choice(arr, size=len(sub), replace=True)
                    sim.extend((sub.tstrand.values == picks).tolist())
                null.append(float(np.mean(sim)))
            null = np.array(null)
            p_emp = (1 + np.sum(null >= obs_all)) / (len(null) + 1)
            pd.DataFrame([dict(group="all_TE_TSS", observed=obs_all,
                               null_mean=float(null.mean()),
                               null_ci_low=float(np.percentile(null, 2.5)),
                               null_ci_high=float(np.percentile(null, 97.5)),
                               p_empirical=float(p_emp))]).to_csv(
                safe_write_path(os.path.join(
                    OUT, "T_permutation_null.tsv")), sep="\t", index=False)

            # logistic model (top families)
            top = list(fams[fams >= 50].index[:10])
            tt["fam_g"] = np.where(tt.t_sub.isin(top), tt.t_sub, "other")
            tt["same_i"] = tt.same.astype(int)
            tt["len_k"] = np.log1p(tt["len"].clip(lower=1))
            try:
                mod = logit("same_i ~ C(order) + C(fam_g) + len_k + C(cls)",
                            data=tt).fit(disp=0)
                coef = pd.DataFrame({"term": mod.params.index,
                                     "coef": mod.params.values,
                                     "pval": mod.pvalues.values})
                coef.to_csv(safe_write_path(os.path.join(
                    OUT, "T_logistic.tsv")), sep="\t", index=False)
            except Exception as e:
                st.note("logistic failed: {}".format(e))
            st.note("TE orientation: all={:.3f} null={:.3f} p={:.3f}; "
                    "LTR={:.3f}".format(
                        obs_all, null.mean(), p_emp,
                        T[T.group == "LTR"].same_frac.iloc[0]))
            st.summary_md("results/strand_controls/step_summaries/47_t.md",
                          "47T TE orientation controls", df_to_md(T, ".3f"))

        elif mode == "C":
            import pysam
            fa = pysam.FastaFile("/media/share/resource/reference/"
                                 "homo_sapiens/GRCh38/genome/GRCh38.CHR.fa")
            rows = []
            with gzip.open(os.path.join(DATA, "fantom",
                                        "hg38_cage_peaks.bed.gz"), "rt") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    s, e = int(f[1]), int(f[2])
                    try:
                        seq = fa.fetch(f[0], s, e).upper()
                    except Exception:
                        seq = ""
                    gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
                    rows.append((f[0], s, e, int(f[4]), f[5], e - s, gc))
            pk = pd.DataFrame(rows, columns=["chrom", "start", "end",
                                             "score", "strand", "len", "gc"])
            # PCG TSS proximity
            pcg = meta[meta.group == "pcg"]
            tss_by = {c: np.sort(sub.tss.values)
                      for c, sub in pcg.groupby("chrom")}
            near = np.zeros(len(pk), dtype=bool)
            for c, sub in pk.groupby("chrom"):
                arr = tss_by.get(c)
                if arr is None:
                    continue
                pos = sub.start.values
                i = np.searchsorted(arr, pos)
                d = np.full(len(pos), 1e9)
                for k in range(len(pos)):
                    cand = []
                    if i[k] < len(arr):
                        cand.append(arr[i[k]])
                    if i[k] > 0:
                        cand.append(arr[i[k] - 1])
                    if cand:
                        d[k] = min(abs(x - pos[k]) for x in cand)
                near[sub.index.values] = d <= 100
            pk["near_pcg_tss"] = near
            rows = []
            for lab, sub in [("all", pk),
                             ("near_PCG_TSS", pk[pk.near_pcg_tss]),
                             ("far_from_PCG_TSS", pk[~pk.near_pcg_tss])]:
                k, n = int((sub.strand == "+").sum()), len(sub)
                rows.append(dict(group=lab, n=n, plus_frac=k / n,
                                 p_vs_half=binomtest(k, n, 0.5).pvalue))
            for feat in ["score", "len", "gc"]:
                pk["q"] = pd.qcut(pk[feat], 5, labels=False, duplicates="drop")
                for q, sub in pk.groupby("q"):
                    k, n = int((sub.strand == "+").sum()), len(sub)
                    rows.append(dict(group="{}_q{}".format(feat, q), n=n,
                                     plus_frac=k / n,
                                     p_vs_half=binomtest(k, n, 0.5).pvalue))
            for c, sub in pk.groupby("chrom"):
                if len(sub) < 1000:
                    continue
                k, n = int((sub.strand == "+").sum()), len(sub)
                rows.append(dict(group="chrom:" + c, n=n,
                                 plus_frac=k / n,
                                 p_vs_half=binomtest(k, n, 0.5).pvalue))
            T = pd.DataFrame(rows)
            T["fdr"] = bh_fdr(T.p_vs_half.values)
            T.to_csv(safe_write_path(os.path.join(
                OUT, "C_cage_controls.tsv")), sep="\t", index=False)
            # logistic
            pk["plus"] = (pk.strand == "+").astype(int)
            pk["log_score"] = np.log1p(pk.score)
            pk["log_len"] = np.log1p(pk["len"])
            try:
                mod = logit("plus ~ log_score + log_len + gc + "
                            "near_pcg_tss + C(chrom)", data=pk).fit(disp=0)
                coef = pd.DataFrame({"term": mod.params.index,
                                     "coef": mod.params.values,
                                     "pval": mod.pvalues.values})
                coef.to_csv(safe_write_path(os.path.join(
                    OUT, "C_logistic.tsv")), sep="\t", index=False)
            except Exception as e:
                st.note("logistic failed: {}".format(e))
            st.note("CAGE plus: " + "; ".join(
                "{} {:.3f}".format(r.group, r.plus_frac)
                for r in T.itertuples() if r.group in
                ("all", "near_PCG_TSS", "far_from_PCG_TSS")))
            st.summary_md("results/strand_controls/step_summaries/47_c.md",
                          "47C CAGE plus-bias controls", df_to_md(T, ".3f"))
        else:
            raise SystemExit("mode S|T|C")


if __name__ == "__main__":
    main()
