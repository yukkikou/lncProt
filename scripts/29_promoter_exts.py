#!/usr/bin/env python
"""29_promoter_exts.py — promoter-level extension analyses (analysis-first).

Modes:
  A -> promoter shape from local FANTOM5 peaks: per promoter the number of
       CAGE peaks, top-peak width (highest score), summit spread; class tests
       and CpG-syntax association (sharp/broad promoter biology).
  B -> TE-donor age x phyloP: per-gene TSS-donating TE subfamily (06 output)
       binned by curated RepName age (12b rules); age vs phylop_mean /
       frac_accel_2 within TE-derived lnc promoters (LTR subset headline).
  C -> class-vs-biotype variance: rank-ANOVA eta^2 for the classic 4-way
       coordinate taxonomy (divergent folded into intergenic) vs the birth
       5-way split, on conservation/activity/composition features.
Outputs -> results/promoter_ext/*.tsv + step summary.
"""
import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "promoter_ext")
DATA = os.path.join(WORKDIR, "data")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]

AGE_RULES = [
    ("ancient", ("L2", "L3", "MIR", "CR1", "MER2", "tRNA")),
    ("primate_young", ("AluY", "L1PA1", "L1PA2", "L1HS", "LTR12", "HERVK",
                       "SVA", "LTR13", "ERVK", "LTR5", "MER51")),
    ("primate_old", ("AluJ", "AluS", "L1PA3", "L1PA4", "L1PA5", "L1PA6",
                     "L1PA7", "L1PB", "THE1", "MST", "Charlie", "MLT1", "MLT2",
                     "ERV1", "ERVL", "LTR7", "LTR2", "MER4", "MER61")),
]
AGE_ORD = {"ancient": 0, "primate_old": 1, "primate_young": 2}


def age_of(sub):
    s = str(sub)
    for age, prefixes in AGE_RULES:
        if any(s.startswith(p) for p in prefixes):
            return age
    return "unclassified"


def sing(*args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}"] + list(args)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(" ".join(cmd) + "\n" + r.stderr[:1000])
    return r.stdout


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def mwu(y, x):
    u, p = mannwhitneyu(y, x, alternative="two-sided")
    return float(2 * u / (len(y) * len(x)) - 1), float(p)


def eta2_rank(groups):
    """One-way ANOVA eta^2 on rank-transformed values (list of arrays)."""
    groups = [np.asarray(g, float) for g in groups if len(g) >= 5]
    allv = np.concatenate(groups)
    ranks = pd.Series(allv).rank().values
    gm = ranks.mean()
    ssb, k = 0.0, 0
    for g in groups:
        r = ranks[k:k + len(g)]
        ssb += len(g) * (r.mean() - gm) ** 2
        k += len(g)
    sst = ((ranks - gm) ** 2).sum()
    return ssb / sst if sst > 0 else np.nan


# ------------------------------------------------------------------ mode A --
def mode_a():
    with step("29_ext_A") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        # flatten peaks once
        pk = os.path.join(DATA, "cage_peaks_flat.bed")
        if not os.path.exists(pk):
            src = gzip.open(
                os.path.join(DATA, "fantom", "hg38_cage_peaks.bed.gz"), "rt")
            with src as fi, open(safe_write_path(pk), "w") as fo:
                for line in fi:
                    fo.write(line)
        ovf = os.path.join(DATA, "shape_ov.tsv")
        cmd = ["singularity", "exec", "--no-mount", "bind-paths",
               "--bind", "/media/share/resource:/media/share/resource:ro",
               "--bind", f"{WORKDIR}:{WORKDIR}", SIF_BT, "bedtools", "intersect",
               "-a", os.path.join(RES, "promoters_feat.bed"), "-b", pk, "-wo"]
        with open(safe_write_path(ovf), "w") as fh:
            subprocess.run(cmd, stdout=fh, check=True)
        cols = ["c", "s", "e", "gene_id", "x", "strand",
                "pc", "ps", "pe", "pid", "score", "pstrand",
                "sum_s", "sum_e", "rgb", "ov"]
        ov = pd.read_csv(ovf, sep="\t", header=None, names=cols,
                         usecols=["gene_id", "ps", "pe", "score", "sum_s"])
        ov["width"] = ov["pe"] - ov["ps"]
        agg = ov.groupby("gene_id").agg(
            n_peaks=("width", "size"),
            top_width=("width", lambda w: w.max()),
            med_width=("width", "median"),
            best_score=("score", "max"),
            spread=("sum_s", lambda s: s.max() - s.min()))
        m = meta.join(agg, how="left")
        m.loc[m.index.difference(agg.index),
              ["n_peaks", "top_width", "med_width", "best_score",
               "spread"]] = 0
        m[["n_peaks", "top_width", "med_width", "spread"]].to_csv(
            safe_write_path(os.path.join(OUT, "A_shape_per_gene.tsv")),
            sep="\t")
        seq = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["cpg_oe"]]
        m = m.join(seq)

        rows = []
        # shape is only defined for CAGE-supported promoters: run the battery
        # twice — all promoters (peak absence inflates zeros) and the
        # peak-bearing subset (the shape question proper)
        for subset, frame in [("all", m), ("peak_bearing", m[m.n_peaks > 0])]:
            pcg = frame[frame.group == "pcg"]
            for feat in ["top_width", "n_peaks", "spread"]:
                x = pcg[feat].dropna().values
                if len(x) < 20:
                    continue
                for cls in CLASSES + ["pseudogene",
                                      "transcribed_pseudogene"]:
                    y = frame.loc[frame.cls == cls, feat].dropna().values
                    if len(y) < 20:
                        continue
                    d, p = mwu(y, x)
                    rows.append(dict(subset=subset, feature=feat, cls=cls,
                                     n=len(y),
                                     median_class=float(np.median(y)),
                                     median_PCG=float(np.median(x)),
                                     cliffs_delta=d, pval=p))
        A = pd.DataFrame(rows)
        A["fdr"] = bh_fdr(A.pval.values)
        A.to_csv(safe_write_path(os.path.join(OUT, "A_shape_by_class.tsv")),
                 sep="\t", index=False)
        ok = m[["cpg_oe", "top_width", "n_peaks"]].dropna()
        r1, p1 = spearmanr(ok.cpg_oe, ok.top_width)
        st.note("A: top_width x cpg_oe rho={:.3f} p={:.1e}; PCG median "
                "width={} n_peaks={}; divergent {} / intergenic {}".format(
                    r1, p1,
                    pcg.top_width.median(), pcg.n_peaks.median(),
                    m[m.cls == "divergent"].top_width.median(),
                    m[m.cls == "intergenic"].top_width.median()))
        st.summary_md("results/promoter_ext/step_summaries/29_a.md",
                      "29A promoter shape (FANTOM peaks)",
                      "## Class vs PCG (promoter shape)\n\n"
                      + df_to_md(A, ".3f", max_rows=30)
                      + "\n\nsharp/broad: top_width = width of the "
                        "highest-scoring CAGE peak in the -500/+200 window; "
                      "n_peaks = number of CAGE peaks; spread = summit span.")


# ------------------------------------------------------------------ mode B --
def mode_b():
    with step("29_ext_B") as st:
        tt = pd.read_csv(os.path.join(DATA, "tss_in_te.tsv"), sep="\t",
                         header=None,
                         names=["c", "s", "e", "gene_id", "tc", "ts", "te",
                                "tid", "t_sub", "tstrand", "order", "fam"],
                         usecols=["gene_id", "t_sub", "order"])
        tt = tt.drop_duplicates("gene_id")
        tt["age"] = tt.t_sub.map(age_of)
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        B = tt.join(ph, on="gene_id", how="inner")
        lnc = B[(B.group == "lnc") & B.age.isin(AGE_ORD)]
        lnc = lnc.assign(age_ord=lnc.age.map(AGE_ORD))

        rows = []
        for label, frame in [("TE_TSS_lnc_all", lnc),
                             ("LTR_TSS_lnc", lnc[lnc.tss_in_ltr == 1]),
                             ("intergenic_TE_TSS", lnc[lnc.cls ==
                                                        "intergenic"])]:
            for feat in ["phylop_mean", "frac_accel_2"]:
                ok = frame[["age_ord", feat]].dropna()
                r, p = spearmanr(ok.age_ord, ok[feat])
                rows.append(dict(scope=label, feature=feat, n=len(ok),
                                 spearman_rho=r, pval=p))
            for feat in ["phylop_mean"]:
                yo = ok = frame[(frame.age == "primate_young")
                                & frame[feat].notna()][feat].values
                an = frame[(frame.age == "ancient")
                           & frame[feat].notna()][feat].values
                if len(yo) >= 20 and len(an) >= 20:
                    d, p = mwu(yo, an)
                    rows.append(dict(scope=label + "_young_v_ancient",
                                     feature=feat, n=len(yo) + len(an),
                                     cliffs_delta=d, pval=p))
        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(OUT, "B_te_age_phylop.tsv")),
                 sep="\t", index=False)
        med = lnc.groupby("age")[["phylop_mean", "frac_accel_2"]].median()
        med.to_csv(safe_write_path(os.path.join(
            OUT, "B_age_medians.tsv")), sep="\t")
        st.note("B: age medians phylop " + "; ".join(
            "{}={:.3f}".format(a, v) for a, v in
            med.phylop_mean.items())
            + "; LTR subset rho(age, phylop)=" + "{:.3f}".format(
                T[(T.scope == "LTR_TSS_lnc") &
                  (T.feature == "phylop_mean")].spearman_rho.iloc[0]))
        st.summary_md("results/promoter_ext/step_summaries/29_b.md",
                      "29B TE-donor age x phyloP",
                      "## Age gradient tests\n\n" + df_to_md(T, ".4g")
                      + "\n\n## Median by donor age (TE-TSS lnc)\n\n"
                      + df_to_md(med.reset_index(), ".4f"))


# ------------------------------------------------------------------ mode C --
def mode_c():
    with step("29_ext_C") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        act = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["tau", "active"]]
        seq = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["cpg_oe"]]
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["te_cov"]]
        m = meta.join(ph[["phylop_mean", "phast_mean"]]).join(act).join(
            seq).join(te)
        lnc = m[m.group == "lnc"].copy()
        # classic 4-way taxonomy: divergent folded back into intergenic
        lnc["cls4"] = lnc.cls.replace({"divergent": "intergenic"})

        feats = ["phylop_mean", "phast_mean", "tau", "cpg_oe", "te_cov"]
        rows = []
        for f in feats:
            d = lnc[[f, "cls", "cls4"]].dropna()
            g5 = [d[d.cls == c][f].values for c in CLASSES]
            g4 = [d[d.cls4 == c][f].values
                  for c in ["antisense", "intronic", "sense_overlapping",
                            "intergenic"]]
            e5, e4 = eta2_rank(g5), eta2_rank(g4)
            rows.append(dict(feature=f, n=len(d), eta2_biotype4=e4,
                             eta2_birthclass5=e5,
                             delta=e5 - e4,
                             rel_gain=(e5 - e4) / e4 if e4 > 0 else np.nan))
        C = pd.DataFrame(rows)
        C.to_csv(safe_write_path(os.path.join(
            OUT, "C_class_vs_biotype.tsv")), sep="\t", index=False)
        st.note("C eta^2 gain from divergent split: " + "; ".join(
            "{}: {:+.4f} ({:+.0%})".format(r.feature, r.delta, r.rel_gain
                                           if r.rel_gain == r.rel_gain else 0)
            for r in C.itertuples()))
        st.summary_md("results/promoter_ext/step_summaries/29_c.md",
                      "29C class-vs-biotype variance (rank-ANOVA eta^2)",
                      "GENCODE v46 labels all 19,258 lncRNAs as a single "
                      "'lncRNA' biotype, so the meaningful comparison is the "
                      "classic 4-way coordinate taxonomy (divergent folded "
                      "into intergenic) vs the birth 5-way split.\n\n"
                      + df_to_md(C, ".4f"))


if __name__ == "__main__":
    mode = (sys.argv[1] if len(sys.argv) > 1 else "").upper()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    if mode == "A":
        mode_a()
    elif mode == "B":
        mode_b()
    elif mode == "C":
        mode_c()
    else:
        raise SystemExit("usage: 29_promoter_exts.py A|B|C")
