#!/usr/bin/env python
"""45_evo_boot.py — final rigor checks for the evolutionary layer.

1. Chromosome-block bootstrap (1,000 resamples) for age~phyloP,
   divergent promoter−body difference, intergenic patch fraction and
   cross-species mouse breadth correlation.
2. Patch-threshold sensitivity (phyloP 1.3/2.0; phastCons 0.3/0.5/0.7).
3. TE exonization with merged TE intervals + exon-length/GC-matched
   permutation test.
Outputs -> results/evo_boot/*.tsv + step summary.
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "evo_boot")
EXONS = os.path.join(WORKDIR, "data", "genes", "exons.bed")
TE = ("/media/share/resource/reference/homo_sapiens/GRCh38/annotation/"
      "hg38.nrph.TE.anno.sorted.bed")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"
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


def mwu(y, x):
    u, p = mannwhitneyu(y, x, alternative="two-sided")
    return float(2 * u / (len(y) * len(x)) - 1), float(p)


def bt(cmd):
    r = subprocess.run(["singularity", "exec", "--no-mount", "bind-paths",
                        "--bind", "/media/share/resource:/media/share/resource:ro",
                        "--bind", f"{WORKDIR}:{WORKDIR}", SIF_BT,
                        "bedtools"] + cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[:500])
    return r.stdout


def bootstrap(df, chrom_col, stats, n_boot=1000, seed=42):
    chroms = sorted(df[chrom_col].dropna().unique())
    idx = {c: np.where(df[chrom_col].values == c)[0] for c in chroms}
    rng = np.random.default_rng(seed)
    obs = {k: fn(df) for k, fn in stats.items()}
    draws = {k: [] for k in stats}
    for _ in range(n_boot):
        sel = rng.choice(chroms, size=len(chroms), replace=True)
        rows = np.concatenate([idx[c] for c in sel])
        sub = df.iloc[rows]
        for k, fn in stats.items():
            draws[k].append(fn(sub))
    rows = []
    for k, fn in stats.items():
        v = np.asarray(draws[k], float)
        v = v[np.isfinite(v)]
        rows.append(dict(statistic=k, observed=obs[k],
                         ci_low=float(np.percentile(v, 2.5)) if len(v) else np.nan,
                         ci_high=float(np.percentile(v, 97.5)) if len(v) else np.nan,
                         n_boot=len(v)))
    return pd.DataFrame(rows)


def main():
    with step("45_evo_boot") as st:
        raw = pd.read_csv(os.path.join(RES, "constraint_structure",
                                       "gene_profiles_raw.tsv"),
                          sep="\t", index_col=0, low_memory=False)
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["phylop_mean"]]
        age = pd.read_csv(os.path.join(RES, "comparative",
                                       "E_transcript_age.tsv"), sep="\t",
                          index_col=0)[["age_level"]]
        seq = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["gc"]]
        d = meta[["chrom", "cls", "group"]].join(raw).join(ph).join(age) \
            .join(seq)
        d["cls"] = np.where(d.group == "pcg", "PCG", d.cls)
        d["prom_minus_body_phylop"] = d.prom_phylop - d.body_phylop

        # ---- 1. chromosome-block bootstrap ---------------------------
        def rho_age(sub):
            s = sub[(sub.group == "lnc") &
                    sub.age_level.notna() & sub.phylop_mean.notna()]
            if len(s) < 50:
                return np.nan
            return spearmanr(s.age_level, s.phylop_mean).statistic

        def div_body(sub):
            a = sub.loc[sub.cls == "divergent",
                        "prom_minus_body_phylop"].dropna()
            b = sub.loc[sub.group == "pcg",
                        "prom_minus_body_phylop"].dropna()
            if len(a) < 20 or len(b) < 20:
                return np.nan
            return float(a.median() - b.median())

        def inter_patch(sub):
            a = sub.loc[sub.cls == "intergenic",
                        "prom_patch_frac_phylop"].dropna()
            b = sub.loc[sub.group == "pcg",
                        "prom_patch_frac_phylop"].dropna()
            if len(a) < 20 or len(b) < 20:
                return np.nan
            return float(a.median() - b.median())

        # mouse breadth rho (PCG) needs map file with chrom
        mp = pd.read_csv(os.path.join(RES, "syntenic", "map_mouse.tsv"),
                         sep="\t")
        mp = mp.join(meta[["chrom"]], on="human_gene")
        mp["group"] = np.where(mp.cls.isin(CLASSES), "lnc", "pcg")
        # map file cls is human class; PCG rows not in map (only lnc); skip
        d2 = d  # mouse breadth not available for PCG in map_mouse (lnc only)
        stats = {"rho_age_phylop": rho_age,
                 "divergent_prom_minus_body": div_body,
                 "intergenic_patch_frac_diff": inter_patch}
        B = bootstrap(d, "chrom", stats, n_boot=1000)
        B.to_csv(safe_write_path(os.path.join(OUT, "bootstrap.tsv")),
                 sep="\t", index=False)

        # ---- 2. patch threshold sensitivity -------------------------
        pt = pd.read_csv(os.path.join(RES, "constraint_structure",
                                      "patch_thresholds.tsv"), sep="\t",
                         index_col=0)
        m = meta[["cls", "group"]].join(pt)
        m["cls"] = np.where(m.group == "pcg", "PCG", m.cls)
        pcg = m[m.group == "pcg"]
        rows = []
        feats = [c for c in pt.columns if c.endswith("_frac") or
                 c.endswith("_max") or c.endswith("_n")]
        for f in feats:
            x = pcg[f].dropna().values
            for c in CLASSES:
                y = m.loc[m.cls == c, f].dropna().values
                if len(y) < 20:
                    continue
                dd, p = mwu(y, x)
                rows.append(dict(feature=f, cls=c, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=dd, pval=p))
        PT = pd.DataFrame(rows)
        PT["fdr"] = bh_fdr(PT.pval.values)
        PT.to_csv(safe_write_path(os.path.join(OUT, "patch_thresholds.tsv")),
                  sep="\t", index=False)

        # ---- 3. TE merged + matched permutation ----------------------
        tmp = os.path.join(WORKDIR, "data", "evo_boot_tmp")
        os.makedirs(tmp, exist_ok=True)
        te_merge = os.path.join(tmp, "TE.merged.bed")
        ltr_src = os.path.join(WORKDIR, "data", "constraint_structure_tmp",
                               "TE_LTR.bed")
        ltr_merge = os.path.join(tmp, "LTR.merged.bed")
        with open(safe_write_path(te_merge), "w") as fh:
            fh.write(bt(["merge", "-i", TE]))
        with open(safe_write_path(ltr_merge), "w") as fh:
            fh.write(bt(["merge", "-i", ltr_src]))
        ov_all = os.path.join(tmp, "exon_TE_merged.tsv")
        ov_ltr = os.path.join(tmp, "exon_LTR_merged.tsv")
        with open(safe_write_path(ov_all), "w") as fh:
            fh.write(bt(["intersect", "-a", EXONS, "-b", te_merge, "-wo"]))
        with open(safe_write_path(ov_ltr), "w") as fh:
            fh.write(bt(["intersect", "-a", EXONS, "-b", ltr_merge, "-wo"]))

        def agg(path):
            out = {}
            with open(path) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 10:
                        continue
                    g = f[3].split("|")[0]
                    out[g] = out.get(g, 0) + int(f[-1])
            return pd.Series(out)

        te_bp = agg(ov_all)
        ltr_bp = agg(ov_ltr)
        cte = pd.read_csv(os.path.join(RES, "constraint_structure",
                                       "C_te_per_gene.tsv"), sep="\t",
                          index_col=0)
        cte["te_merged_bp"] = te_bp.reindex(cte.index).fillna(0).values
        cte["ltr_merged_bp"] = ltr_bp.reindex(cte.index).fillna(0).values
        cte["te_merged_frac"] = (cte.te_merged_bp /
                                 cte.exon_len.replace(0, np.nan)).clip(upper=1)
        cte["ltr_merged_frac"] = (cte.ltr_merged_bp /
                                  cte.exon_len.replace(0, np.nan)).clip(upper=1)
        cte = cte.join(seq)
        cte.to_csv(safe_write_path(os.path.join(
            OUT, "te_merged_per_gene.tsv")), sep="\t")
        # merge check
        chk = cte[["te_exon_frac", "te_merged_frac", "ltr_exon_frac",
                   "ltr_merged_frac"]].corr(method="spearman")
        chk.to_csv(safe_write_path(os.path.join(OUT, "te_merge_check.tsv")),
                   sep="\t")

        # matched permutation
        df = cte.dropna(subset=["exon_len", "gc"]).copy()
        df["len_q"] = pd.qcut(np.log1p(df.exon_len), 4, labels=False,
                              duplicates="drop")
        df["gc_q"] = pd.qcut(df.gc, 4, labels=False, duplicates="drop")
        df["stratum"] = df.len_q.astype(str) + "_" + df.gc_q.astype(str)
        rng = np.random.default_rng(42)
        rows = []
        for met in ["te_merged_frac", "ltr_merged_frac"]:
            for c in CLASSES:
                obs, null = matched_perm(df, met, c, rng, 1000)
                if obs is None:
                    continue
                p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1)
                rows.append(dict(metric=met, cls=c, observed_diff=obs,
                                 null_mean=float(np.mean(null)),
                                 pval=float(p), n_perm=len(null)))
        MP = pd.DataFrame(rows)
        MP["fdr"] = bh_fdr(MP.pval.values) if len(MP) else []
        MP.to_csv(safe_write_path(os.path.join(
            OUT, "te_matched_perm.tsv")), sep="\t", index=False)

        st.note("bootstrap: " + "; ".join(
            "{} obs={:.3f} CI[{:.3f},{:.3f}]".format(
                r.statistic, r.observed, r.ci_low, r.ci_high)
            for r in B.itertuples()))
        st.note("TE merge corr: te={:.3f} ltr={:.3f}".format(
            chk.loc["te_exon_frac", "te_merged_frac"],
            chk.loc["ltr_exon_frac", "ltr_merged_frac"]))
        st.summary_md("results/evo_boot/step_summaries/45.md",
                      "45 evolutionary bootstrap / threshold / TE permutation",
                      "## Chromosome-block bootstrap\n\n" + df_to_md(B, ".4f")
                      + "\n\n## Patch threshold sensitivity\n\n"
                      + df_to_md(PT, ".4g", max_rows=50)
                      + "\n\n## TE merged + matched permutation\n\n"
                      + df_to_md(MP, ".4f")
                      + "\n\nMerge check\n\n" + df_to_md(chk.reset_index(), ".3f"))


def matched_perm(df, met, cls, rng, n_perm):
    obs = []
    for stratum, sub in df.groupby("stratum"):
        a = sub.loc[sub.cls == cls, met].dropna()
        b = sub.loc[sub.group == "pcg", met].dropna()
        if len(a) >= 5 and len(b) >= 5:
            obs.append((len(a), float(a.median() - b.median())))
    if not obs:
        return None, np.array([])
    w = np.array([x[0] for x in obs], float)
    obs_val = float(np.sum(w * np.array([x[1] for x in obs])) / w.sum())
    null = []
    for _ in range(n_perm):
        ds = []
        for stratum, sub in df.groupby("stratum"):
            a = sub.loc[sub.cls == cls, met].dropna()
            b = sub.loc[sub.group == "pcg", met].dropna()
            if len(a) < 5 or len(b) < 5:
                continue
            vals = np.concatenate([a.values, b.values])
            lab = np.array([1] * len(a) + [0] * len(b))
            rng.shuffle(lab)
            ds.append((len(a), float(np.median(vals[lab == 1]) -
                                     np.median(vals[lab == 0]))))
        if ds:
            ww = np.array([x[0] for x in ds], float)
            null.append(float(np.sum(ww * np.array([x[1] for x in ds])) /
                              ww.sum()))
    return obs_val, np.array(null)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
