#!/usr/bin/env python
"""80_reviewer_checks.py — targeted reanalyses requested by the panel.

Four checks that can be run from existing tables.

  A. Is the activity hierarchy independent of CAGE TSS refinement? Compare
     the CAGE-active fraction per class with the CAGE-refinement rate, and
     with an independent RNA-seq detection fraction from GTEx.
  B. How much do the two promoter windows of a divergent pair overlap, and
     does the 62% versus 98% CAGE asymmetry survive among pairs whose windows
     do not overlap?
  C. Is the LTR-TSS enrichment above the local LTR density? Shuffle each
     intergenic TSS within +/- 50 kb and compare the point-in-LTR rate.
  D. Is the GBM LTR de-repression Fisher test significant after correction
     across the 18 cancers?

Outputs
  results/revision_checks/A_activity_circularity.tsv
  results/revision_checks/B_window_overlap.tsv
  results/revision_checks/C_ltr_local_null.tsv
  results/revision_checks/D_gbm_multiplicity.tsv
  results/step_summaries/80_reviewer_checks.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "revision_checks")
TE_BED = ("/media/share/resource/reference/homo_sapiens/GRCh38/annotation/"
          "hg38.nrph.TE.anno.sorted.bed")
CLASSES = ["PCG", "divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic", "pseudogene", "transcribed_pseudogene"]
NSHUF = 20
WIN = 50_000
SEED = 42


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def load_ltr():
    out = {}
    with open(TE_BED) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[7] != "LTR":
                continue
            c = f[0].replace("chr", "")
            out.setdefault(c, []).append((int(f[1]), int(f[2])))
    merged = {}
    for c, v in out.items():
        v.sort()
        s, e = [], []
        for a, b in v:
            if s and a <= e[-1]:
                e[-1] = max(e[-1], b)
            else:
                s.append(a); e.append(b)
        merged[c] = (np.array(s), np.array(e))
    return merged


def in_ltr(iv, chrom, pos):
    t = iv.get(chrom)
    if t is None:
        return False
    s, e = t
    j = np.searchsorted(s, pos, "right") - 1
    return bool(j >= 0 and e[j] >= pos)


def main():
    with step("80_reviewer_checks") as st:
        os.makedirs(OUT, exist_ok=True)
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["base"] = m.index.astype(str).str.split(".").str[0]
        m["tss"] = pd.to_numeric(m.tss, errors="coerce")
        m["chrom_c"] = m.chrom.astype(str).str.replace("chr", "", regex=False)

        # ---------------- A. activity versus refinement --------------------
        ca = pd.read_csv(os.path.join(RES, "cage_activity_by_class.tsv"),
                         sep="\t").set_index("cls")
        gt = pd.read_csv(os.path.join(
            RES, "tissue_disease", "tables", "gtex_gene_tau.tsv"), sep="\t")
        gt = gt.rename(columns={"Unnamed: 0": "base"})
        gt["base"] = gt.base.astype(str).str.split(".").str[0]
        gt = gt.drop_duplicates("base").set_index("base")
        rows = []
        for c in CLASSES:
            sub = m[m.cls == c]
            if not len(sub):
                continue
            refined = float((sub.cage_refined == 1).mean())
            active = float(ca.loc[c, "active_frac"]) if c in ca.index else np.nan
            g = gt.reindex(sub.base).max_tissue_mean.dropna()
            expr = float((g > 1).mean()) if len(g) else np.nan
            rows.append(dict(cls=c, n=len(sub), cage_refined=refined,
                             cage_active=active, gtex_expressed=expr,
                             n_gtex=len(g)))
        A = pd.DataFrame(rows)
        A.to_csv(safe_write_path(os.path.join(
            OUT, "A_activity_circularity.tsv")), sep="\t", index=False)

        # ---------------- B. window overlap and CAGE subset -----------------
        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        tss = m.set_index("base").tss
        dp["ltss"] = dp.lnc_gid.astype(str).str.split(".").str[0].map(tss)
        dp["ptss"] = dp.pcg_gid.astype(str).str.split(".").str[0].map(tss)
        dp = dp.dropna(subset=["ltss", "ptss"])
        # windows are lnc -500/+200 on the lncRNA strand; use genomic span
        ov = np.clip(np.minimum(dp.ltss + 200, dp.ptss + 200)
                     - np.maximum(dp.ltss - 500, dp.ptss - 500), 0, None)
        dp["overlap_bp"] = ov
        dp["jaccard"] = ov / (700 * 2 - ov)
        cage = pd.read_csv(os.path.join(RES, "cage_activity.tsv"),
                           sep="\t", low_memory=False)
        cage["base"] = cage.gene_id.astype(str).str.split(".").str[0]
        cage = cage.drop_duplicates("base").set_index("base")
        dp["lnc_active"] = dp.lnc_gid.astype(str).str.split(".").str[0].map(
            cage.active)
        dp["pcg_active"] = dp.pcg_gid.astype(str).str.split(".").str[0].map(
            cage.active)
        dp = dp.dropna(subset=["lnc_active", "pcg_active"])
        rows = []
        for lab, sub in (("all pairs", dp),
                         ("windows disjoint (>700 bp)", dp[dp.overlap_bp == 0]),
                         ("windows overlap", dp[dp.overlap_bp > 0])):
            if len(sub) < 20:
                continue
            rows.append(dict(subset=lab, n=len(sub),
                             median_overlap_bp=float(sub.overlap_bp.median()),
                             median_jaccard=float(sub.jaccard.median()),
                             frac_lnc_active=float(sub.lnc_active.mean()),
                             frac_pcg_active=float(sub.pcg_active.mean())))
        B = pd.DataFrame(rows)
        B.to_csv(safe_write_path(os.path.join(
            OUT, "B_window_overlap.tsv")), sep="\t", index=False)
        st.note("window overlap: median {} bp, Jaccard {:.3f}; disjoint pairs "
                "CAGE {:.1%} v {:.1%}".format(
                    int(dp.overlap_bp.median()), float(dp.jaccard.median()),
                    float(B.loc[B.subset.str.startswith("windows disjoint"),
                                "frac_lnc_active"].iloc[0]),
                    float(B.loc[B.subset.str.startswith("windows disjoint"),
                                "frac_pcg_active"].iloc[0])))

        # ---------------- C. LTR local shuffle ------------------------------
        ltr = load_ltr()
        ig = m[(m.cls == "intergenic") & m.tss.notna()].copy()
        rng = np.random.default_rng(SEED)
        obs = []
        null = []
        for r in ig.itertuples():
            if r.chrom_c not in ltr:
                continue
            pos = int(r.tss)
            obs.append(in_ltr(ltr, r.chrom_c, pos))
            qi = rng.integers(-WIN, WIN + 1, NSHUF)
            null.append([in_ltr(ltr, r.chrom_c, max(0, pos + int(q)))
                         for q in qi])
        obs = np.array(obs, dtype=bool)
        null = np.array(null, dtype=bool)
        obsf = float(obs.mean())
        nullf = float(null.mean())
        p = binomtest(int(obs.sum()), len(obs), nullf,
                      alternative="greater").pvalue
        C = pd.DataFrame([dict(n=len(obs), observed_frac=obsf,
                               local_null_frac=nullf,
                               ratio=(obsf / nullf if nullf > 0 else np.nan),
                               pval=p)])
        C.to_csv(safe_write_path(os.path.join(
            OUT, "C_ltr_local_null.tsv")), sep="\t", index=False)

        # ---------------- D. GBM multiplicity -------------------------------
        d = pd.read_csv(os.path.join(
            RES, "tissue_disease", "tables", "tcga_ltr_derepression.tsv"),
            sep="\t")
        d["fdr"] = bh_fdr(d.fisher_p.values)
        d["bonferroni"] = np.minimum(d.fisher_p * len(d), 1)
        D = d[["cancer", "diff", "fisher_p", "fdr", "bonferroni"]].sort_values(
            "fisher_p")
        D.to_csv(safe_write_path(os.path.join(
            OUT, "D_gbm_multiplicity.tsv")), sep="\t", index=False)

        body = "\n\n".join([
            "## A. activity versus TSS refinement\n\n" + df_to_md(A, ".3f"),
            "## B. window overlap and CAGE asymmetry\n\n" + df_to_md(B, ".3f"),
            "## C. LTR local shuffle\n\n" + df_to_md(C, ".4f"),
            "## D. LTR de-repression across 18 cancers\n\n" + df_to_md(D, ".3g"),
        ])
        st.summary_md(os.path.join(
            RES, "step_summaries", "80_reviewer_checks.md"),
            "80 reviewer-requested reanalyses", body)
        st.note("LTR intergenic observed {:.1%} v local null {:.1%} "
                "(ratio {:.2f}, P={:.2g})".format(
                    obsf, nullf, obsf / nullf if nullf else np.nan, p))
        st.note("GBM Fisher P={:.2g}, FDR={:.2g}, Bonferroni={:.2g}".format(
            float(d.loc[d.cancer == "TCGA-GBM", "fisher_p"].iloc[0]),
            float(d.loc[d.cancer == "TCGA-GBM", "fdr"].iloc[0]),
            float(d.loc[d.cancer == "TCGA-GBM", "bonferroni"].iloc[0])))


if __name__ == "__main__":
    main()
