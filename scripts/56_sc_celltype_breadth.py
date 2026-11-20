#!/usr/bin/env python
"""56_sc_celltype_breadth.py — is class ordering constant across cell types?

The GTEx analysis asks how widely each class is expressed across tissues. The
single-cell compendium allows the same question at cell-type resolution, on a
quantity that does not move with sequencing depth, because each class is ranked
against the others inside the same cell type. A class whose rank is fixed
across cell types has the same relative preference everywhere. A class whose
rank moves is read out differently by different lineages.

Three metrics per sample and class
  rank_sd    standard deviation of the within-cell-type rank (depth free)
  tau_ratio  tissue-specificity index of the class-to-PCG detection ratio. Both
             numerator and denominator move with the cell's depth, so the ratio
             cancels most of it. This asks how far the lncRNA share of the
             coding output shifts between lineages (primary)
  tau        tissue-specificity index of the raw detection rate (reported as a
             contrast, it mostly tracks per-cell depth)

Outputs
  results/scbasecount/celltype_breadth_by_sample.tsv
  results/scbasecount/celltype_breadth_tests.tsv
  results/scbasecount/step_summaries/56.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, spearmanr, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

OUT = os.path.join(WORKDIR, "results", "scbasecount")
RES = os.path.join(WORKDIR, "results")
CLS = ["PCG", "sense_overlapping", "divergent", "antisense", "intergenic",
       "intronic"]
MIN_CELLS = 100          # a cell type must be well sampled to carry a rank
MIN_TYPES = 4            # a sample needs a few cell types to have a spread


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def tau(x):
    """Tissue-specificity index. 0 = flat across cell types, 1 = one type only."""
    x = np.asarray(x, float)
    mx = np.nanmax(x)
    if not np.isfinite(mx) or mx <= 0:
        return np.nan
    return float(np.nanmean(1 - x / mx))


def per_sample(S):
    rows = []
    S = S[S.n_cells >= MIN_CELLS]
    det_cols = ["{}_det".format(c) for c in CLS]
    for srx, sub in S.groupby("srx"):
        sub = sub.dropna(subset=det_cols)
        if len(sub) < MIN_TYPES:
            continue
        M = sub[det_cols].values                      # cell types x classes
        # Rank inside each cell type. Depth cancels here because every class in
        # a cell type is divided by the same cell count and the same per-cell
        # gene count, so a rank is invariant to both.
        R = np.apply_along_axis(
            lambda v: pd.Series(v).rank(ascending=False).values, 1, M)
        row = dict(srx=srx, group=sub.group.iloc[0], arm=sub.arm.iloc[0],
                   n_types=len(sub), cells=int(sub.n_cells.sum()))
        pcg = M[:, 0]
        for j, c in enumerate(CLS):
            row["{}_rank_sd".format(c)] = float(np.std(R[:, j], ddof=1))
            row["{}_rank_range".format(c)] = float(np.ptp(R[:, j]))
            row["{}_tau".format(c)] = tau(M[:, j])
            row["{}_top_frac".format(c)] = float((R[:, j] <= 2).mean())
            # Class share of coding output. Dividing by the PCG rate in the
            # same cell type removes the cell type's depth, which is what makes
            # a raw tau uninterpretable here.
            with np.errstate(invalid="ignore", divide="ignore"):
                share = M[:, j] / pcg
            row["{}_share_med".format(c)] = float(np.nanmedian(share))
            row["{}_tau_ratio".format(c)] = tau(share)
        rows.append(row)
    return pd.DataFrame(rows)


def gtex_breadth():
    """Class-level tissue specificity from the GTEx table, for comparison.

    frac_expr is the fraction of genes in a class detected in a tissue, so the
    tau over tissues describes how evenly a class is represented in the body.
    """
    g = pd.read_csv(os.path.join(RES, "tissue_disease", "tables",
                                 "gtex_class_by_tissue.tsv"), sep="\t")
    out = {}
    for c in CLS:
        col = "{}_frac_expr".format(c)
        if col in g:
            out[c] = tau(g[col].values)
    return out


def tests(B):
    rows = []
    metrics = ["rank_sd", "tau_ratio", "tau"]
    for m in metrics:
        cols = ["{}_{}".format(c, m) for c in CLS]
        D = B[cols].dropna()
        if len(D) < 10:
            continue
        stat, p = friedmanchisquare(*[D[c].values for c in cols])
        rows.append(dict(metric=m, test="Friedman across six entities",
                         n=len(D), statistic=float(stat), pval=float(p),
                         comparison="all classes", median_diff=np.nan))
        ref = "divergent_{}".format(m)
        for c in CLS:
            if c == "divergent":
                continue
            col = "{}_{}".format(c, m)
            d = D[col] - D[ref]
            if not np.isfinite(d).any() or np.allclose(d, 0):
                continue
            w, pw = wilcoxon(d, zero_method="wilcox")
            rows.append(dict(metric=m, test="Wilcoxon signed-rank",
                             n=len(d),
                             statistic=float(w), pval=float(pw),
                             comparison="{} vs divergent".format(c),
                             median_diff=float(np.median(d))))
    T = pd.DataFrame(rows)
    T["fdr"] = np.nan
    for m, sub in T.groupby("metric"):
        T.loc[sub.index, "fdr"] = bh_fdr(sub.pval.values)
    return T


def depth_diagnostic(S):
    """How much of each metric is just per-cell depth, across cell types.

    A correlation near one means the quantity is a depth readout and not a
    lineage property. It is reported so that tau and tau_ratio are not read as
    equivalent evidence.
    """
    S = S[S.n_cells >= MIN_CELLS]
    rows = []
    g = S.groupby("cell_type")
    depth = g.mean_genes_per_cell.mean()
    for c in CLS[1:]:
        det = g["{}_det".format(c)].mean()
        share = det / g["PCG_det"].mean()
        r1 = spearmanr(det.reindex(depth.index), depth)
        r2 = spearmanr(share.reindex(depth.index), depth)
        rows.append(dict(klass=c, n_types=len(depth),
                         rho_det_vs_depth=float(r1[0]),
                         rho_share_vs_depth=float(r2[0])))
    return pd.DataFrame(rows)


def main():
    with step("56_sc_celltype_breadth") as st:
        S = pd.read_csv(os.path.join(OUT, "sample_celltype_detection.tsv"),
                        sep="\t")
        B = per_sample(S)
        B.to_csv(safe_write_path(os.path.join(
            OUT, "celltype_breadth_by_sample.tsv")), sep="\t", index=False)
        T = tests(B)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "celltype_breadth_tests.tsv")), sep="\t", index=False)

        summ = pd.DataFrame({
            "class": CLS,
            "rank_sd_median": [B["{}_rank_sd".format(c)].median()
                               for c in CLS],
            "rank_span_median": [B["{}_rank_range".format(c)].median()
                                 for c in CLS],
            "share_median": [B["{}_share_med".format(c)].median()
                             for c in CLS],
            "tau_ratio_median": [B["{}_tau_ratio".format(c)].median()
                                 for c in CLS],
            "tau_median": [B["{}_tau".format(c)].median() for c in CLS],
            "top2_frac_median": [B["{}_top_frac".format(c)].median()
                                 for c in CLS]})
        gt = gtex_breadth()
        summ["gtex_tissue_tau"] = summ["class"].map(gt)
        st.note("samples with >= {} cell types of >= {} cells: {}".format(
            MIN_TYPES, MIN_CELLS, len(B)))
        st.note("cell types kept: {}".format(
            int(S[S.n_cells >= MIN_CELLS].cell_type.nunique())))
        Dd = depth_diagnostic(S)
        Dd.to_csv(safe_write_path(os.path.join(
            OUT, "celltype_breadth_depth_diagnostic.tsv")), sep="\t",
            index=False)
        st.note("depth diagnostic (rho vs per-cell genes): "
                + "; ".join("{0}: det {1:+.2f} share {2:+.2f}".format(
                    r.klass, r.rho_det_vs_depth, r.rho_share_vs_depth)
                    for r in Dd.itertuples()))
        st.summary_md(os.path.join(OUT, "step_summaries", "56.md"),
                      "56 cell-type breadth",
                      df_to_md(summ, ".3f")
                      + "\n\nClass share of coding output and its depth "
                      "sensitivity\n\n" + df_to_md(Dd, ".3f")
                      + "\n\n" + df_to_md(T, ".4g", max_rows=40))


if __name__ == "__main__":
    main()
