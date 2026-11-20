#!/usr/bin/env python
"""61_cross_species_hierarchy.py — class hierarchy at scale, depth matched.

Module 55 compared 120 mouse samples with 520 human samples and reported that
the coarse order held in 90% versus 95%. Those two numbers were not comparable,
because mouse and human samples differ in per-cell depth and every class
detection rate moves with depth. This module expands the mouse set to over 600
samples across 121 tissues and puts the comparison on a matched footing.

It does three things
  1. per-sample class detection for the expanded mouse design
  2. concordance of the coarse order as a function of per-cell depth, in mouse
     and in human, so the two species are compared inside the same depth band
  3. a tissue-level summary over the 121 mouse tissues

Outputs
  results/scbasecount/mouse_sample_detection_expand.tsv
  results/scbasecount/hierarchy_depth_bins.tsv
  results/scbasecount/hierarchy_tissue_summary.tsv
  results/scbasecount/step_summaries/61.md
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

_spec = importlib.util.spec_from_file_location(
    "m55", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "55_scbasecount_mouse.py"))
m55 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m55)
m52 = m55.m52

OUT = os.path.join(WORKDIR, "results", "scbasecount")
ROOT = "/media/share/exchange/1610305236/scBaseCount"
MAN = os.path.join(ROOT, "manifest_mouse_expand.tsv")
DET = ["PCG_det", "sense_overlapping_det", "divergent_det", "antisense_det",
       "intergenic_det", "intronic_det"]
BINS = [0, 1000, 2000, 3000, 4000, 6000, 10000, 1e9]


def order_metrics(S):
    """Coarse-order indicators used throughout the manuscript."""
    R = S[DET].rank(axis=1, ascending=False)
    S = S.copy()
    S["top3_ok"] = ((R.sense_overlapping_det < R.divergent_det)
                    & (R.divergent_det < R.antisense_det))
    S["last2_ok"] = ((R[["intergenic_det", "intronic_det"]].max(axis=1) == 6)
                     & (R[["intergenic_det", "intronic_det"]].min(axis=1) == 5))
    S["intergenic_first"] = R.intergenic_det < R.intronic_det
    return S


def depth_bins(S, label):
    S = S.dropna(subset=["mean_genes_per_cell"])
    b = pd.cut(S.mean_genes_per_cell, BINS, right=False)
    rows = []
    for k, g in S.groupby(b, observed=True):
        if len(g) < 20:
            continue
        rows.append(dict(species=label, depth_bin=str(k), n=len(g),
                         median_depth=float(g.mean_genes_per_cell.median()),
                         top3_ok=float(g.top3_ok.mean()),
                         last2_ok=float(g.last2_ok.mean()),
                         intergenic_first=float(g.intergenic_first.mean())))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    with step("61_cross_species_hierarchy") as st:
        man = pd.read_csv(MAN, sep="\t")
        if args.limit:
            man = man.head(args.limit)
        # tolerate a partially finished download
        man = man[[os.path.exists(p) for p in man.dest]]
        st.note("samples={} cells={}".format(len(man), int(man.cells.sum())))
        code_map = m55.mouse_codes()
        rows = [dict(r) for _, r in man.iterrows()]
        if args.workers > 1:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(args.workers) as ex:
                res = list(ex.map(m52.one_sample, rows,
                                  [code_map] * len(rows), chunksize=4))
        else:
            res = [m52.one_sample(r, code_map) for r in rows]
        S = pd.DataFrame(res)
        S["tissue"] = S.group.str.replace("mouse_", "", regex=False)
        S.to_csv(safe_write_path(os.path.join(
            OUT, "mouse_sample_detection_expand.tsv")), sep="\t", index=False)
        st.note("scanned {} samples, {} cells".format(
            len(S), int(S.n_cells_file.sum())))

        # the 120-sample pilot is a subset of these tissues but a disjoint set
        # of accessions, so both are used
        old = pd.read_csv(os.path.join(OUT, "mouse_sample_detection.tsv"),
                          sep="\t")
        old["tissue"] = old.group.str.replace("mouse_", "", regex=False)
        A = pd.concat([S, old], ignore_index=True).drop_duplicates("srx")
        A = order_metrics(A)
        st.note("combined mouse samples={} cells={}".format(
            len(A), int(A.n_cells_file.sum())))
        st.note("top3 order holds in {:.1%}, last two classes last in "
                "{:.1%}".format(A.top3_ok.mean(), A.last2_ok.mean()))

        hu = pd.read_csv(os.path.join(OUT, "sample_class_detection.tsv"),
                         sep="\t")
        hu = order_metrics(hu)
        st.note("human n={} top3 {:.1%} last2 {:.1%}".format(
            len(hu), hu.top3_ok.mean(), hu.last2_ok.mean()))

        B = pd.concat([depth_bins(A, "mouse"), depth_bins(hu, "human")],
                      ignore_index=True)
        B.to_csv(safe_write_path(os.path.join(
            OUT, "hierarchy_depth_bins.tsv")), sep="\t", index=False)

        # matched-depth statement: restrict both species to the overlapping
        # interquartile range of per-cell depth
        lo = max(A.mean_genes_per_cell.quantile(0.25),
                 hu.mean_genes_per_cell.quantile(0.25))
        hi = min(A.mean_genes_per_cell.quantile(0.75),
                 hu.mean_genes_per_cell.quantile(0.75))
        am = A[(A.mean_genes_per_cell >= lo) & (A.mean_genes_per_cell <= hi)]
        hm = hu[(hu.mean_genes_per_cell >= lo) & (hu.mean_genes_per_cell <= hi)]
        st.note("matched depth {:.0f}-{:.0f} genes/cell: mouse n={} top3 "
                "{:.1%} last2 {:.1%}; human n={} top3 {:.1%} last2 {:.1%}"
                .format(lo, hi, len(am), am.top3_ok.mean(), am.last2_ok.mean(),
                        len(hm), hm.top3_ok.mean(), hm.last2_ok.mean()))

        T = A.groupby("tissue").agg(
            n=("srx", "size"),
            top3_ok=("top3_ok", "mean"),
            last2_ok=("last2_ok", "mean"),
            median_depth=("mean_genes_per_cell", "median")).reset_index()
        T.to_csv(safe_write_path(os.path.join(
            OUT, "hierarchy_tissue_summary.tsv")), sep="\t", index=False)
        st.note("tissues: {} with >=4 samples: {}".format(
            len(T), int((T.n >= 4).sum())))

        # Do tissues really differ, or is the spread just binomial sampling at
        # n=4? The observed between-tissue variance is compared with a null in
        # which every sample keeps its own count but shares the pooled rate.
        Tt = T[T.n >= 4]
        if len(Tt) >= 10:
            obs = float(np.var(Tt.top3_ok, ddof=1))
            rng = np.random.default_rng(42)
            p0 = float(A.top3_ok.mean())
            null = []
            for _ in range(2000):
                sim = [rng.binomial(int(n), p0) / n for n in Tt.n]
                null.append(np.var(sim, ddof=1))
            pv = float(np.mean(np.array(null) >= obs))
            st.note("between-tissue variance in top3 order {:.4f} versus "
                    "binomial null {:.4f} (P={:.3f}, {} tissues)".format(
                        obs, float(np.mean(null)), pv, len(Tt)))

        st.summary_md(os.path.join(OUT, "step_summaries", "61.md"),
                      "61 cross-species class hierarchy at scale",
                      "Concordance by per-cell depth\n\n" + df_to_md(B, ".3f")
                      + "\n\nTissues with the least consistent order\n\n"
                      + df_to_md(T.sort_values("top3_ok").head(12), ".3f"))


if __name__ == "__main__":
    main()
