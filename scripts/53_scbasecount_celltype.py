#!/usr/bin/env python
"""53_scbasecount_celltype.py — cell-type-matched tumour versus normal contrast.

The sample-level comparison of module 52 is confounded, because tumour and
normal samples differ in per-cell depth and in tissue composition. This module
repeats the scan but keeps cells separated by annotated cell type, so the same
lineage can be compared across arms. A comparison is only reported when the
cell type occurs in at least five samples of each arm within a cancer group,
and every model adjusts for the mean number of genes detected per cell.

Outputs
  results/scbasecount/sample_celltype_detection.tsv
  results/scbasecount/celltype_matched_tests.tsv
  results/scbasecount/step_summaries/53.md
"""
import argparse
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "m52", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "52_scbasecount_classes.py"))
m52 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m52)
CODES, NC = m52.CODES, m52.NC

OUT = os.path.join(WORKDIR, "results", "scbasecount")
ROOT = "/media/share/exchange/1610305236/scBaseCount"


def one_sample(row, code_map):
    dest = row["dest"]
    t0 = time.time()
    with h5py.File(dest, "r") as f:
        var = np.array([x.decode() if isinstance(x, bytes) else str(x)
                        for x in f["var"]["_index"][:]])
        gcode = code_map.reindex(var).fillna(0).astype(np.int8).values
        nvar = len(var)
        sizes = np.bincount(gcode, minlength=NC).astype(float)
        g = f["obs"]["cell_type"]
        cats = [x.decode() if isinstance(x, bytes) else str(x)
                for x in g["categories"][:]]
        ctype = np.asarray(g["codes"][:])
        indptr = f["X"]["indptr"][:]
        idx_all = f["X"]["indices"]
        shape = [int(x) for x in f["X"].attrs["shape"]]
        gene_major = (len(indptr) == nvar + 1) and shape[1] == nvar
        n = shape[0]
        acc = np.zeros(n * NC)
        per_cell = np.zeros(n, dtype=np.int64)
        if gene_major:
            BLK = 2000
            for g0 in range(0, nvar, BLK):
                g1 = min(g0 + BLK, nvar)
                lo, hi = int(indptr[g0]), int(indptr[g1])
                if hi <= lo:
                    continue
                cid = idx_all[lo:hi].astype(np.int64)
                rep = np.repeat(gcode[g0:g1], np.diff(indptr[g0:g1 + 1]))
                per_cell += np.bincount(cid, minlength=n)
                keep = rep > 0
                if keep.any():
                    acc += np.bincount(cid[keep] * NC + rep[keep],
                                       minlength=n * NC)
        else:
            raise SystemExit("unexpected CSR layout in " + dest)
    acc = acc.reshape(n, NC)
    ngene = per_cell.astype(float)
    rows = []
    for t, name in enumerate(cats):
        sel = ctype == t
        k = int(sel.sum())
        if k < 50:
            continue
        rows.append(dict(srx=row["srx"], group=row["group"], arm=row["arm"],
                         cell_type=name, n_cells=k,
                         mean_genes_per_cell=float(ngene[sel].mean()),
                         seconds=round(time.time() - t0, 1)))
        r = rows[-1]
        overall = float((ngene[sel] / nvar).mean())
        r["overall_det"] = overall
        for nm, i in CODES.items():
            det = float((acc[sel, i] / sizes[i]).mean()) if sizes[i] else np.nan
            r["{}_det".format(nm)] = det
            r["{}_adj".format(nm)] = det - overall
    return rows


def matched_tests(S):
    rows = []
    S = S[~S.cell_type.isin(["", "unknown", "nan"])]
    for (grp, ct), sub in S.groupby(["group", "cell_type"]):
        a = sub[sub.arm == "tumour"]
        b = sub[sub.arm == "normal"]
        if len(a) < 5 or len(b) < 5:
            continue
        for m in ["PCG_adj", "divergent_adj", "antisense_adj", "intronic_adj",
                  "sense_overlapping_adj", "intergenic_adj"]:
            d = sub.dropna(subset=[m]).copy()
            if len(d[d.arm == "tumour"]) < 5 or len(d[d.arm == "normal"]) < 5:
                continue
            d["tum"] = (d.arm == "tumour").astype(int)
            mod = smf.ols("{} ~ tum + mean_genes_per_cell".format(m),
                          data=d).fit()
            rows.append(dict(group=grp, cell_type=ct, metric=m,
                             n_tumour=int((d.arm == "tumour").sum()),
                             n_normal=int((d.arm == "normal").sum()),
                             beta=float(mod.params["tum"]),
                             pval=float(mod.pvalues["tum"])))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    with step("53_scbasecount_celltype") as st:
        os.makedirs(OUT, exist_ok=True)
        man = pd.read_csv(os.path.join(ROOT, "manifest_design.tsv"), sep="\t")
        if args.limit:
            man = man.head(args.limit)
        code_map = m52.class_codes()
        rows = [dict(r) for _, r in man.iterrows()]
        if args.workers > 1:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(args.workers) as ex:
                res = list(ex.map(one_sample, rows,
                                  [code_map] * len(rows), chunksize=4))
        else:
            res = [one_sample(r, code_map) for r in rows]
        S = pd.DataFrame([x for sub in res for x in sub])
        S.to_csv(safe_write_path(os.path.join(
            OUT, "sample_celltype_detection.tsv")), sep="\t", index=False)
        T = matched_tests(S)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "celltype_matched_tests.tsv")), sep="\t", index=False)
        st.note("rows={} comparisons={}".format(len(S), len(T)))
        st.summary_md(os.path.join(OUT, "step_summaries", "53.md"),
                      "53 cell-type-matched contrast",
                      df_to_md(T, ".4f", max_rows=80) if len(T) else "")


if __name__ == "__main__":
    main()
