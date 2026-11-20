#!/usr/bin/env python
"""52_scbasecount_classes.py — birth-class detection in the scBaseCount subset.

Scans the 520 human matrices downloaded by module 51 and expresses, for every
sample, how often each birth class is detected per cell. Detection is made
comparable across cells of different depth by subtracting the fraction of all
genes detected in the same cell (the "adjusted" measure used for the glioma
atlas), and samples are the unit of analysis.

Outputs
  results/scbasecount/sample_class_detection.tsv   per sample means
  results/scbasecount/group_tests.tsv              tumour vs normal per group
  results/scbasecount/step_summaries/52.md
"""
import argparse
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

OUT = os.path.join(WORKDIR, "results", "scbasecount")
ROOT = "/media/share/exchange/1610305236/scBaseCount"
LNC5 = ["divergent", "antisense", "intronic", "sense_overlapping",
        "intergenic"]
# 0 = untracked, then PCG, the five classes, pseudogene and LTR-derived lncRNA
CODES = {"PCG": 1, "divergent": 2, "antisense": 3, "intronic": 4,
         "sense_overlapping": 5, "intergenic": 6, "pseudogene": 7,
         "LTR_lnc": 8}
NC = 9


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def class_codes():
    meta = pd.read_csv(os.path.join(WORKDIR, "results",
                                    "gene_master_table.tsv"),
                       sep="\t", index_col=0, low_memory=False)
    base = meta.index.str.split(".").str[0]
    meta = meta.assign(base=base.values).drop_duplicates("base")
    te = pd.read_csv(os.path.join(WORKDIR, "results",
                                  "te_promoter_stats.tsv"),
                     sep="\t", index_col=0, low_memory=False)[["tss_in_ltr"]]
    te.index = te.index.str.split(".").str[0]
    te = te[~te.index.duplicated()]
    m = meta.set_index("base").join(te, on="base")
    cls = m["cls"].astype(str)
    code = pd.Series(0, index=m.index, dtype=np.int8)
    code[cls.isin(LNC5)] = cls[cls.isin(LNC5)].map(CODES).astype(np.int8)
    code[cls == "pseudogene"] = CODES["pseudogene"]
    code[m["group"] == "pcg"] = CODES["PCG"]
    ltr = cls.isin(LNC5) & (m["tss_in_ltr"] == 1)
    code[ltr] = CODES["LTR_lnc"]
    return code


def one_sample(row, code_map):
    """Class detection for one sample.

    The release stores X as a CSC matrix whose shape is (cells, genes), so the
    compressed axis is genes. Non-zeros of one gene are its detected cells, and
    the per-cell counts are accumulated by grouping on the index array. A CSR
    layout is still handled in case a sample is ever stored that way.
    """
    dest, cells = row["dest"], int(row["cells"])
    t0 = time.time()
    with h5py.File(dest, "r") as f:
        var = np.array([x.decode() if isinstance(x, bytes) else str(x)
                        for x in f["var"]["_index"][:]])
        gcode = code_map.reindex(var).fillna(0).astype(np.int8).values
        nvar = len(var)
        sizes = np.bincount(gcode, minlength=NC).astype(float)
        indptr = f["X"]["indptr"][:]
        idx_all = f["X"]["indices"]
        shape = [int(x) for x in f["X"].attrs["shape"]]
        gene_major = (len(indptr) == nvar + 1) and shape[1] == nvar
        n = shape[0] if gene_major else shape[1]
        acc = np.zeros(NC)
        per_cell = np.zeros(n, dtype=np.int64)
        if gene_major:
            # walk gene blocks so the repeated index array stays small
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
                                       minlength=n * NC
                                       ).reshape(n, NC).sum(axis=0)
        else:
            trk = gcode > 0
            CH = 50000
            for s in range(0, n, CH):
                e = min(s + CH, n)
                lo, hi = int(indptr[s]), int(indptr[e])
                idx = idx_all[lo:hi]
                rows = np.repeat(np.arange(e - s, dtype=np.int64),
                                 np.diff(indptr[s:e + 1]))
                per_cell[s:e] = np.bincount(rows, minlength=e - s)
                keep = trk[idx]
                bins = rows[keep] * NC + gcode[idx[keep]]
                acc += np.bincount(bins, minlength=(e - s) * NC).reshape(
                    e - s, NC).sum(axis=0)
        nnz = per_cell.astype(float)
    out = {"srx": row["srx"], "group": row["group"], "arm": row["arm"],
           "n_cells_file": n, "n_cells_manifest": cells,
           "n_genes": nvar, "mean_genes_per_cell": float(nnz.mean()),
           "seconds": round(time.time() - t0, 1)}
    overall = float((nnz / nvar).mean())
    out["overall_det"] = overall
    for name, i in CODES.items():
        det = float((acc[i] / sizes[i] / n).mean()) if sizes[i] else np.nan
        out["{}_det".format(name)] = det
        out["{}_adj".format(name)] = det - overall
    return out


def tests(samples):
    rows = []
    metrics = ["overall_det"] + ["{}_adj".format(k) for k in CODES]
    for grp, sub in samples.groupby("group"):
        a = sub[sub.arm == "tumour"]
        b = sub[sub.arm == "normal"]
        if len(a) < 5 or len(b) < 5:
            continue
        ps = []
        for m in metrics:
            x, y = a[m].dropna(), b[m].dropna()
            if len(x) < 5 or len(y) < 5:
                continue
            u, p = mannwhitneyu(x, y, alternative="two-sided")
            d = 2 * u / (len(x) * len(y)) - 1
            rows.append(dict(group=grp, metric=m, n_tumour=len(x),
                             n_normal=len(y), median_tumour=float(x.median()),
                             median_normal=float(y.median()),
                             cliffs_delta=float(d), pval=float(p)))
    T = pd.DataFrame(rows)
    if len(T):
        T["fdr"] = np.nan
        for grp, sub in T.groupby("group"):
            T.loc[sub.index, "fdr"] = bh_fdr(sub.pval.values)
    return T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()

    with step("52_scbasecount_classes") as st:
        os.makedirs(OUT, exist_ok=True)
        man = pd.read_csv(os.path.join(ROOT, "manifest_design.tsv"), sep="\t")
        if args.limit:
            man = man.head(args.limit)
        st.note("samples={} cells={}".format(len(man), int(man.cells.sum())))
        code_map = class_codes()
        rows = [dict(r) for _, r in man.iterrows()]
        if args.workers > 1:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(args.workers) as ex:
                res = list(ex.map(one_sample, rows,
                                  [code_map] * len(rows), chunksize=4))
        else:
            res = []
            for i, r in enumerate(rows):
                res.append(one_sample(r, code_map))
                if i % 20 == 0:
                    st.note("done {}/{} ({:.0f}s)".format(
                        i + 1, len(rows), res[-1]["seconds"]))
        S = pd.DataFrame(res)
        S.to_csv(safe_write_path(os.path.join(
            OUT, "sample_class_detection.tsv")), sep="\t", index=False)
        T = tests(S)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "group_tests.tsv")), sep="\t", index=False)
        st.note("groups tested={}".format(T.group.nunique() if len(T) else 0))
        if len(T):
            sig = T[(T.fdr < 0.05) & T.metric.str.endswith("_adj")]
            st.note("significant class shifts: {}".format(len(sig)))
        st.summary_md(os.path.join(OUT, "step_summaries", "52.md"),
                      "52 scBaseCount class detection",
                      "Per-sample class detection, then tumour versus normal "
                      "within each cancer group.\n\n"
                      + df_to_md(T, ".4f", max_rows=60) if len(T) else "")


if __name__ == "__main__":
    main()
