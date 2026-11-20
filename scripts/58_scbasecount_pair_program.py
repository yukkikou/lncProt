#!/usr/bin/env python
"""58_scbasecount_pair_program.py — external replication of the pair programme.

Module 57 showed, in one glioma atlas, that a divergent lncRNA's cell-type
profile tracks its partner PCG more closely than three matched controls. That
result rests on a single cohort with a single harmonised annotation, so this
module repeats it on the 520-sample scBaseCount compendium, where the cell-type
labels come from many studies and only partially overlap.

The design differences forced by the data
  - labels are per sample and redundant, so a cell type enters only if it is
    annotated in at least MIN_SAMPLES samples
  - pooling cells across samples would confound a cell type with the studies
    that happen to carry it, so the profile is computed inside each sample and
    averaged over samples with equal weight
  - a percentile is taken inside each sample and cell type, so depth cancels
    within each sample as well as within each cell type

Outputs
  results/scbasecount/pair_program.tsv
  results/scbasecount/pair_program_tests.tsv
  results/scbasecount/pair_program_coverage.tsv
  results/scbasecount/step_summaries/58.md
"""
import argparse
import importlib.util
import os
import sys

import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

_spec = importlib.util.spec_from_file_location(
    "m57", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "57_sc_pair_program.py"))
m57 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m57)

OUT = os.path.join(WORKDIR, "results", "scbasecount")
ROOT = "/media/share/exchange/1610305236/scBaseCount"
MAN = os.path.join(ROOT, "manifest_design.tsv")
MIN_CELLS = 100     # cells of a cell type within a sample
MIN_SAMPLES = 20    # samples that must annotate a cell type
MIN_TYPES = 6       # cell types a sample must contribute for its profile


def gene_rates(row):
    """Phase 1. Gene-major CSC makes the per-gene count the indptr gap.

    No index array is read here, so the whole design costs one small read per
    file before the expensive pass.
    """
    with h5py.File(row["dest"], "r") as f:
        indptr = f["X"]["indptr"][:]
        nvar = int(f["X"].attrs["shape"][1])
    if len(indptr) != nvar + 1:
        raise SystemExit("unexpected layout in " + row["dest"])
    return np.diff(indptr).astype(np.int64)


def gene_by_type(row, genes):
    """Phase 2. Detection count per gene and cell type in one sample.

    Returns the percentile of each requested gene inside each annotated cell
    type, computed against all genes of that cell type, plus the cell count of
    every retained type.
    """
    with h5py.File(row["dest"], "r") as f:
        var = np.array([x.decode() if isinstance(x, bytes) else str(x)
                        for x in f["var"]["_index"][:]])
        nvar = len(var)
        g = f["obs"]["cell_type"]
        cats = [x.decode() if isinstance(x, bytes) else str(x)
                for x in g["categories"][:]]
        ctype = np.asarray(g["codes"][:], dtype=np.int64)
        nc = np.bincount(ctype, minlength=len(cats))
        keep = np.where(nc >= MIN_CELLS)[0]
        if len(keep) < 2:
            return None
        remap = -np.ones(len(cats), dtype=np.int64)
        remap[keep] = np.arange(len(keep))
        ct = remap[ctype]
        nk = len(keep)
        indptr = f["X"]["indptr"][:]
        idx_all = f["X"]["indices"]
        C = np.zeros(nvar * nk, dtype=np.int64)
        BLK = 2000
        for g0 in range(0, nvar, BLK):
            g1 = min(g0 + BLK, nvar)
            lo, hi = int(indptr[g0]), int(indptr[g1])
            if hi <= lo:
                continue
            cid = idx_all[lo:hi].astype(np.int64)
            sub = ct[cid]
            ok = sub >= 0
            if not ok.any():
                continue
            cid, sub = cid[ok], sub[ok]
            rep = np.repeat(np.arange(g0, g1, dtype=np.int64),
                            np.diff(indptr[g0:g1 + 1]))[ok]
            C += np.bincount(rep * nk + sub, minlength=nvar * nk)
        C = C.reshape(nvar, nk).astype(np.float64)
        C /= nc[keep][None, :]
        P = np.empty((len(genes), nk), dtype=np.float32)
        for j in range(nk):
            r = rankdata(C[:, j]) / nvar
            P[:, j] = r[genes]
    names = [cats[t] for t in keep]
    return P, names, nc[keep]


def worker(args):
    row, genes = args
    try:
        res = gene_by_type(row, genes)
    except Exception as exc:            # a damaged file must not kill the run
        return row["srx"], "error: " + str(exc), None
    # None without an error means the sample carries no usable annotation,
    # which is expected and distinct from a broken file.
    return row["srx"], None, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    with step("58_scbasecount_pair_program") as st:
        man = pd.read_csv(MAN, sep="\t")
        if args.limit:
            man = man.head(args.limit)
        rows = [dict(r) for _, r in man.iterrows()]
        st.note("samples={} cells={}".format(len(man), int(man.cells.sum())))

        # ---- phase 1: global detection rate per gene -----------------------
        from concurrent.futures import ProcessPoolExecutor
        counts = None
        with ProcessPoolExecutor(args.workers) as ex:
            for i, c in enumerate(ex.map(gene_rates, rows, chunksize=8)):
                counts = c if counts is None else counts + c
                if i % 100 == 0:
                    st.note("rates {}/{}".format(i + 1, len(rows)))
        nvar = len(counts)
        ncell = int(man.cells.sum())
        rate = pd.Series(counts / ncell, index=[
            "g%d" % i for i in range(nvar)])
        st.note("genes={} mean cells={}".format(nvar, ncell))

        # var order is shared by every file in the release, so the vector above
        # is aligned across samples. Pair construction needs the identifiers.
        with h5py.File(rows[0]["dest"], "r") as f:
            var = np.array([x.decode() if isinstance(x, bytes) else str(x)
                            for x in f["var"]["_index"][:]])
        rate.index = var
        pairs, var = m57.build_pairs(var, counts / ncell)
        st.note("testable pairs={}".format(len(pairs)))
        need = sorted({int(x) for t in pairs.itertuples()
                       for x in (t.a, t.b, t.c1, t.c3, t.c4) if x >= 0})
        st.note("genes needed={}".format(len(need)))
        gene_pos = {g: i for i, g in enumerate(need)}

        # ---- phase 2: per-sample, per-cell-type percentile -----------------
        acc = {}
        failures, unannotated = [], 0
        with ProcessPoolExecutor(args.workers) as ex:
            for i, (srx, err, res) in enumerate(
                    ex.map(worker, [(r, np.array(need)) for r in rows],
                           chunksize=4)):
                if res is None:
                    if err:
                        failures.append((srx, err))
                    else:
                        unannotated += 1
                    continue
                P, names, nc = res
                for j, nm in enumerate(names):
                    if nc[j] < MIN_CELLS:
                        continue
                    d = acc.setdefault(nm, [np.zeros(len(need)), 0])
                    d[0] += P[:, j]
                    d[1] += 1
                if i % 50 == 0:
                    st.note("profiles {}/{}".format(i + 1, len(rows)))
        st.note("samples without usable annotation={}".format(unannotated))
        st.note("files failed={}".format(len(failures)))
        for srx, err in failures[:5]:
            st.note("  {} {}".format(srx, err))

        cov = pd.DataFrame([dict(cell_type=k, n_samples=v[1])
                            for k, v in acc.items()]).sort_values(
            "n_samples", ascending=False)
        cov.to_csv(safe_write_path(os.path.join(
            OUT, "pair_program_coverage.tsv")), sep="\t", index=False)
        keep_types = [k for k, v in acc.items() if v[1] >= MIN_SAMPLES]
        st.note("cell types >= {} samples: {}".format(
            MIN_SAMPLES, len(keep_types)))
        if len(keep_types) < MIN_TYPES:
            st.note("too few cell types, stopping")
            return
        # Average the within-sample percentiles, one sample one vote.
        P = np.column_stack([acc[k][0] / max(acc[k][1], 1)
                             for k in keep_types])

        def prof(i):
            return P[gene_pos[i]]

        rows_out = []
        for t in pairs.itertuples():
            r = dict(pair="{}|{}".format(var[t.a], var[t.b]),
                     distance=t.distance,
                     rho_partner=spearmanr(prof(t.a), prof(t.b))[0],
                     rho_pcg_control=spearmanr(prof(t.a), prof(t.c3))[0],
                     rho_lnc_control=spearmanr(prof(t.c1), prof(t.b))[0],
                     rho_neighbour=(spearmanr(prof(t.a), prof(t.c4))[0]
                                    if t.c4 >= 0 else np.nan))
            rows_out.append(r)
        R = pd.DataFrame(rows_out)
        R.to_csv(safe_write_path(os.path.join(OUT, "pair_program.tsv")),
                 sep="\t", index=False)
        T = m57.run_tests(R, "scBaseCount {}".format(len(keep_types)))
        T.to_csv(safe_write_path(os.path.join(
            OUT, "pair_program_tests.tsv")), sep="\t", index=False)
        for r in T.itertuples():
            st.note("{} n={} d={:+.3f} [{:+.3f},{:+.3f}] frac+={:.2f} "
                    "p={:.3g}".format(r.comparison, r.n, r.median_diff,
                                      r.ci_low, r.ci_high,
                                      r.frac_left_higher, r.pval))
        st.summary_md(os.path.join(OUT, "step_summaries", "58.md"),
                      "58 scBaseCount pair programme",
                      df_to_md(T, ".4g") + "\n\nCell types used\n\n"
                      + df_to_md(cov.head(60), "4"))


if __name__ == "__main__":
    main()
