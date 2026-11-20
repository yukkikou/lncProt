#!/usr/bin/env python
"""57_sc_pair_program.py — do divergent pairs share a cell-type programme?

Every previous pair test asked whether two genes are detected in the same cell.
That quantity is dominated by per-cell depth, and the excess over a
distance-matched neighbour was null twice. This module changes the unit of
analysis from the cell to the cell type.

Each gene is compressed into a profile, its detection percentile in each
annotated cell type. Percentiles are taken within cell type, so a cell type's
sequencing depth cancels. The question is whether a divergent lncRNA's profile
is closer to its partner PCG than to the controls below.

Controls, all matched to the pair on a different axis
  c1  another lncRNA of matched detection rate, paired with the same partner
      (is the lncRNA itself specific, or would any lncRNA track the partner)
  c3  a PCG of matched detection rate at an unrelated locus, paired with the
      same lncRNA (is the signal just expression level)
  c4  a PCG at the same distance from the lncRNA TSS, on the same strand as the
      partner, within 100 kb (is the signal just the local neighbourhood)

Outputs
  results/sc_atlas/pair_celltype_program.tsv
  results/sc_atlas/pair_celltype_program_tests.tsv
  results/sc_atlas/step_summaries/57.md
"""
import os
import sys

import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
H5 = os.path.join(WORKDIR, "data", "scrna_atlas", "GBmap.h5ad")
OUT = os.path.join(RES, "sc_atlas")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
PAIR_MAX = 400
SEED = 42
MIN_CELLS = 100      # primary profile uses only well-sampled cell types
MIN_CELLS_LOOSE = 20  # sensitivity run keeps everything


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def read_cat(f, key):
    g = f["obs"][key]
    cats = [x.decode() if isinstance(x, bytes) else str(x)
            for x in g["categories"][:]]
    return np.asarray(g["codes"][:], dtype=np.int64), cats


def gene_by_celltype(f, nvar, ncat):
    """Count, per gene and cell type, how many cells detect the gene.

    One pass over the CSR matrix in cell blocks. The bincount key is
    gene * ncat * 2 + cell type * 2 + donor half, so each block touches
    nvar * ncat * 2 bins at most and the split-half replication below costs no
    extra read of the matrix.
    """
    ctype, cats = read_cat(f, "cell_type")
    donor, _ = read_cat(f, "donor_id")
    # Donors, not cells, define the halves, so a single large donor cannot sit
    # on both sides.
    half = (donor % 2).astype(np.int64)
    n = int(f["X"].attrs["shape"][0])
    indptr = f["X"]["indptr"][:]
    indices = f["X"]["indices"]
    G = np.zeros(nvar * ncat * 2, dtype=np.int64)
    CH = 20000
    for s in range(0, n, CH):
        e = min(s + CH, n)
        lo, hi = int(indptr[s]), int(indptr[e])
        idx = indices[lo:hi]
        rows = np.repeat(np.arange(e - s, dtype=np.int64),
                         np.diff(indptr[s:e + 1]))
        ct = ctype[s:e][rows] * 2 + half[s:e][rows]
        G += np.bincount(idx * (ncat * 2) + ct, minlength=nvar * ncat * 2)
    G = G.reshape(nvar, ncat, 2)
    nc = np.zeros((ncat, 2))
    for h in (0, 1):
        nc[:, h] = np.bincount(ctype[half == h], minlength=ncat)
    return G, nc, ctype, cats


def profiles(C, nc, keep):
    """Detection percentile of every gene inside every retained cell type.

    C is a genes x cell types count matrix and nc the matching cell counts.
    """
    nvar = C.shape[0]
    P = np.empty((nvar, len(keep)), dtype=np.float32)
    for j, t in enumerate(keep):
        P[:, j] = rankdata(C[:, t] / max(nc[t], 1)) / nvar
    return P


def build_pairs(var, D):
    meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                       sep="\t", index_col=0, low_memory=False)
    m2 = meta.assign(base=meta.index.str.split(".").str[0])
    m2 = m2.drop_duplicates("base").set_index("base").reindex(var)
    cls = m2["cls"]
    vidx = {g: i for i, g in enumerate(var)}
    dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
    rate = D

    cand = []
    for r in dp.itertuples():
        a = vidx.get(str(r.lnc_gid).split(".")[0])
        b = vidx.get(str(r.pcg_gid).split(".")[0])
        if a is not None and b is not None:
            cand.append((a, b, float(r.distance)))

    lnc_idx = np.where(cls.isin(CLASSES).values)[0]
    lnc_rate = rate[lnc_idx]
    order = np.argsort(lnc_rate)
    pcg_idx = np.where((m2["group"] == "pcg").values)[0]
    pcg_rate = rate[pcg_idx]
    porder = np.argsort(pcg_rate)

    def matched(target, exclude, pool, prate, pord):
        j = np.searchsorted(prate[pord], target)
        for k in range(1, 200):
            for jj in (j - k, j + k, j):
                if 0 <= jj < len(pord):
                    c = int(pool[pord[jj]])
                    if c not in exclude:
                        return c
        return int(pool[pord[0]])

    pcg_mask = (m2["group"] == "pcg").values
    chrom = m2["chrom"].astype(str).values
    tss = pd.to_numeric(m2["tss"], errors="coerce").values
    strand = m2["strand"].astype(str).values
    by_chr = {}
    for c in np.unique(chrom[pcg_mask]):
        ids = np.where(pcg_mask & (chrom == c))[0]
        by_chr[c] = ids[np.argsort(tss[ids])]

    def neighbour(a, b, d):
        ids = by_chr.get(chrom[a])
        if ids is None or not np.isfinite(tss[a]):
            return None, np.nan
        dist = np.abs(tss[ids] - tss[a])
        ok = (ids != b) & (strand[ids] == strand[b]) & (dist <= 1e5)
        cands, dd = ids[ok], dist[ok]
        if len(cands) == 0:
            return None, np.nan
        rr = np.abs(np.log2(rate[cands] / max(rate[b], 1e-6)))
        keep = rr < 1.0
        if keep.any():
            cands, dd = cands[keep], dd[keep]
        j = int(np.argmin(np.abs(dd - d)))
        return int(cands[j]), float(dd[j])

    out = []
    for a, b, d in cand[:PAIR_MAX]:
        if rate[a] < 0.005 or rate[b] < 0.005:
            continue
        c1 = matched(rate[a], {a, b}, lnc_idx, lnc_rate, order)
        c3 = matched(rate[b], {a, b, c1}, pcg_idx, pcg_rate, porder)
        c4, d4 = neighbour(a, b, d) if d > 0 else (None, np.nan)
        out.append(dict(a=a, b=b, c1=c1, c3=c3,
                        c4=-1 if c4 is None else c4, distance=d,
                        distance_control=d4))
    return pd.DataFrame(out), var


def run_tests(R, tag):
    rng = np.random.default_rng(SEED)
    rows = []

    def one(name, x, y, note):
        d = (x - y).dropna()
        if len(d) < 10:
            return
        w, p = wilcoxon(d, zero_method="wilcox")
        boots = [np.median(rng.choice(d.values, len(d))) for _ in range(1000)]
        rows.append(dict(run=tag, comparison=name, n=len(d),
                         median_rho_left=float(np.nanmedian(x)),
                         median_rho_right=float(np.nanmedian(y)),
                         median_diff=float(np.median(d)),
                         ci_low=float(np.percentile(boots, 2.5)),
                         ci_high=float(np.percentile(boots, 97.5)),
                         frac_left_higher=float((d > 0).mean()),
                         pval=float(p), note=note))

    one("partner vs distance-matched neighbour", R.rho_partner,
        R.rho_neighbour, "same locus, same strand, matched detection")
    one("partner vs rate-matched PCG", R.rho_partner, R.rho_pcg_control,
        "matched detection, unrelated locus")
    one("partner vs rate-matched lncRNA", R.rho_partner, R.rho_lnc_control,
        "same partner, another lncRNA of matched detection")
    R2 = R.assign(partner=R.pair.str.split("|").str[1])
    agg = R2.groupby("partner")[["rho_partner", "rho_neighbour",
                                 "rho_pcg_control", "rho_lnc_control"]].mean()
    one("partner-level mean, vs neighbour", agg.rho_partner,
        agg.rho_neighbour, "one row per partner PCG")
    T = pd.DataFrame(rows)
    T["fdr"] = np.nan
    for r, sub in T.groupby("run"):
        T.loc[sub.index, "fdr"] = bh_fdr(sub.pval.values)
    return T


def main():
    with step("57_sc_pair_program") as st:
        os.makedirs(OUT, exist_ok=True)
        f = h5py.File(H5, "r")
        var = np.array([x.decode() if isinstance(x, bytes) else str(x)
                        for x in f["var"]["_index"][:]])
        n = int(f["X"].attrs["shape"][0])
        nvar = int(f["X"].attrs["shape"][1])
        ncat = len(f["obs"]["cell_type"]["categories"])
        st.note("cells={} genes={} cell types={}".format(n, nvar, ncat))
        G3, nc2, ctype, cats = gene_by_celltype(f, nvar, ncat)
        f.close()
        G = G3.sum(axis=2)
        nc = nc2.sum(axis=1)
        D = G.sum(axis=1) / n          # overall detection rate per gene
        pairs, var = build_pairs(var, D)
        st.note("testable pairs={}".format(len(pairs)))
        st.note("donor halves: {} / {} cells".format(
            int(nc2.sum(axis=0)[0]), int(nc2.sum(axis=0)[1])))

        keep_hi = np.where(nc >= MIN_CELLS)[0]
        keep_all = np.where(nc >= MIN_CELLS_LOOSE)[0]
        st.note("cell types kept: {} (>= {} cells), {} (>= {} cells)".format(
            len(keep_hi), MIN_CELLS, len(keep_all), MIN_CELLS_LOOSE))

        tables = []
        for tag, keep in (("strict", keep_hi), ("loose", keep_all)):
            P = profiles(G, nc, keep)
            rows = []
            for r in pairs.itertuples():
                pa, pb = P[r.a], P[r.b]
                rows.append(dict(
                    pair="{}|{}".format(var[r.a], var[r.b]),
                    distance=r.distance,
                    distance_control=r.distance_control,
                    rho_partner=spearmanr(pa, pb)[0],
                    rho_pcg_control=spearmanr(pa, P[r.c3])[0],
                    rho_lnc_control=spearmanr(P[r.c1], pb)[0],
                    rho_neighbour=(spearmanr(pa, P[r.c4])[0]
                                   if r.c4 >= 0 else np.nan),
                    rate_lnc=float(D[r.a]), rate_partner=float(D[r.b]),
                    rate_pcg_control=float(D[r.c3]),
                    rate_neighbour=(float(D[r.c4]) if r.c4 >= 0 else np.nan),
                    # Spread of the partner and lncRNA profiles. A flat partner
                    # cannot be tracked by anything, so pairs are stratified on
                    # it before the correlation is read as shared regulation.
                    spread_partner=float(np.subtract(
                        *np.percentile(pb, [75, 25]))),
                    spread_lnc=float(np.subtract(
                        *np.percentile(pa, [75, 25])))))
            Rt = pd.DataFrame(rows)
            Rt.insert(0, "cell_types", len(keep))
            tables.append((tag, Rt))
        R = pd.concat([t for _, t in tables], ignore_index=True)
        R.to_csv(safe_write_path(os.path.join(
            OUT, "pair_celltype_program.tsv")), sep="\t", index=False)

        T = pd.concat([run_tests(t, tag) for tag, t in tables],
                      ignore_index=True)
        # Pairs whose partner changes between lineages are the only ones on
        # which a shared programme can be observed, so they are also tested on
        # their own. Pairs with a flat partner are retained above and dilute
        # every effect toward zero.
        strict = tables[0][1]
        cut = strict.spread_partner.median()
        top = strict[strict.spread_partner >= cut]
        T = pd.concat([T, run_tests(top, "strict, variable partner")],
                      ignore_index=True)
        T["partner_spread_cut"] = np.nan
        T.loc[T.run == "strict, variable partner", "partner_spread_cut"] = cut
        T.to_csv(safe_write_path(os.path.join(
            OUT, "pair_celltype_program_tests.tsv")), sep="\t", index=False)

        # Diagnostics. A difference that grows with a mismatch between the
        # controls and the true pair would mean the matching failed.
        diags = []
        for tag, t in tables:
            d = (t.rho_partner - t.rho_neighbour).dropna()
            sub = t.loc[d.index]
            diags.append(dict(
                run=tag, n_pairs=len(t),
                rho_partner_median=float(t.rho_partner.median()),
                rho_neighbour_median=float(t.rho_neighbour.median()),
                rho_pcg_median=float(t.rho_pcg_control.median()),
                rho_lnc_median=float(t.rho_lnc_control.median()),
                rho_dist_equal=float(spearmanr(
                    d, sub.distance_control - sub.distance)[0]),
                rho_rate_mismatch=float(spearmanr(
                    d, np.abs(np.log2(sub.rate_neighbour /
                                      sub.rate_partner)))[0]),
                rho_lnc_rate=float(spearmanr(d, sub.rate_lnc)[0])))
        Dg = pd.DataFrame(diags)
        Dg.to_csv(safe_write_path(os.path.join(
            OUT, "pair_celltype_program_diagnostics.tsv")), sep="\t",
            index=False)
        st_ = strict.copy()
        st_["delta"] = st_.rho_partner - st_.rho_neighbour
        st_["q"] = pd.qcut(st_.spread_partner, 4, labels=False)
        strat = st_.groupby("q").agg(
            n=("delta", "count"), spread=("spread_partner", "median"),
            delta=("delta", "median"),
            partner_rho=("rho_partner", "median"),
            neighbour_rho=("rho_neighbour", "median")).reset_index()
        strat.to_csv(safe_write_path(os.path.join(
            OUT, "pair_celltype_program_by_partner_spread.tsv")), sep="\t",
            index=False)

        # Split-half replication over donors. The same pair should show the
        # same profile agreement when the cell types are measured in a
        # disjoint set of donors, otherwise the per-pair correlation is noise.
        both = np.where((nc2[:, 0] >= 50) & (nc2[:, 1] >= 50))[0]
        P0 = profiles(G3[:, :, 0], nc2[:, 0], both)
        P1 = profiles(G3[:, :, 1], nc2[:, 1], both)
        sh = []
        for r in pairs.itertuples():
            row = dict(pair="{}|{}".format(var[r.a], var[r.b]),
                       half0=spearmanr(P0[r.a], P0[r.b])[0],
                       half1=spearmanr(P1[r.a], P1[r.b])[0])
            if r.c4 >= 0:
                row["delta0"] = row["half0"] - spearmanr(P0[r.a], P0[r.c4])[0]
                row["delta1"] = row["half1"] - spearmanr(P1[r.a], P1[r.c4])[0]
            sh.append(row)
        SH = pd.DataFrame(sh)
        SH.to_csv(safe_write_path(os.path.join(
            OUT, "pair_celltype_program_split_half.tsv")), sep="\t",
            index=False)
        rep = spearmanr(SH.half0, SH.half1)
        dh = SH[["delta0", "delta1"]].dropna()
        w0, p0 = wilcoxon(dh.delta0)
        w1, p1 = wilcoxon(dh.delta1)
        st.note("split-half: {} cell types in both halves, {} pairs; "
                "per-pair rho agreement rho={:.2f} p={:.2g}".format(
                    len(both), len(SH), rep[0], rep[1]))
        st.note("delta vs neighbour, half0 median {:.3f} (P={:.3g}), "
                "half1 median {:.3f} (P={:.3g}), n={}".format(
                    dh.delta0.median(), p0, dh.delta1.median(), p1, len(dh)))

        # Leave one cell type out. The atlas is dominated by malignant cells,
        # so the effect must not depend on that one category.
        loo = []
        for t in keep_hi:
            kk = keep_hi[keep_hi != t]
            Pk = profiles(G, nc, kk)
            d = [spearmanr(Pk[r.a], Pk[r.b])[0] -
                 spearmanr(Pk[r.a], Pk[r.c4])[0]
                 for r in pairs.itertuples() if r.c4 >= 0]
            loo.append(dict(held_out=cats[t], n_types=len(kk), n=len(d),
                            median_delta=float(np.median(d))))
        LOO = pd.DataFrame(loo).sort_values("median_delta")
        LOO.to_csv(safe_write_path(os.path.join(
            OUT, "pair_celltype_program_loo.tsv")), sep="\t", index=False)
        st.note("leave-one-out delta range {:.3f} ({}) to {:.3f} ({})".format(
            LOO.median_delta.iloc[0], LOO.held_out.iloc[0],
            LOO.median_delta.iloc[-1], LOO.held_out.iloc[-1]))

        for r in T.itertuples():
            st.note("[{}] {} n={} d={:+.3f} [{:+.3f},{:+.3f}] frac+={:.2f} "
                    "p={:.3g} fdr={:.3g}".format(
                        r.run, r.comparison, r.n, r.median_diff, r.ci_low,
                        r.ci_high, r.frac_left_higher, r.pval, r.fdr))
        st.summary_md(os.path.join(OUT, "step_summaries", "57.md"),
                      "57 divergent-pair cell-type programme",
                      df_to_md(T, ".4g")
                      + "\n\nDiagnostics (a large absolute value would mean "
                      "the matching failed)\n\n" + df_to_md(Dg, ".3f")
                      + "\n\nPair difference by quartile of partner profile "
                      "spread (strict run)\n\n" + df_to_md(strat, ".3f"))


if __name__ == "__main__":
    main()
