#!/usr/bin/env python
"""72_sc_tissue_disease_specificity.py — is single-cell specificity the same
axis as GTEx tissue specificity and TCGA disease specificity?

The manuscript now uses three specificity readouts that were collected in
different systems. Single cells give cell-type specificity, GTEx gives tissue
specificity and TCGA gives disease (up-regulation) breadth. Nobody has checked
whether, gene by gene, they describe the same property, and if they do, whether
birth class or transcript age adds anything beyond detection rate.

This module computes per-gene cell-type specificity from the public glioma
atlas, joins it to the GTEx tau, the pan-cancer up count and the evolutionary
annotations, and asks four questions.

  1. Do genes that are cell-type specific in single cells also have high GTEx
     tau? Spearman per gene, then partialled on detection rate.
  2. Do the birth classes differ in cell-type specificity, and does that
     difference survive matching each lncRNA to a PCG of the same detection
     rate? tau rises as detection falls, so the raw class order cannot be read
     as biology.
  3. Within lncRNAs, does cell-type specificity track transcript age or
     promoter constraint once detection is held constant?
  4. Does cell-type specificity add anything to GTEx expression breadth in
     predicting pan-cancer up-regulation?

Outputs
  results/sc_atlas/gene_celltype_specificity.tsv
  results/sc_atlas/cross_system_specificity_tests.tsv
  results/step_summaries/72_sc_tissue_disease_specificity.md
"""
import os
import sys

import h5py
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
H5 = os.path.join(WORKDIR, "data", "scrna_atlas", "GBmap.h5ad")
OUT = os.path.join(RES, "sc_atlas")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
MIN_CELLS = 100
CALIPER = 0.1
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


def read_cat(f, key):
    g = f["obs"][key]
    cats = [x.decode() if isinstance(x, bytes) else str(x)
            for x in g["categories"][:]]
    return np.asarray(g["codes"][:], dtype=np.int64), cats


def gene_by_celltype(f, nvar, ncat):
    """Count detected cells per gene and cell type in one pass.

    Same construction as 57_sc_pair_program.py, without the donor halves.
    """
    ctype, cats = read_cat(f, "cell_type")
    n = int(f["X"].attrs["shape"][0])
    indptr = f["X"]["indptr"][:]
    indices = f["X"]["indices"]
    G = np.zeros(nvar * ncat, dtype=np.int64)
    CH = 20000
    for s in range(0, n, CH):
        e = min(s + CH, n)
        lo, hi = int(indptr[s]), int(indptr[e])
        idx = indices[lo:hi]
        rows = np.repeat(np.arange(e - s, dtype=np.int64),
                         np.diff(indptr[s:e + 1]))
        ct = ctype[s:e][rows]
        G += np.bincount(idx * ncat + ct, minlength=nvar * ncat)
    G = G.reshape(nvar, ncat)
    nc = np.bincount(ctype, minlength=ncat)
    return G, nc, cats


def tau_specificity(F):
    """Gini-style specificity over cell types, 0 uniform and 1 one cell type."""
    mx = F.max(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (1.0 - F / mx[:, None]).sum(axis=1) / (F.shape[1] - 1)
    t[~np.isfinite(mx) | (mx <= 0)] = np.nan
    return np.clip(t, 0, 1)


def partial_spearman(x, y, z):
    """Spearman on rank residuals so detection does not drive the association."""
    d = pd.DataFrame(dict(x=x, y=y, z=z)).dropna()
    if len(d) < 20:
        return np.nan, np.nan, len(d)
    rx = d.x.rank().values
    ry = d.y.rank().values
    rz = d.z.rank().values
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    ex = rx - np.polyval(bx, rz)
    ey = ry - np.polyval(by, rz)
    r, p = spearmanr(ex, ey)
    return float(r), float(p), len(d)


def match_one_to_one(target, pool, caliper):
    used = np.zeros(len(pool), dtype=bool)
    order = np.argsort(pool)
    ps = pool[order]
    out = np.full(len(target), -1, dtype=int)
    for i in np.argsort(target):
        j = np.searchsorted(ps, target[i])
        best, bd = -1, np.inf
        for k in (j - 1, j, j + 1, j - 2, j + 2):
            if 0 <= k < len(ps):
                pi = order[k]
                if used[pi]:
                    continue
                dd = abs(ps[k] - target[i])
                if dd < bd and dd <= caliper:
                    best, bd = pi, dd
        if best >= 0:
            used[best] = True
            out[i] = best
    return out


def main():
    with step("72_sc_tissue_disease_specificity") as st:
        os.makedirs(OUT, exist_ok=True)
        f = h5py.File(H5, "r")
        var = np.array([x.decode() if isinstance(x, bytes) else str(x)
                        for x in f["var"]["_index"][:]])
        n = int(f["X"].attrs["shape"][0])
        nvar = int(f["X"].attrs["shape"][1])
        ncat = len(f["obs"]["cell_type"]["categories"])
        G, nc, cats = gene_by_celltype(f, nvar, ncat)
        f.close()
        keep = np.where(nc >= MIN_CELLS)[0]
        st.note("cells={} genes={} cell types={} kept={} (>= {} cells)".format(
            n, nvar, ncat, len(keep), MIN_CELLS))

        F = G[:, keep] / np.maximum(nc[keep], 1)[None, :]
        rate = G.sum(axis=1) / n
        tau_c = tau_specificity(F)
        breadth_c = (F > 0).sum(axis=1)
        top = np.array(cats)[keep][np.argmax(F, axis=1)] if len(keep) else \
            np.array([""] * nvar)

        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        m2 = meta.assign(base=meta.index.str.split(".").str[0])
        m2 = m2.drop_duplicates("base").set_index("base")
        gtex = pd.read_csv(os.path.join(
            RES, "tissue_disease", "tables", "gtex_gene_tau.tsv"), sep="\t")
        gtex = gtex.rename(columns={"Unnamed: 0": "gene_id"})
        cb = pd.read_csv(os.path.join(
            RES, "cancer_bridge", "C_bridge_pergene.tsv"), sep="\t",
            index_col=0)
        cb["base"] = cb.base_id.astype(str)
        cb = cb.drop_duplicates("base").set_index("base")
        ev = pd.read_csv(os.path.join(
            RES, "comparative", "E_transcript_age.tsv"), sep="\t",
            index_col=0)
        ev["base"] = ev.index.astype(str).str.split(".").str[0]
        ev = ev.drop_duplicates("base").set_index("base")
        ph = pd.read_csv(os.path.join(
            RES, "phylop", "phylop_per_gene.tsv"), sep="\t", index_col=0)
        ph["base"] = ph.index.astype(str).str.split(".").str[0]
        ph = ph.drop_duplicates("base").set_index("base")

        V = pd.DataFrame(dict(
            gene_id=var,
            base=[str(g).split(".")[0] for g in var],
            det_rate=rate, n_celltypes=len(keep),
            tau_celltype=tau_c, breadth_celltype=breadth_c,
            top_celltype=top))
        V = V.drop_duplicates("base").set_index("base")
        V = V.join(m2[["cls", "group"]])
        g2 = gtex.assign(base=gtex.gene_id.astype(str).str.split(".").str[0])
        g2 = g2.drop_duplicates("base").set_index("base")
        V = V.join(g2[["tau_gtex", "max_tissue_mean"]])
        V = V.join(cb[["n_up", "n_de", "breadth"]].rename(
            columns={"breadth": "gtex_breadth_count"}))
        V = V.join(ev[["age_level", "age_class"]])
        V = V.join(ph[["phylop_mean", "phast_mean", "tss_in_ltr"]])
        V.to_csv(safe_write_path(os.path.join(
            OUT, "gene_celltype_specificity.tsv")), sep="\t")

        rows = []

        def add(family, test, subset, r, p, n_, note):
            rows.append(dict(family=family, test=test, subset=subset,
                             statistic=r, pval=p, n=n_, note=note))

        # ---- 1. does cell-type specificity agree with GTEx tissue tau? -----
        for grp in ("all", "lnc", "pcg"):
            d = V if grp == "all" else V[V.group == grp]
            r, p = spearmanr(d.tau_celltype, d.tau_gtex, nan_policy="omit")
            add("concordance", "Spearman tau_celltype vs tau_gtex", grp,
                r, p, int(d[["tau_celltype", "tau_gtex"]].dropna().shape[0]),
                "raw")
            r2, p2, n2 = partial_spearman(d.tau_celltype, d.tau_gtex,
                                          d.det_rate)
            add("concordance", "partial Spearman controlling detection", grp,
                r2, p2, n2, "rank residuals on detection rate")

        # ---- 2. class differences, raw and detection matched ---------------
        ln = V[V.group == "lnc"].dropna(subset=["tau_celltype"])
        groups = [ln[ln.cls == c].tau_celltype.values for c in CLASSES
                  if (ln.cls == c).sum() >= 20]
        if len(groups) >= 3:
            h, pk = kruskal(*groups)
            add("class", "Kruskal-Wallis tau_celltype across lncRNA classes",
                "lncRNA", h, pk, int(len(ln)),
                "raw, not detection controlled")
        pcg = V[(V.group == "pcg") & V.tau_celltype.notna()].copy()
        pcg["lg"] = np.log10(pcg.det_rate.clip(lower=1e-5))
        for c in CLASSES:
            sub = ln[ln.cls == c].copy()
            if len(sub) < 20:
                continue
            sub["lg"] = np.log10(sub.det_rate.clip(lower=1e-5))
            mm = match_one_to_one(sub.lg.values, pcg.lg.values, CALIPER)
            ok = mm >= 0
            if ok.sum() < 20:
                continue
            dtau = sub.tau_celltype.values[ok] - pcg.tau_celltype.values[mm[ok]]
            w, pw = mannwhitneyu(sub.tau_celltype.values[ok],
                                 pcg.tau_celltype.values[mm[ok]])
            add("class", "tau_celltype vs detection-matched PCG", c,
                float(np.median(dtau)), pw, int(ok.sum()),
                "median lncRNA minus matched PCG, caliper 0.1 log10")

        # ---- 3. evolution within lncRNAs -----------------------------------
        ln2 = ln.dropna(subset=["tau_celltype", "det_rate"])
        r, p, n2 = partial_spearman(ln2.tau_celltype, ln2.phylop_mean,
                                    ln2.det_rate)
        add("evolution", "partial Spearman tau_celltype vs phyloP", "lncRNA",
            r, p, n2, "controls detection")
        r, p, n2 = partial_spearman(ln2.tau_celltype, ln2.phast_mean,
                                    ln2.det_rate)
        add("evolution", "partial Spearman tau_celltype vs phastCons",
            "lncRNA", r, p, n2, "controls detection")
        r, p = spearmanr(ln2.tau_celltype, ln2.age_level, nan_policy="omit")
        add("evolution", "Spearman tau_celltype vs transcript age level",
            "lncRNA", r, p, int(len(ln2)), "0 human-specific to 3 eutherian")
        r, p, n2 = partial_spearman(ln2.tau_celltype, ln2.age_level,
                                    ln2.det_rate)
        add("evolution", "partial Spearman tau_celltype vs age level", "lncRNA",
            r, p, n2, "controls detection")
        # LTR-derived promoters are rarer, so the raw comparison is repeated
        # against non-LTR lncRNAs of matched detection.
        ltr = ln2[ln2.tss_in_ltr == 1].copy()
        non = ln2[ln2.tss_in_ltr == 0].copy()
        if len(ltr) > 20 and len(non) > 20:
            u, pu = mannwhitneyu(ltr.tau_celltype, non.tau_celltype)
            add("evolution", "MWU tau_celltype LTR-derived vs other", "lncRNA",
                float(np.median(ltr.tau_celltype) - np.median(non.tau_celltype)),
                pu, int(len(ltr) + len(non)), "median difference, raw")
            ltr["lg"] = np.log10(ltr.det_rate.clip(lower=1e-5))
            non["lg"] = np.log10(non.det_rate.clip(lower=1e-5))
            mm = match_one_to_one(ltr.lg.values, non.lg.values, CALIPER)
            ok = mm >= 0
            if ok.sum() > 20:
                dl = ltr.tau_celltype.values[ok] - non.tau_celltype.values[mm[ok]]
                u2, pu2 = mannwhitneyu(ltr.tau_celltype.values[ok],
                                      non.tau_celltype.values[mm[ok]])
                add("evolution", "MWU tau_celltype LTR vs non-LTR, matched",
                    "lncRNA", float(np.median(dl)), pu2, int(ok.sum()),
                    "detection-matched, caliper 0.1 log10")

        # ---- 4. disease: does cell-type specificity add to tissue breadth --
        dz = V[V.n_de >= 1].dropna(subset=["tau_celltype"])
        r, p = spearmanr(dz.tau_celltype, dz.n_up, nan_policy="omit")
        add("disease", "Spearman tau_celltype vs pan-cancer up count",
            "lncRNA with DE call", r, p, int(len(dz)), "raw")
        r2, p2, n2 = partial_spearman(dz.tau_celltype, dz.n_up,
                                      dz.gtex_breadth_count)
        add("disease", "partial Spearman controlling GTEx expression breadth",
            "lncRNA with DE call", r2, p2, n2, "rank residuals")
        r3, p3, n3 = partial_spearman(dz.tau_celltype, dz.n_up, dz.det_rate)
        add("disease", "partial Spearman controlling detection",
            "lncRNA with DE call", r3, p3, n3, "rank residuals")

        T = pd.DataFrame(rows)
        for fam, sub in T.groupby("family"):
            T.loc[sub.index, "fdr"] = bh_fdr(sub.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "cross_system_specificity_tests.tsv")), sep="\t",
            index=False)

        for r in T.itertuples():
            st.note("{} | {} | {} | stat={:.3f} n={} FDR={:.3g}".format(
                r.family, r.test, r.subset, r.statistic, r.n, r.fdr))
        st.note("cell-type tau medians: " + ", ".join(
            "{}={:.3f}".format(c, ln[ln.cls == c].tau_celltype.median())
            for c in CLASSES if (ln.cls == c).sum()))
        st.summary_md(os.path.join(
            RES, "step_summaries", "72_sc_tissue_disease_specificity.md"),
            "72 single-cell specificity against tissue and disease specificity",
            df_to_md(T, ".4f"))


if __name__ == "__main__":
    main()
