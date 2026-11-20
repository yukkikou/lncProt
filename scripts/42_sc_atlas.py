#!/usr/bin/env python
"""42_sc_atlas.py — public glioma single-cell atlas analysis (GBmap).

Modes:
  inspect  — h5ad structure (shape, obs/var keys, category values)
  analyze  — per cell-type / disease detection and expression of birth
             classes (5 lnc classes + PCG + LTR subset) in the atlas.
Outputs -> results/sc_atlas/*.tsv + step summary.
"""
import os
import sys

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import (binomtest, mannwhitneyu, rankdata, spearmanr,
                         wilcoxon)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
H5 = os.path.join(WORKDIR, "data", "scrna_atlas", "GBmap.h5ad")
OUT = os.path.join(RES, "sc_atlas")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
OBS_CANDIDATES = ["cell_type", "celltype", "cell_type_final", "annotation",
                  "cell_state", "celltype_major"]
DIS_CANDIDATES = ["disease", "condition", "sample_type", "tissue_type"]


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
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = [x.decode() if isinstance(x, bytes) else str(x)
                for x in g["categories"][:]]
        codes = g["codes"][:]
        return pd.Categorical.from_codes(codes, cats)
    vals = g[:]
    return pd.Series([x.decode() if isinstance(x, bytes) else str(x)
                      for x in vals])


def inspect():
    with step("42_sc_inspect") as st:
        f = h5py.File(H5, "r")
        print("shape", f["X"].attrs.get("shape"),
              "encoding", f["X"].attrs.get("encoding-type"))
        print("obs keys", list(f["obs"].keys()))
        print("var keys", list(f["var"].keys()))
        for k in f["obs"].keys():
            if k.startswith("_"):
                continue
            try:
                s = read_cat(f, k)
                if s.nunique() <= 40:
                    print("OBS", k, "nunique", s.nunique(),
                          dict(pd.Series(s).value_counts().head(15)))
            except Exception as e:
                pass
        vi = f["var"]["_index"][:5]
        print("var index head", [x.decode() if isinstance(x, bytes)
                                 else str(x) for x in vi])
        f.close()
        st.note("inspect done")



def mwu(y, x):
    u, p = mannwhitneyu(y, x, alternative="two-sided")
    return float(2 * u / (len(y) * len(x)) - 1), float(p)


def analyze():
    with step("42_sc_analyze") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        base = pd.Series(meta.index.str.split(".").str[0].values,
                         index=meta.index)
        m2 = meta.copy()
        m2["base"] = base.values
        m2 = m2.drop_duplicates("base")
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["tss_in_ltr"]]
        te.index = te.index.str.split(".").str[0]
        te = te[~te.index.duplicated()]
        m2 = m2.join(te, on="base")

        f = h5py.File(H5, "r")
        var = [x.decode() if isinstance(x, bytes) else str(x)
               for x in f["var"]["_index"][:]]
        cls = m2.set_index("base")["cls"].reindex(var)
        grp = m2.set_index("base")["group"].reindex(var)
        ltr = m2.set_index("base")["tss_in_ltr"].reindex(var)
        sets = {}
        for c in CLASSES:
            sets[c] = np.where((cls == c).values)[0]
        sets["PCG"] = np.where((grp == "pcg").values)[0]
        sets["pseudogene"] = np.where((grp != "pcg") &
                                      (cls == "pseudogene").values)[0]
        lnc_mask = cls.isin(CLASSES).values
        sets["LTR_lnc"] = np.where(lnc_mask & (ltr == 1).values)[0]
        sets["other_lnc"] = np.where(lnc_mask & ~(ltr == 1).values)[0]
        st.note("atlas genes: {}; class sizes: {}".format(
            len(var), {k: len(v) for k, v in sets.items()}))

        n = int(f["X"].attrs["shape"][0])
        nv = int(f["X"].attrs["shape"][1])
        indptr = f["X"]["indptr"][:]
        data = f["X"]["data"]
        indices = f["X"]["indices"]
        obs = {}
        for k in ["cell_type", "annotation_level_1", "donor_id", "stage",
                  "tissue", "disease"]:
            try:
                obs[k] = np.asarray(read_cat(f, k).astype(str))
            except Exception:
                obs[k] = np.array(["NA"] * n)
        cols = {}
        CH = 10000
        for s in range(0, n, CH):
            e = min(s + CH, n)
            lo, hi = int(indptr[s]), int(indptr[e])
            X = csr_matrix((data[lo:hi], indices[lo:hi],
                            indptr[s:e + 1] - lo), shape=(e - s, nv))
            overall = X.getnnz(axis=1) / nv
            cols.setdefault("overall_det", []).extend(overall.tolist())
            for name, idx in sets.items():
                if len(idx) == 0:
                    continue
                sub = X[:, idx]
                det = sub.getnnz(axis=1) / len(idx)
                mean = np.asarray(sub.mean(axis=1)).ravel()
                cols.setdefault(name + "_det", []).extend(det.tolist())
                cols.setdefault(name + "_adj", []).extend(
                    (det - overall).tolist())
                cols.setdefault(name + "_mean", []).extend(mean.tolist())
        f.close()
        df = pd.DataFrame(cols)
        for k in obs:
            df[k] = obs[k]
        df.to_csv(safe_write_path(os.path.join(
            OUT, "per_cell_scores.tsv.gz")), sep="\t", index=False,
            compression="gzip")

        # cell-type table
        det_cols = [c + "_det" for c in sets] + ["overall_det"]
        adj_cols = [c + "_adj" for c in sets]
        ct = df.groupby("cell_type")[det_cols + adj_cols].mean()
        ct.to_csv(safe_write_path(os.path.join(
            OUT, "class_by_celltype.tsv")), sep="\t")

        # donor-level means by neoplastic status
        dm = df.groupby(["donor_id", "annotation_level_1"])[
            adj_cols + [c + "_mean" for c in sets]].mean().reset_index()
        rows = []
        for metric in adj_cols:
            a = dm.loc[dm.annotation_level_1 == "Neoplastic", metric].dropna()
            b = dm.loc[dm.annotation_level_1 == "Non-neoplastic",
                       metric].dropna()
            if len(a) >= 5 and len(b) >= 5:
                d, p = mwu(a.values, b.values)
                rows.append(dict(metric=metric, n_neoplastic=len(a),
                                 n_nonneoplastic=len(b),
                                 mean_neoplastic=float(a.mean()),
                                 mean_non=float(b.mean()),
                                 cliffs_delta=d, pval=p))
        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values) if len(T) else []
        T.to_csv(safe_write_path(os.path.join(
            OUT, "malignant_vs_non.tsv")), sep="\t", index=False)

        # stage within neoplastic donors
        stg = df[df.annotation_level_1 == "Neoplastic"].groupby(
            ["donor_id", "stage"])[adj_cols].mean().reset_index()
        rows = []
        for metric in adj_cols:
            a = stg.loc[stg.stage == "Primary", metric].dropna()
            b = stg.loc[stg.stage == "Recurrent", metric].dropna()
            if len(a) >= 5 and len(b) >= 5:
                d, p = mwu(a.values, b.values)
                rows.append(dict(metric=metric, n_primary=len(a),
                                 n_recurrent=len(b),
                                 mean_primary=float(a.mean()),
                                 mean_recurrent=float(b.mean()),
                                 cliffs_delta=d, pval=p))
        S = pd.DataFrame(rows)
        S["fdr"] = bh_fdr(S.pval.values) if len(S) else []
        S.to_csv(safe_write_path(os.path.join(
            OUT, "stage_tests.tsv")), sep="\t", index=False)

        st.note("malignant vs non: " + "; ".join(
            "{} d={:.3f} p={:.1e}".format(r.metric, r.cliffs_delta, r.pval)
            for r in T.itertuples() if r.metric in
            ("divergent_adj", "intergenic_adj", "LTR_lnc_adj")))
        st.summary_md("results/sc_atlas/step_summaries/42_sc.md",
                      "42 GBmap single-cell atlas",
                      "## Class detection by cell type\n\n"
                      + df_to_md(ct, ".4f", max_rows=40)
                      + "\n\n## Malignant vs non-neoplastic (donor-level)\n\n"
                      + df_to_md(T, ".3f")
                      + "\n\n## Primary vs recurrent (neoplastic)\n\n"
                      + df_to_md(S, ".3f"))


PAIR_MAX = 400


def load_per_cell():
    return pd.read_csv(os.path.join(OUT, "per_cell_scores.tsv.gz"), sep="\t")


def paired_donors():
    """Within-donor contrast of malignant versus non-malignant cells.

    malignant_vs_non.tsv contrasts 100 neoplastic with 108 non-neoplastic
    donors. Most donors contribute both cell states, so a paired test removes
    between-patient variation from the comparison.
    """
    with step("42_sc_paired") as st:
        d = load_per_cell()
        mets = [c for c in d.columns if c.endswith("_adj")]
        dm = d.groupby(["donor_id", "annotation_level_1"])[mets].mean()
        wide = dm.unstack("annotation_level_1")
        rows = []
        for m in mets:
            if (m, "Neoplastic") not in wide.columns:
                continue
            pair = pd.concat([wide[(m, "Neoplastic")],
                              wide[(m, "Non-neoplastic")]],
                             axis=1, keys=["neo", "non"]).dropna()
            if len(pair) < 5:
                continue
            delta = (pair.neo - pair.non).values
            p = float(wilcoxon(delta).pvalue)
            rng = np.random.default_rng(42)
            boot = np.array([np.median(rng.choice(delta, len(delta),
                                                  replace=True))
                             for _ in range(1000)])
            rows.append(dict(metric=m, n_pairs=len(delta),
                             median_delta=float(np.median(delta)),
                             ci_low=float(np.percentile(boot, 2.5)),
                             ci_high=float(np.percentile(boot, 97.5)),
                             frac_positive=float((delta > 0).mean()),
                             pval=p))
        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values) if len(T) else []
        T.to_csv(safe_write_path(os.path.join(
            OUT, "malignant_paired.tsv")), sep="\t", index=False)
        st.note("paired donors n={}; ".format(int(T.n_pairs.max()))
                + "; ".join("{} d={:+.4f} p={:.1e}".format(
                    r.metric, r.median_delta, r.pval)
                    for r in T.itertuples()
                    if r.metric in ("PCG_adj", "divergent_adj",
                                    "intergenic_adj")))
        st.summary_md(os.path.join(
            OUT, "step_summaries/42_sc_paired.md"),
            "42 paired malignant versus non-malignant",
            df_to_md(T, ".4f"))


def hierarchy():
    """Rank consistency of the class-detection hierarchy across cell types.

    Detection rates differ by cell type simply because cells differ in total
    RNA content. What matters for the paper is whether the ordering of the six
    classes is the same everywhere, which a rank correlation answers.
    """
    with step("42_sc_hierarchy") as st:
        ct = pd.read_csv(os.path.join(OUT, "class_by_celltype.tsv"), sep="\t")
        mets = ["PCG_det", "divergent_det", "antisense_det",
                "intronic_det", "sense_overlapping_det", "intergenic_det"]
        m = ct[mets].to_numpy(float)
        ranks = np.apply_along_axis(rankdata, 1, m)
        r = np.corrcoef(ranks)
        iu = np.triu_indices_from(r, 1)
        off = r[iu]
        out = ct[["cell_type"]].copy()
        for j, c in enumerate(mets):
            out[c.replace("_det", "_rank")] = rankdata(
                -m[:, j], method="min").astype(int)
        out.to_csv(safe_write_path(os.path.join(
            OUT, "celltype_hierarchy.tsv")), sep="\t", index=False)
        st.note("cell types={}; rank Spearman min/median/max = "
                "{:.3f}/{:.3f}/{:.3f}".format(
                    len(ct), off.min(), np.median(off), off.max()))
        st.summary_md(os.path.join(
            OUT, "step_summaries/42_sc_hierarchy.md"),
            "42 cell-type rank consistency",
            "Pairwise Spearman over {} cell types: min {:.3f}, median {:.3f}, "
            "max {:.3f}.\n\n".format(len(ct), off.min(), np.median(off),
                                     off.max())
            + df_to_md(out, "4"))


def pairs():
    """Divergent-pair co-detection in the atlas, with matched-control pairs.

    Repeats the 765-cell analysis of module 23 at 338,564 cells. Detection
    rates for every gene are counted in a first pass, then only the pair genes
    and their detection-matched controls are extracted in a second pass so the
    dense matrix stays small.
    """
    with step("42_sc_pairs") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        base = pd.Series(meta.index.str.split(".").str[0].values,
                         index=meta.index)
        m2 = meta.copy()
        m2["base"] = base.values
        m2 = m2.drop_duplicates("base")
        f = h5py.File(H5, "r")
        var = np.array([x.decode() if isinstance(x, bytes) else str(x)
                        for x in f["var"]["_index"][:]])
        vidx = {g: i for i, g in enumerate(var)}
        cls = m2.set_index("base")["cls"].reindex(var)
        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        cand = []
        for r in dp.itertuples():
            a = vidx.get(str(r.lnc_gid).split(".")[0])
            b = vidx.get(str(r.pcg_gid).split(".")[0])
            if a is not None and b is not None:
                cand.append((a, b))
        n = int(f["X"].attrs["shape"][0])
        nv = int(f["X"].attrs["shape"][1])
        indptr = f["X"]["indptr"][:]
        indices = f["X"]["indices"]
        counts = np.zeros(nv, dtype=np.int64)
        CH = 20000
        for s in range(0, n, CH):
            e = min(s + CH, n)
            lo, hi = int(indptr[s]), int(indptr[e])
            counts += np.bincount(indices[lo:hi], minlength=nv)
        rate = counts / n
        st.note("pairs in matrix={}; atlas genes={}".format(len(cand), nv))

        lnc_idx = np.where(cls.isin(CLASSES).values)[0]
        lnc_rate = rate[lnc_idx]
        order = np.argsort(lnc_rate)
        pcg_idx = np.where((m2.set_index("base")["group"]
                            .reindex(var) == "pcg").values)[0]
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

        gm = m2.set_index("base").reindex(var)
        pcg_mask = (gm["group"] == "pcg").values
        chrom_arr = gm["chrom"].astype(str).values
        tss_arr = pd.to_numeric(gm["tss"], errors="coerce").values
        strand_arr = gm["strand"].astype(str).values
        by_chr = {}
        for c in np.unique(chrom_arr[pcg_mask]):
            ids = np.where(pcg_mask & (chrom_arr == c))[0]
            by_chr[c] = ids[np.argsort(tss_arr[ids])]

        def neighbour(a, b, d):
            """Distance-matched PCG on the partner strand near the lncRNA TSS.

            This mirrors the bulk CAGE control, where a nearby gene rather than
            an arbitrary gene of matching expression is the null.
            """
            ids = by_chr.get(chrom_arr[a])
            if ids is None or not np.isfinite(tss_arr[a]):
                return None
            dist = np.abs(tss_arr[ids] - tss_arr[a])
            ok = (ids != b) & (strand_arr[ids] == strand_arr[b]) & (dist <= 1e5)
            cands = ids[ok]
            dd = dist[ok]
            if len(cands) == 0:
                return None
            rr = np.abs(np.log2(rate[cands] / max(rate[b], 1e-6)))
            keep = rr < 1.0
            if keep.any():
                cands, dd = cands[keep], dd[keep]
            return int(cands[np.argmin(np.abs(dd - d))])

        dist_lookup = {}
        for r in dp.itertuples():
            a = vidx.get(str(r.lnc_gid).split(".")[0])
            b = vidx.get(str(r.pcg_gid).split(".")[0])
            if a is not None and b is not None:
                dist_lookup[(a, b)] = r.distance

        test = []
        for a, b in cand[:PAIR_MAX]:
            if rate[a] < 0.005 or rate[b] < 0.005:
                continue
            c1 = matched(rate[a], {a, b}, lnc_idx, lnc_rate, order)
            c2 = matched(rate[b], {a, b, c1}, lnc_idx, lnc_rate, order)
            # second control keeps the pair composition (lnc plus PCG) and
            # replaces only the partner, so a shared cell-type programme cannot
            # masquerade as pair coupling.
            c3 = matched(rate[b], {a, b, c1, c2}, pcg_idx, pcg_rate, porder)
            d = dist_lookup.get((a, b), -1)
            c4 = neighbour(a, b, float(d)) if d > 0 else None
            test.append((a, b, c1, c2, c3, c4))
        st.note("pairs testable={}".format(len(test)))
        cols = np.array(sorted({x for t in test for x in t if x is not None}),
                        dtype=np.int64)
        local = {int(c): i for i, c in enumerate(cols)}
        ia = np.array([local[t[0]] for t in test])
        ib = np.array([local[t[1]] for t in test])
        ic = np.array([local[t[2]] for t in test])
        id_ = np.array([local[t[3]] for t in test])
        ie = np.array([local[t[4]] for t in test])
        has4 = np.array([t[5] is not None for t in test])
        i4 = np.array([local[t[5]] if t[5] is not None else local[t[0]]
                       for t in test])
        k_both = np.zeros(len(test))
        k_ctrl = np.zeros(len(test))
        k_ctrlp = np.zeros(len(test))
        k_neigh = np.zeros(len(test))
        for s in range(0, n, CH):
            e = min(s + CH, n)
            lo, hi = int(indptr[s]), int(indptr[e])
            idx = indices[lo:hi]
            rows_e = np.repeat(np.arange(e - s, dtype=np.int32),
                               np.diff(indptr[s:e + 1]))
            keep = np.isin(idx, cols)
            if not keep.any():
                continue
            B = np.zeros((e - s, len(cols)), dtype=bool)
            B[rows_e[keep], np.searchsorted(cols, idx[keep])] = True
            k_both += (B[:, ia] & B[:, ib]).sum(axis=0)
            k_ctrl += (B[:, ic] & B[:, id_]).sum(axis=0)
            k_ctrlp += (B[:, ia] & B[:, ie]).sum(axis=0)
            k_neigh[has4] += (B[:, ia[has4]] & B[:, i4[has4]]).sum(axis=0)
        f.close()

        rows = []
        for j, t in enumerate(test):
            k, kc, kp = k_both[j], k_ctrl[j], k_ctrlp[j]
            pa, pb = rate[t[0]], rate[t[1]]
            exp = pa * pb * n
            pc = rate[t[2]] * rate[t[3]] * n
            pp = rate[t[0]] * rate[t[4]] * n
            if k < 5:
                continue
            p = binomtest(int(k), n, pa * pb, alternative="greater").pvalue
            row = dict(pair="{0}|{1}".format(var[t[0]], var[t[1]]),
                       distance=int(dist_lookup.get((t[0], t[1]), -1)),
                       k_both=int(k), expected=exp,
                       ratio=k / max(exp, 1e-9), p_binom=float(p),
                       control_k=int(kc), control_expected=pc,
                       control_ratio=kc / max(pc, 1e-9),
                       pcg_control_k=int(kp), pcg_control_expected=pp,
                       pcg_control_ratio=kp / max(pp, 1e-9))
            if has4[j]:
                kn = k_neigh[j]
                pn = rate[t[0]] * rate[t[5]] * n
                row.update(neighbour_k=int(kn), neighbour_expected=pn,
                           neighbour_ratio=kn / max(pn, 1e-9))
            rows.append(row)
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "pair_codetection_atlas.tsv")), sep="\t", index=False)
        if len(T):
            above = float((T.ratio > T.control_ratio).mean())
            abovep = float((T.ratio > T.pcg_control_ratio).mean())
            if "neighbour_ratio" in T:
                nb = T.dropna(subset=["neighbour_ratio"])
                aboven = float((nb.ratio > nb.neighbour_ratio).mean())
                nbmed = float(nb.neighbour_ratio.median())
            else:
                aboven, nbmed = float("nan"), float("nan")
            rho_d = spearmanr(T.ratio, T.distance).statistic
            st.note("pairs reported={}; median obs/exp={:.3f}; lnc control="
                    "{:.3f} (above {:.2f}); PCG control={:.3f} (above {:.2f});"
                    " neighbour control={:.3f} (above {:.2f}, n={}); "
                    "rho(ratio, distance)={:.3f}".format(
                        len(T), T.ratio.median(), T.control_ratio.median(),
                        above, T.pcg_control_ratio.median(), abovep, nbmed,
                        aboven, len(T) - T.neighbour_ratio.isna().sum()
                        if "neighbour_ratio" in T else 0, rho_d))
            st.summary_md(os.path.join(
                OUT, "step_summaries/42_sc_pairs.md"),
                "42 divergent-pair co-detection (atlas)",
                "Pairs reported {}. Median observed/expected ratio {:.3f}. "
                "Lnc-matched control {:.3f}, fraction of pairs above it "
                "{:.2f}. PCG-matched control {:.3f}, fraction above it "
                "{:.2f}. Distance-matched neighbour control {:.3f}, fraction "
                "above it {:.2f}. Spearman of ratio with intergenic distance "
                "{:.3f}.\n\n".format(
                    len(T), T.ratio.median(), T.control_ratio.median(),
                    above, T.pcg_control_ratio.median(), abovep, nbmed,
                    aboven, rho_d)
                + df_to_md(T.head(25), ".3f"))


def main():
    mode = (sys.argv[1] if len(sys.argv) > 1 else "inspect").lower()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    if mode == "inspect":
        inspect()
    elif mode == "analyze":
        analyze()
    elif mode == "paired":
        paired_donors()
    elif mode == "hierarchy":
        hierarchy()
    elif mode == "pairs":
        pairs()
    else:
        raise SystemExit(
            "usage: 42_sc_atlas.py inspect|analyze|paired|hierarchy|pairs")


if __name__ == "__main__":
    main()
