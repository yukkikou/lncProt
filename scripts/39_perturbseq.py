#!/usr/bin/env python
"""39_perturbseq.py — Replogle et al. 2022 genome-scale Perturb-seq.

Class-level response of lncRNA birth classes / LTR-derived lncRNAs to
perturbation of shared TFs and cofactors (CTCF, MAZ, E2F6, EP300).

Data (Figshare, CC BY 4.0):
  processed Perturb-seq AnnData (article 20029387):
    $pop_{raw,normalized}_{bulk,singlecell}_01.h5ad
  We use the *pseudobulk* files only (374 MB for K562_gwps); the single-cell
  files are 8.7-65.8 GB each and are never downloaded by this script.

Modes
  inspect  — file inventory, obs/var schema, perturbation/control labels,
             cell counts and target availability (no expression loaded).
  analyze  — pseudobulk log2 CPM vs non-targeting controls -> per-gene
             response -> class-level medians, within-target gene-set tests,
             and target-vs-background empirical ranks over all perturbations.

Examples
  python -I scripts/39_perturbseq.py inspect
  python -I scripts/39_perturbseq.py analyze --dataset K562_gwps
  python -I scripts/39_perturbseq.py analyze --file data/perturb/replogle2022/K562_gwps_raw_bulk_01.h5ad

Outputs -> results/perturbseq/
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
PD = os.path.join(WORKDIR, "data", "perturb", "replogle2022")
OUT = os.path.join(RES, "perturbseq")
GMT = os.path.join(RES, "gene_master_table.tsv")
TE = os.path.join(RES, "te_promoter_stats.tsv")

TARGETS = ["CTCF", "MAZ", "E2F6", "EP300"]
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
CONTROL_RE = re.compile(
    r"^(NTC|NT|NON[-_ ]?TARGETING|CONTROL|CTRL|SCRAMBLE[D]?)([_\-\s].*)?$",
    re.I)
# controls may be embedded in a composite label, e.g.
# "10747_non-targeting_non-targeting_non-targeting"
CTRL_SEARCH = re.compile(
    r"(^|[_\-\s])(NTC|NT|NON[-_ ]?TARGETING|CONTROL|CTRL|SCRAMBLE[D]?)"
    r"([_\-\s]|$)", re.I)
# Replogle pseudobulk obs index: "<n>_<SYMBOL>_<P1P2|ENST...>_<ENSG...>"
REPLOGLE_IDX = re.compile(r"^\d+_([^_]+)_")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def safe_median(v):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else np.nan


def target_re(t):
    return re.compile(r"^" + re.escape(t) + r"([_\-\s].*)?$", re.I)


def guess_pert_col(obs):
    """Pick the obs column holding the perturbation / target label."""
    named = ["perturbation", "target", "gene", "perturbation_name",
             "condition", "guide_target", "target_gene"]
    for c in named:
        if c in obs.columns:
            return c
    # fall back: string/categorical column with a control-like value
    best, best_hits = None, -1
    for c in obs.columns:
        s = obs[c]
        if s.dtype.kind not in "OUS" and str(s.dtype) != "category":
            continue
        vals = s.astype(str)
        hits = int(vals.map(lambda x: bool(CONTROL_RE.match(x))).sum())
        if hits > best_hits:
            best, best_hits = c, hits
    return best


def is_control(vals):
    return np.asarray([bool(CTRL_SEARCH.search(str(v))) for v in vals])


def parse_pert_label(s):
    """Extract the target symbol from a Replogle composite obs-index label."""
    s = str(s)
    if CTRL_SEARCH.search(s):
        return "NTC"
    m = REPLOGLE_IDX.match(s)
    if m:
        return m.group(1)
    m2 = re.search(r"(ENSG\d+)", s)
    return m2.group(1) if m2 else s


def get_pert_labels(obs, pert_col=None):
    """Return parsed perturbation labels, or None if none can be found."""
    if pert_col:
        return pd.Series(obs[pert_col].astype(str).values, index=obs.index)
    idx = pd.Index(obs.index.astype(str))
    if len(idx) == 0:
        return None
    frac = np.mean([bool(REPLOGLE_IDX.match(x)) or
                    bool(CTRL_SEARCH.search(x)) for x in idx])
    if frac > 0.5:
        return pd.Series([parse_pert_label(x) for x in idx], index=obs.index)
    return None


def n_cells_col(obs):
    for c in ["n_cells", "ncells", "n_cell", "cells", "cell_count",
              "n_cells_perturbation", "num_cells_filtered",
              "num_cells_unfiltered"]:
        if c in obs.columns:
            return c
    return None


def gene_identifiers(adata):
    """Return (ids, kind) where kind is 'ensembl' or 'symbol'."""
    idx = pd.Index(adata.var_names.astype(str))
    ens_frac = idx.str.upper().str.startswith("ENSG").mean()
    if ens_frac > 0.5:
        return idx, "ensembl"
    for c in ["gene_ids", "gene_id", "ensembl_id", "gene_symbols", "symbol",
              "gene_name", "gene"]:
        if c in adata.var.columns:
            v = pd.Index(adata.var[c].astype(str))
            if v.str.upper().str.startswith("ENSG").mean() > 0.5:
                return v, "ensembl"
            if c in ("gene_symbols", "symbol", "gene_name", "gene"):
                return v, "symbol"
    return idx, "symbol"


def load_annotation():
    """gene_master_table + te_promoter_stats -> annotation indexed two ways."""
    gmt = pd.read_csv(GMT, sep="\t", low_memory=False)
    te = pd.read_csv(TE, sep="\t", usecols=["gene_id", "tss_in_ltr",
                                            "te_frac_LTR"], low_memory=False)
    ann = gmt.merge(te, on="gene_id", how="left", suffixes=("", "_te"))
    phy_path = os.path.join(RES, "phylop", "phylop_per_gene.tsv")
    if os.path.exists(phy_path):
        phy = pd.read_csv(phy_path, sep="\t",
                          usecols=["gene_id", "phylop_mean", "phast_mean",
                                   "gc"], low_memory=False)
        ann = ann.merge(phy, on="gene_id", how="left")
    ann["base"] = ann["gene_id"].astype(str).str.split(".").str[0]
    ann["tss_in_ltr"] = pd.to_numeric(ann["tss_in_ltr"], errors="coerce")
    ann["is_lnc"] = ann["group"].eq("lnc") & ann["cls"].isin(CLASSES)
    ann["is_ltr"] = ann["tss_in_ltr"].fillna(0).eq(1)
    ann["group3"] = np.where(
        ann["is_lnc"] & ann["is_ltr"], "LTR_lnc",
        np.where(ann["is_lnc"], "lnc_other",
                 np.where(ann["group"].eq("pcg"), "PCG", "other")))
    by_ens = ann.drop_duplicates("base").set_index("base", drop=False)
    by_sym = ann.drop_duplicates("gene_name").set_index("gene_name", drop=False)
    return ann, by_ens, by_sym


def map_genes(ids, kind, by_ens, by_sym):
    """Map a mixed vector of Ensembl IDs and gene symbols, one by one."""
    ids = pd.Index(ids).astype(str)
    ens = np.asarray(ids.str.upper().str.startswith("ENSG"), dtype=bool)
    frames = []
    if ens.any():
        e = by_ens.reindex(ids[ens].str.split(".").str[0]).copy()
        e.index = ids[ens]
        frames.append(e)
    if (~ens).any():
        s = by_sym.reindex(ids[~ens]).copy()
        s.index = ids[~ens]
        frames.append(s)
    if frames:
        return pd.concat(frames).reindex(ids)
    return pd.DataFrame(index=ids, columns=by_ens.columns)


def class_masks(mapped):
    cls = mapped["cls"].astype(str).values
    grp = mapped["group"].astype(str).values
    ltr = mapped["tss_in_ltr"].fillna(0).values.astype(float) == 1
    lnc = (grp == "lnc")
    masks = {}
    for c in CLASSES:
        masks[c] = lnc & (cls == c)
    masks["LTR_lnc"] = lnc & ltr
    masks["lnc_other"] = lnc & ~ltr
    masks["all_lnc"] = lnc
    masks["PCG"] = grp == "pcg"
    masks["other_lnc"] = lnc & (cls != "divergent")
    masks["expressed_all"] = np.ones(len(cls), dtype=bool)
    return masks


def partner_pairs(mapped):
    """Index pairs (divergent lncRNA, its partner PCG) present in the matrix."""
    base_to_pos = {}
    for pos, b in enumerate(mapped["base"].values):
        if isinstance(b, str) and b and b not in base_to_pos:
            base_to_pos[b] = pos
    li, pi = [], []
    div = mapped[(mapped["cls"].astype(str) == "divergent")
                 & mapped["partner_pcg"].notna()
                 & mapped["partner_pcg"].astype(str).ne("")]
    for _, r in div.iterrows():
        b = str(r["base"])
        p = str(r["partner_pcg"]).split(".")[0]
        if b in base_to_pos and p in base_to_pos:
            li.append(base_to_pos[b])
            pi.append(base_to_pos[p])
    return np.asarray(li, int), np.asarray(pi, int)


def partner_stats(resp, det, li, pi):
    if li.size < 20:
        return dict(partner_n=0, partner_sign_conc=np.nan,
                    partner_pearson=np.nan)
    keep = det[li] & det[pi]
    if int(keep.sum()) < 20:
        return dict(partner_n=int(keep.sum()), partner_sign_conc=np.nan,
                    partner_pearson=np.nan)
    a = resp[li[keep]]
    b = resp[pi[keep]]
    same = float(np.mean((a > 0) & (b > 0)) + np.mean((a < 0) & (b < 0)))
    if np.std(a) > 0 and np.std(b) > 0:
        pear = float(np.corrcoef(a, b)[0, 1])
    else:
        pear = np.nan
    return dict(partner_n=int(keep.sum()), partner_sign_conc=same,
                partner_pearson=pear)


def response_metrics(resp, det, masks, bin_id, n_bins):
    """Class-level summaries for one perturbation response vector."""
    out = {}
    det_resp = resp[det]
    out["global_median"] = safe_median(det_resp)
    out["global_mean_abs"] = (float(np.mean(np.abs(det_resp)))
                              if det_resp.size else np.nan)
    # expression-matched null: median response per control-expression bin
    binmed = np.full(n_bins, np.nan)
    for b in range(n_bins):
        v = resp[det & (bin_id == b)]
        if v.size >= 5:
            binmed[b] = np.median(v)
    for name, m in masks.items():
        mm = m & det
        n = int(mm.sum())
        if n >= 5:
            v = resp[mm]
            valid = bin_id[mm] >= 0
            resid = (v[valid] - binmed[bin_id[mm][valid]]
                     if valid.any() else np.array([]))
            out[name] = dict(
                n=n, median=safe_median(v), mean=float(np.mean(v)),
                mean_abs=float(np.mean(np.abs(v))),
                frac_up=float(np.mean(v > 1)), frac_down=float(np.mean(v < -1)),
                matched=safe_median(resid))
        else:
            out[name] = dict(n=n, median=np.nan, mean=np.nan, mean_abs=np.nan,
                             frac_up=np.nan, frac_down=np.nan, matched=np.nan)
    out["ltr_minus_other"] = (out["LTR_lnc"]["median"]
                              - out["lnc_other"]["median"])
    out["divergent_minus_otherlnc"] = (out["divergent"]["median"]
                                       - out["other_lnc"]["median"])
    return out


ALL_GROUPS = CLASSES + ["LTR_lnc", "lnc_other", "PCG", "all_lnc",
                        "other_lnc"]


def flatten_metrics(met):
    """Turn nested per-class dicts into scalar columns (for DataFrames)."""
    out = {k: v for k, v in met.items() if not isinstance(v, dict)}
    for c in ALL_GROUPS:
        d = met.get(c, {})
        for m in ["n", "median", "mean", "mean_abs", "frac_up", "frac_down",
                  "matched"]:
            out[f"{c}_{m}"] = d.get(m, np.nan)
    return out


def mwu_row(a, b, label, target):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 20 or b.size < 20:
        return None
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return dict(target=target, contrast=label, n_case=a.size, n_ctrl=b.size,
                median_case=float(np.median(a)),
                median_ctrl=float(np.median(b)),
                cliffs=float(2 * u / (a.size * b.size) - 1), pval=float(p))


# --------------------------------------------------------------------------- #
# inspect
# --------------------------------------------------------------------------- #
def cmd_inspect(args):
    import anndata as ad
    files = sorted(glob.glob(os.path.join(PD, "*.h5ad")))
    os.makedirs(OUT, exist_ok=True)
    inv, avail = [], []
    with step("39_perturbseq_inspect") as st:
        if not files:
            st.note("no h5ad files found in " + PD)
            return
        for f in files:
            base = os.path.basename(f)
            size = os.path.getsize(f)
            try:
                a = ad.read_h5ad(f, backed="r")
            except Exception as e:  # corrupt / partial download
                inv.append(dict(file=base, size_bytes=size, status="ERROR",
                                error=str(e)[:200]))
                st.note(f"{base}: unreadable ({e})")
                continue
            obs = a.obs
            pert_col = guess_pert_col(obs)
            labels = get_pert_labels(obs, pert_col)
            if labels is None:
                labels = pd.Series(["?"] * a.n_obs, index=obs.index)
                pert_disp = "none"
            else:
                pert_disp = (str(pert_col) if pert_col
                             else "index:gene_transcript(parsed)")
            vals = labels
            ctrl = is_control(vals)
            nc = n_cells_col(obs)
            cells_total = (pd.to_numeric(obs[nc].astype(str),
                                         errors="coerce").sum()
                           if nc else a.n_obs)
            inv.append(dict(
                file=base, size_bytes=size, status="ok",
                n_obs=int(a.n_obs), n_vars=int(a.n_vars),
                x_type=type(a.X).__name__, pert_col=pert_disp,
                n_perturbations=int(vals[~ctrl].nunique()),
                n_control_labels=int(vals[ctrl].nunique()),
                control_labels=";".join(sorted(vals[ctrl].unique())[:6]),
                cells_total=float(cells_total) if nc else int(a.n_obs),
                n_cells_col=str(nc) if nc else "",
                obs_cols=";".join(map(str, obs.columns))[:300],
                var_cols=";".join(map(str, a.var.columns))[:300]))
            for t in TARGETS:
                m = vals.map(lambda x: bool(target_re(t).match(x)))
                if m.any():
                    rows = int(m.sum())
                    cells = (float(pd.to_numeric(obs.loc[m, nc].astype(str),
                                                 errors="coerce").sum())
                             if nc else rows)
                    labels = ";".join(sorted(vals[m].unique())[:4])
                    avail.append(dict(file=base, target=t, present=1,
                                      n_rows=rows, n_cells=cells,
                                      labels=labels))
                else:
                    avail.append(dict(file=base, target=t, present=0,
                                      n_rows=0, n_cells=0, labels=""))
            st.note(f"{base}: {a.n_obs} obs x {a.n_vars} vars, "
                    f"pert_col={pert_col}, ctrl={sorted(vals[ctrl].unique())[:3]}")
            del a
        I = pd.DataFrame(inv)
        A = pd.DataFrame(avail)
        I.to_csv(safe_write_path(os.path.join(OUT, "inventory.tsv")),
                 sep="\t", index=False)
        if len(A):
            A.to_csv(safe_write_path(
                os.path.join(OUT, "target_availability_all.tsv")),
                sep="\t", index=False)
        st.summary_md("results/perturbseq/step_summaries/39_inspect.md",
                      "39 Perturb-seq file inventory",
                      "## Files\n\n" + df_to_md(I, max_rows=30)
                      + "\n\n## Target availability\n\n"
                      + df_to_md(A, max_rows=60))
        st.note("inventory -> " + os.path.join(OUT, "inventory.tsv"))
        if len(A):
            st.note("targets -> " + os.path.join(
                OUT, "target_availability_all.tsv"))


# --------------------------------------------------------------------------- #
# analyze
# --------------------------------------------------------------------------- #
def pick_file(dataset):
    cands = sorted(glob.glob(os.path.join(PD, "*_bulk_*.h5ad")))
    if dataset:
        pref = [c for c in cands if dataset.lower() in
                os.path.basename(c).lower()]
        cands = pref
    raw = [c for c in cands if "raw_bulk" in os.path.basename(c)]
    norm = [c for c in cands if "normalized_bulk" in os.path.basename(c)]
    if raw:
        return raw[0]
    if norm:
        return norm[0]
    return None


def cmd_analyze(args):
    import anndata as ad
    f = args.file or pick_file(args.dataset)
    if not f or not os.path.exists(f):
        raise SystemExit("no pseudobulk h5ad found; run inspect first or pass "
                         "--file")
    os.makedirs(OUT, exist_ok=True)
    with step("39_perturbseq_analyze") as st:
        st.note("file=" + os.path.basename(f))
        a = ad.read_h5ad(f)
        obs = a.obs
        pert_col = args.pert_col or guess_pert_col(obs)
        labels = get_pert_labels(obs, pert_col)
        if labels is None:
            raise SystemExit("cannot identify perturbation labels")
        st.note("perturbation labels: " + (str(pert_col) if pert_col
                                           else "obs index parsed"))
        vals = labels.astype(str).values
        ctrl = is_control(vals)
        if not ctrl.any():
            raise SystemExit("no non-targeting control rows found")
        nc = n_cells_col(obs)
        row_cells = (pd.to_numeric(obs[nc].astype(str),
                                   errors="coerce").fillna(1).values
                     if nc else np.ones(a.n_obs))

        # expression matrix -> log2 CPM
        X = a.X
        if sparse.issparse(X):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)
        np.nan_to_num(X, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        lib = X.sum(axis=1, dtype=np.float64)
        logcpm = np.log2(X / np.maximum(lib, 1.0)[:, None] * 1e6 + 1.0)
        logcpm = logcpm.astype(np.float32)
        del X

        ids, kind = gene_identifiers(a)
        ann, by_ens, by_sym = load_annotation()
        mapped = map_genes(ids, kind, by_ens, by_sym)
        n_mapped = int(mapped["cls"].notna().sum())
        st.note(f"genes={len(ids)} mapped={n_mapped} ({n_mapped/len(ids):.1%}) "
                f"kind={kind}")
        masks = class_masks(mapped)
        pair_li, pair_pi = partner_pairs(mapped)
        sym_to_pos = {}
        for pos, s in enumerate(mapped["gene_name"].astype(str).values):
            if s and s not in ("nan", "None") and s.upper() not in sym_to_pos:
                sym_to_pos[s.upper()] = pos
        st.note(f"divergent partner pairs in matrix: {pair_li.size}")

        # descriptive class annotation (LTR / TE / phyloP / GC)
        ltr_vec = mapped["tss_in_ltr"].fillna(0).values.astype(float) == 1
        ann_rows = []
        for c in ALL_GROUPS:
            m = masks[c]
            n = int(m.sum())
            row = dict(cls=c, n=n, n_ltr=int((m & ltr_vec).sum()),
                       frac_ltr=float((m & ltr_vec).sum() / max(n, 1)))
            for col, out in [("phylop_mean", "median_phylop"),
                             ("phast_mean", "median_phast"),
                             ("gc", "median_gc"),
                             ("te_frac_LTR", "median_te_frac_LTR")]:
                row[out] = (safe_median(pd.to_numeric(
                    mapped.loc[m, col], errors="coerce").values)
                    if col in mapped.columns else np.nan)
            ann_rows.append(row)
        CA = pd.DataFrame(ann_rows)
        CA.to_csv(safe_write_path(os.path.join(
            OUT, "class_annotation.tsv")), sep="\t", index=False)

        # control mean profile
        ctrl_mean = logcpm[ctrl].mean(axis=0)
        # expression bins from control profile (10 quantile bins)
        n_bins = args.bins
        expr_ok = ctrl_mean >= 1.0
        bin_id = np.full(len(ids), -1, dtype=int)
        if expr_ok.sum() >= n_bins * 10:
            q = np.quantile(ctrl_mean[expr_ok],
                            np.linspace(0, 1, n_bins + 1))
            q[0] -= 1e-9
            q[-1] += 1e-9
            bin_id[expr_ok] = np.clip(np.digitize(ctrl_mean[expr_ok], q) - 1,
                                      0, n_bins - 1)

        # group rows by perturbation
        groups = {}
        for i, v in enumerate(vals):
            groups.setdefault(v, []).append(i)
        tgt_sets = {t: [k for k in groups if target_re(t).match(k)]
                    for t in TARGETS}
        ctrl_keys = [k for k in groups if CONTROL_RE.match(k)]

        def group_profile(idx):
            w = row_cells[idx]
            if np.sum(w) > 0:
                return np.average(logcpm[idx], axis=0, weights=w)
            return logcpm[idx].mean(axis=0)

        def evaluate(label, idx, force=False):
            cells = float(np.sum(row_cells[idx]))
            if cells < args.min_cells and not force:
                return None
            pm = group_profile(idx)
            resp = pm - ctrl_mean
            det = ((pm + ctrl_mean) / 2.0) >= 1.0
            met = response_metrics(resp, det, masks, bin_id, n_bins)
            met.update(perturbation=label, n_rows=len(idx), n_cells=cells,
                       is_target=int(any(target_re(t).match(label)
                                         for t in TARGETS)),
                       is_control=int(bool(CONTROL_RE.match(label))))
            met.update(partner_stats(resp, det, pair_li, pair_pi))
            sp = sym_to_pos.get(str(label).upper())
            met["self_response"] = float(resp[sp]) if sp is not None else np.nan
            return met

        # ---- all perturbations (background) ----------------------------- #
        bg_rows = []
        for label, idx in groups.items():
            if CONTROL_RE.match(label):
                continue
            met = evaluate(label, idx)
            if met is not None:
                bg_rows.append(flatten_metrics(met))
        B = pd.DataFrame(bg_rows)
        if len(B) == 0:
            raise SystemExit("no perturbation groups passed --min-cells")
        st.note(f"background perturbations={len(B)} "
                f"median cells={B.n_cells.median():.0f}")

        # ---- targets ---------------------------------------------------- #
        tgt_met, test_rows, target_avail = [], [], []
        for t in TARGETS:
            keys = tgt_sets[t]
            if not keys:
                target_avail.append(dict(target=t, present=0, labels="",
                                         n_rows=0, n_cells=0))
                continue
            idx = [i for k in keys for i in groups[k]]
            met = evaluate(t, idx, force=True)
            tgt_met.append(met)
            target_avail.append(dict(target=t, present=1,
                                     labels=";".join(sorted(keys)),
                                     n_rows=len(idx),
                                     n_cells=float(np.sum(row_cells[idx]))))
            # same response vector for gene-set tests
            pm = group_profile(idx)
            resp = pm - ctrl_mean
            det = ((pm + ctrl_mean) / 2.0) >= 1.0
            def vals_of(name):
                return resp[masks[name] & det]
            pairs = [
                ("LTR_lnc", "lnc_other"),
                ("divergent", "other_lnc"),
                ("divergent", "antisense"),
                ("divergent", "intronic"),
                ("divergent", "sense_overlapping"),
                ("divergent", "intergenic"),
                ("LTR_lnc", "all_lnc"),
            ]
            for a_name, b_name in pairs:
                r = mwu_row(vals_of(a_name), vals_of(b_name),
                            f"{a_name}_vs_{b_name}", t)
                if r:
                    test_rows.append(r)
            # Fisher on up fraction: LTR vs other
            for a_name, b_name in [("LTR_lnc", "lnc_other"),
                                   ("divergent", "other_lnc")]:
                va, vb = vals_of(a_name), vals_of(b_name)
                if va.size >= 20 and vb.size >= 20:
                    a1 = int(np.sum(va > 1)); b1 = va.size - a1
                    c1 = int(np.sum(vb > 1)); d1 = vb.size - c1
                    orr, pf = fisher_exact([[a1, b1], [c1, d1]])
                    test_rows.append(dict(
                        target=t, contrast=f"{a_name}_vs_{b_name}_UP",
                        n_case=va.size, n_ctrl=vb.size,
                        median_case=a1 / va.size, median_ctrl=c1 / vb.size,
                        cliffs=float(np.log2(orr)) if orr > 0 else np.nan,
                        pval=float(pf)))

        T = pd.DataFrame([flatten_metrics(m) for m in tgt_met])
        TST = pd.DataFrame(test_rows)
        if len(TST):
            TST["qval"] = bh_fdr(TST.pval.values)
        AV = pd.DataFrame(target_avail)
        AV.to_csv(safe_write_path(os.path.join(
            OUT, "target_availability.tsv")), sep="\t", index=False)
        if T.empty:
            st.note("none of the requested targets found in this file")
            st.summary_md("results/perturbseq/step_summaries/39_perturbseq.md",
                          "39 Perturb-seq class-level response (Replogle 2022)",
                          "## Target availability\n\n"
                          + df_to_md(AV, max_rows=20))
            return

        # ---- target vs background -------------------------------------- #
        bg = B[~B.is_target.astype(bool)]
        vb_rows = []
        metrics = (["global_median"]
                   + [f"{c}_{m}" for c in CLASSES + ["LTR_lnc", "lnc_other",
                                                     "PCG", "all_lnc",
                                                     "other_lnc"]
                      for m in ["median", "mean_abs", "matched"]]
                   + ["ltr_minus_other", "divergent_minus_otherlnc",
                      "partner_sign_conc", "partner_pearson"])
        for t in TARGETS:
            if t not in set(T.perturbation):
                continue
            tr = T[T.perturbation == t].iloc[0]
            for met in metrics:
                if met in tr.index and met in bg.columns:
                    bv = pd.to_numeric(bg[met], errors="coerce").dropna()
                    tv = float(tr[met])
                    if not np.isfinite(tv) or len(bv) < 20:
                        continue
                    p_enrich = (1 + int(np.sum(bv >= tv))) / (1 + len(bv))
                    p_deplete = (1 + int(np.sum(bv <= tv))) / (1 + len(bv))
                    vb_rows.append(dict(
                        target=t, metric=met, target_value=tv,
                        bg_median=float(bv.median()),
                        bg_sd=float(bv.std()),
                        percentile=float((bv < tv).mean() * 100),
                        z=float((tv - bv.mean()) / bv.std())
                        if bv.std() > 0 else np.nan,
                        p_enrich=float(p_enrich), p_deplete=float(p_deplete)))
        VB = pd.DataFrame(vb_rows)
        if len(VB):
            VB["q_enrich"] = bh_fdr(VB.p_enrich.values)
            VB["q_deplete"] = bh_fdr(VB.p_deplete.values)

        # ---- all-perturbation ranking ---------------------------------- #
        B["rank_divergent_minus_otherlnc"] = B[
            "divergent_minus_otherlnc"].rank(ascending=False, method="min")
        B["rank_ltr_minus_other"] = B["ltr_minus_other"].rank(
            ascending=False, method="min")
        B["is_target"] = B["is_target"].astype(int)
        B.to_csv(safe_write_path(os.path.join(
            OUT, "perturbation_scores.tsv")), sep="\t", index=False)

        # ---- write outputs --------------------------------------------- #
        AV.to_csv(safe_write_path(os.path.join(
            OUT, "target_availability.tsv")), sep="\t", index=False)

        # long class response table for targets
        long_rows = []
        for met in tgt_met:
            t = met["perturbation"]
            for c in ALL_GROUPS:
                d = met.get(c, {})
                long_rows.append(dict(
                    target=t, cls=c,
                    **{k: d.get(k, np.nan) for k in
                       ["n", "median", "mean", "mean_abs", "frac_up",
                        "frac_down", "matched"]}))
            long_rows.append(dict(
                target=t, cls="GLOBAL", n=np.nan, median=met["global_median"],
                mean=np.nan, mean_abs=met["global_mean_abs"], frac_up=np.nan,
                frac_down=np.nan, matched=np.nan))
        L = pd.DataFrame(long_rows)
        L.to_csv(safe_write_path(os.path.join(
            OUT, "class_response_by_target.tsv")), sep="\t", index=False)
        if len(TST):
            TST.to_csv(safe_write_path(os.path.join(
                OUT, "class_tests_within_target.tsv")), sep="\t", index=False)
        if len(VB):
            VB.to_csv(safe_write_path(os.path.join(
                OUT, "target_vs_background.tsv")), sep="\t", index=False)
        qc_cols = [c for c in ["perturbation", "n_cells", "n_rows",
                               "self_response", "global_median", "partner_n",
                               "partner_sign_conc", "partner_pearson"]
                   if c in T.columns]
        QC = T[qc_cols].copy()
        QC.to_csv(safe_write_path(os.path.join(OUT, "target_qc.tsv")),
                  sep="\t", index=False)

        # ---- summary ---------------------------------------------------- #
        key = VB[VB.metric.isin(["divergent_median", "LTR_lnc_median",
                                 "ltr_minus_other",
                                 "divergent_minus_otherlnc",
                                 "partner_sign_conc", "partner_pearson",
                                 "global_median"])] if len(VB) else VB
        keys = key.sort_values(["target", "metric"]) if len(key) else key
        TSTs = TST.sort_values(["target", "pval"]) if len(TST) else TST
        top = B.sort_values("divergent_minus_otherlnc", ascending=False).head(12)
        body = (
            f"File: `{os.path.basename(f)}`\n\n"
            f"Perturbations: {len(B)} "
            f"(non-target background {len(bg)}); controls: "
            f"{int(pd.Series(ctrl).sum())} pseudobulk rows; "
            f"genes mapped {n_mapped}/{len(ids)}.\n\n"
            "## Target availability\n\n" + df_to_md(AV, max_rows=20)
            + "\n\n## Class annotation (genes present in matrix)\n\n"
            + df_to_md(CA, ".4g", max_rows=20)
            + "\n\n## Target QC (self-knockdown, cells, partner concordance)\n\n"
            + df_to_md(QC, ".3f", max_rows=20)
            + "\n\n## Class response per target\n\n"
            + df_to_md(L, ".4g", max_rows=60)
            + "\n\n## Within-target gene-set tests\n\n"
            + df_to_md(TSTs, ".4g", max_rows=40)
            + "\n\n## Target vs all-perturbation background\n\n"
            + df_to_md(keys, ".4g", max_rows=60)
            + "\n\n## Top perturbations by divergent-vs-other-lnc score\n\n"
            + df_to_md(top[["perturbation", "n_cells", "is_target",
                            "divergent_median", "other_lnc_median",
                            "divergent_minus_otherlnc",
                            "rank_divergent_minus_otherlnc"]], ".3f", 12))
        st.summary_md("results/perturbseq/step_summaries/39_perturbseq.md",
                      "39 Perturb-seq class-level response (Replogle 2022)",
                      body)
        st.note("targets: " + "; ".join(
            f"{r.perturbation} cells={r.n_cells:.0f} "
            f"divΔ={r.divergent_minus_otherlnc:.3f}"
            for _, r in T.iterrows()))
        st.note("outputs -> " + OUT)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    pi = sub.add_parser("inspect")
    pi.set_defaults(func=cmd_inspect)
    pa = sub.add_parser("analyze")
    pa.add_argument("--dataset", default="K562_gwps",
                    help="substring of bulk h5ad (default K562_gwps)")
    pa.add_argument("--file", default=None, help="explicit h5ad path")
    pa.add_argument("--pert-col", default=None)
    pa.add_argument("--min-cells", type=int, default=20)
    pa.add_argument("--bins", type=int, default=10)
    pa.set_defaults(func=cmd_analyze)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
