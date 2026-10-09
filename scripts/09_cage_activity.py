#!/usr/bin/env python
"""09_cage_activity.py — promoter activity & tissue specificity from FANTOM CAGE.

Primary: per-sample TPM matrix (fair+new) keyed by peak id (844MB osc matrix).
Fallback: peaks-only binary activity from peak score (TPM matrix unavailable).
Also: divergent-pair activity correlation (Spearman, vs same-strand control).

Outputs: results/cage_activity.tsv, results/divergent_pairs_activity.tsv,
         results/step_summaries/09_cage_activity.md
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR, df_to_md

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data")
TPM = os.path.join(DATA, "fantom", "hg38_cage_tpm_ann.osc.txt.gz")
TPM_MIN = 1.0


def read_osc_matrix(path):
    """Parse FANTOM .osc.txt.gz: ## metadata lines, then TSV matrix whose first
    column is the peak id. Returns (peak_ids, sample_names, matrix ndarray)."""
    meta_lines, header, rows = [], None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                meta_lines.append(line)
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            rows.append(f)
    hdr0 = header[0]
    names = header[1:]
    ids = [r[0] for r in rows]
    # keep only numeric sample columns (ann variant appends annotation cols)
    ncol = []
    keep_idx = []
    for j, nm in enumerate(names, start=1):
        try:
            float(rows[0][j]) if len(rows[0]) > j else None
            keep_idx.append(j)
            ncol.append(nm)
        except ValueError:
            continue
    M = np.array([[float(r[j]) for j in keep_idx] for r in rows], dtype=np.float32)
    return ids, ncol, M


def tau_specificity(mat):
    """tau tissue-specificity per row (max over columns as reference)."""
    x = mat / np.maximum(mat.max(axis=1, keepdims=True), 1e-9)
    n = mat.shape[1]
    return ((1 - x).sum(axis=1) / (n - 1))


def main():
    with step("09_cage_activity") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                           index_col=0)

        if os.path.exists(TPM) and os.path.getsize(TPM) > 1e6:
            st.note("loading TPM matrix (large, may take a while)")
            ids, names, M = read_osc_matrix(TPM)
            st.note(f"matrix: {len(ids)} peaks x {len(names)} samples")
            pid2row = {p: i for i, p in enumerate(ids)}
            active = (M > TPM_MIN).any(axis=1)
            maxtpm = M.max(axis=1)
            meantpm = M.mean(axis=1)
            breadth = (M > TPM_MIN).sum(axis=1)
            tau = tau_specificity(M)
            idx = meta["cage_peak_id"].map(pid2row)
            ok = idx.notna()
            act = pd.DataFrame(index=meta.index)
            act["active"] = 0
            act.loc[ok, "active"] = active[idx[ok].astype(int).values].astype(int)
            act["max_tpm"] = np.nan
            act.loc[ok, "max_tpm"] = maxtpm[idx[ok].astype(int).values]
            act["mean_tpm"] = np.nan
            act.loc[ok, "mean_tpm"] = meantpm[idx[ok].astype(int).values]
            act["breadth"] = np.nan
            act.loc[ok, "breadth"] = breadth[idx[ok].astype(int).values]
            act["tau"] = np.nan
            act.loc[ok, "tau"] = tau[idx[ok].astype(int).values]
            st.note(f"promoters linked to TPM rows: {int(ok.sum())} ({ok.mean():.0%})")

            # ---- divergent pair activity correlation
            dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
            rhos, n_used = [], 0
            prof_cache = {}
            for _, r in dp.iterrows():
                a, b = r["lnc_gid"], r["pcg_gid"]
                ia, ib = meta.at[a, "cage_peak_id"], meta.at[b, "cage_peak_id"]
                if ia not in pid2row or ib not in pid2row:
                    continue
                for gid, pid in ((a, ia), (b, ib)):
                    if pid not in prof_cache:
                        prof_cache[pid] = M[pid2row[pid]].astype(np.float64)
                va, vb = prof_cache[ia], prof_cache[ib]
                if va.std() == 0 or vb.std() == 0:
                    continue
                rho = spearmanr(va, vb).statistic
                if not np.isnan(rho):
                    rhos.append(rho)
                    n_used += 1
            # control: random same-strand pairs matched by distance
            rng = np.random.default_rng(7)
            pcg_ids = meta[(meta.group == "pcg") & meta.cage_peak_id.notna()]
            lnc_ids = meta[(meta.group == "lnc") & meta.cage_peak_id.notna()]
            ctrl = []
            cand = list(zip(lnc_ids.index, lnc_ids["chrom"], lnc_ids["tss"]))
            pcg_list = list(zip(pcg_ids.index, pcg_ids["chrom"], pcg_ids["tss"]))
            from collections import defaultdict
            bychrom = defaultdict(list)
            for g, c, t in pcg_list:
                bychrom[c].append((g, t))
            for g, c, t in cand:
                for g2, t2 in bychrom.get(c, []):
                    d = abs(t - t2)
                    if 200 <= d <= 2000:
                        ia, ib = meta.at[g, "cage_peak_id"], meta.at[g2, "cage_peak_id"]
                        if ia in prof_cache and ib in prof_cache:
                            va, vb = prof_cache[ia], prof_cache[ib]
                            if va.std() > 0 and vb.std() > 0:
                                r = spearmanr(va, vb).statistic
                                if not np.isnan(r):
                                    ctrl.append(r)
                        break
            pa = pd.DataFrame({
                "rho_divergent": pd.Series(rhos),
                "rho_control_same_strand": pd.Series(ctrl)})
            pa.to_csv(safe_write_path(
                os.path.join(RES, "divergent_pairs_activity.tsv")), sep="\t",
                index=False)
            st.note(f"divergent pairs used={n_used} "
                    f"median rho={np.median(rhos) if rhos else float('nan'):.2f}; "
                    f"control median rho={np.median(ctrl) if ctrl else float('nan'):.2f}")
            mode = "TPM"
        else:
            st.note("TPM matrix unavailable — peaks-only binary activity mode")
            act = pd.DataFrame(index=meta.index)
            act["active"] = (meta["cage_refined"] == 1).astype(int)
            act["max_tpm"] = meta["cage_peak_score"]
            act["mean_tpm"] = np.nan
            act["breadth"] = np.nan
            act["tau"] = np.nan
            mode = "peaks-only"

        out = meta.join(act)
        out.to_csv(safe_write_path(os.path.join(RES, "cage_activity.tsv")), sep="\t")

        groups = ["PCG"] + ["divergent", "antisense", "intronic",
                            "sense_overlapping", "intergenic",
                            "pseudogene", "transcribed_pseudogene"]
        rows = []
        for g in groups:
            sub = out[(out["group"] == "pcg")] if g == "PCG" else out[out["cls"] == g]
            if not len(sub):
                continue
            rows.append(dict(cls=g, n=len(sub),
                             active_frac=sub["active"].mean(),
                             log10_maxtpm=np.log10(sub["max_tpm"].fillna(0) + 1).mean(),
                             tau=sub["tau"].mean() if mode == "TPM" else np.nan,
                             breadth=sub["breadth"].mean() if mode == "TPM" else np.nan,
                             enhancer_frac=sub["enhancer_overlap"].mean()))
        sm = pd.DataFrame(rows)
        sm.to_csv(safe_write_path(os.path.join(RES, "cage_activity_by_class.tsv")),
                  sep="\t", index=False)
        st.note("activity by class: " + "; ".join(
            f"{r.cls}: active={r.active_frac:.2f}" for r in sm.itertuples()))
        st.summary_md("results/step_summaries/09_cage_activity.md",
                      f"09 promoter activity ({mode})", df_to_md(sm, ".2f"))


if __name__ == "__main__":
    main()
