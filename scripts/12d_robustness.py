#!/usr/bin/env python
"""12d_robustness.py — reviewer-proofing controls:

 1. Divergent-window sensitivity: class counts at 500/1000/2000 bp
 2. Activity-matched control: do class differences survive within
    CAGE-active promoters only?
 3. GC partial correlation: class -> phastCons controlling GC content
    (rank-residual approach)
Outputs: results/robustness_sensitivity.tsv, results/robustness_active.tsv,
         results/robustness_partial_corr.tsv + summary
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")


def partial_spearman(x, y, z):
    """Spearman partial correlation of x,y controlling z (rank residuals)."""
    rx, ry, rz = pd.Series(x).rank(), pd.Series(y).rank(), pd.Series(z).rank()
    def resid(a):
        A = np.c_[rz, np.ones(len(rz))]
        b, *_ = np.linalg.lstsq(A, a.values, rcond=None)
        return a.values - A @ b
    ex, ey = resid(rx), resid(ry)
    r, p = spearmanr(ex, ey)
    return r, p


def main():
    with step("12d_robustness") as st:
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        seq = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["gc"]]
        cons = pd.read_csv(os.path.join(RES, "conservation.tsv"), sep="\t",
                           index_col=0, low_memory=False)[["phast_mean"]]
        act = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["active"]]
        ch = pd.read_csv(os.path.join(RES, "chromhmm_states.tsv"), sep="\t",
                         index_col=0, low_memory=False)
        m = m.join(seq).join(cons).join(act)
        m["prom_like"] = ch.eq("TssA").mean(axis=1).reindex(m.index).values

        # ---- 1. divergent window sensitivity
        lnc = m[m.group == "lnc"]
        counts = {"500bp": int(lnc["div_500"].sum()),
                  "1000bp (primary)": int((lnc["cls"] == "divergent").sum()),
                  "2000bp": int(lnc["div_2000"].sum())}
        sens = pd.DataFrame([counts])
        sens.to_csv(safe_write_path(
            os.path.join(RES, "robustness_sensitivity.tsv")), sep="\t",
            index=False)
        st.note(f"divergent sensitivity: {counts}")

        # ---- 2. active-only control
        rows = []
        active = m[(m.active == 1)]
        for feat, label in [("phast_mean", "phastCons"),
                            ("prom_like", "prom-like chromatin frac"),
                            ("gc", "GC content")]:
            pcg_v = active.loc[active.group == "pcg", feat].dropna()
            for cls in ["divergent", "antisense", "intronic", "intergenic"]:
                v = active.loc[active.cls == cls, feat].dropna()
                if len(v) < 20 or len(pcg_v) < 20:
                    continue
                rows.append(dict(feature=label, cls=cls, n_active_class=len(v),
                                 mean_active_class=float(v.mean()),
                                 mean_active_PCG=float(pcg_v.mean()),
                                 diff=float(v.mean() - pcg_v.mean())))
        at = pd.DataFrame(rows)
        at.to_csv(safe_write_path(os.path.join(RES, "robustness_active.tsv")),
                  sep="\t", index=False)
        st.note(f"active-only rows={len(at)} (n active PCG="
                f"{int((active.group=='pcg').sum())})")

        # ---- 3. partial correlation class-conservation controlling GC
        prows = []
        sub = m[m.group.isin(["pcg", "lnc"]) & m.phast_mean.notna() & m.gc.notna()]
        # encode each lnc class vs PCG as 0/1, correlate with conservation
        for cls in ["divergent", "antisense", "intronic", "sense_overlapping",
                    "intergenic"]:
            s = sub[(sub.group == "pcg") | (sub.cls == cls)].copy()
            x = (s["cls"] == cls).astype(int).values
            r0, p0 = spearmanr(x, s["phast_mean"].values)
            r, p = partial_spearman(pd.Series(x), s["phast_mean"], s["gc"])
            prows.append(dict(cls=cls, n=len(s), rho_raw=r0, p_raw=p0,
                              rho_GC_controlled=r, p_controlled=p))
        pt = pd.DataFrame(prows)
        pt.to_csv(safe_write_path(
            os.path.join(RES, "robustness_partial_corr.tsv")), sep="\t",
            index=False)
        st.note("GC-controlled conservation: " + "; ".join(
            f"{r.cls}: {r.rho_raw:.3f}->{r.rho_GC_controlled:.3f}"
            for r in pt.itertuples()))

        body = ("## Divergent window sensitivity\n\n" + df_to_md(sens, ".0f")
                + "\n\n## Active-promoter-only comparison\n\n"
                + df_to_md(at, ".4f", max_rows=30)
                + "\n\n## Conservation vs class, GC-controlled (partial Spearman)\n\n"
                + df_to_md(pt, ".3f"))
        st.summary_md("results/step_summaries/12d_robustness.md",
                      "12d robustness controls", body)


if __name__ == "__main__":
    main()
