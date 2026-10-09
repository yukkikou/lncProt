#!/usr/bin/env python
"""13c_gtex_tissue.py — tissue dimension (GTEx): how does each birth class
behave across 30+ normal tissues?

Per tissue (n>=20): class expression fraction, median log2TPM, and
"route usage" = share of expressed lncRNAs by birth class + LTR share.
Class-level tissue specificity (tau across tissues, GTEx-based).

Outputs -> results/tissue_disease/tables/:
  gtex_class_by_tissue.tsv, gtex_route_usage.tsv,
  gtex_class_specificity.tsv, gtex_top_per_tissue.tsv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

TD = os.path.join(WORKDIR, "data", "tissue_disease")
OUT = os.path.join(WORKDIR, "results", "tissue_disease", "tables")
MIN_N = 20
EXPR = 1.0          # log2(TPM+1) >= 1  <=>  TPM >= 1
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def main():
    with step("13c_gtex_tissue") as st:
        os.makedirs(OUT, exist_ok=True)
        d = np.nan_to_num(np.load(os.path.join(TD, "toil_class_expr.npz"))["M"])
        genes = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                            index_col=0)
        samp = pd.read_csv(os.path.join(TD, "toil_samples.tsv"), sep="\t")
        gi = genes  # rows aligned with M by construction
        gtex = samp[(samp.cohort == "GTEx") & samp.tissue.notna()]
        tissues = [t for t, sub in gtex.groupby("tissue") if len(sub) >= MIN_N]
        st.note(f"tissues used (n>={MIN_N}): {len(tissues)}")
        if not tissues:
            raise RuntimeError("no GTEx tissues with n>=20 — check cohort/tissue annotation")

        cls_idx = {c: (gi.cls == c).values for c in CLASSES}
        pcg_idx = (gi.cls == "PCG").values
        ltr_idx = (gi.tss_in_ltr == 1).values
        rows = []
        for t in tissues:
            cols = gtex[gtex.tissue == t].index.values
            sub = d[:, cols]
            row = dict(tissue=t, n_samples=len(cols))
            for c in CLASSES + ["PCG"]:
                idx = pcg_idx if c == "PCG" else cls_idx[c]
                vals = sub[idx]
                row[f"{c}_frac_expr"] = float((vals >= EXPR).any(axis=1).mean())
                med = np.median(vals[vals >= EXPR]) if (vals >= EXPR).any() else 0
                row[f"{c}_med_expr"] = float(med)
            # route usage among EXPRESSED lncRNAs
            lnc_all = np.zeros(len(gi), bool)
            for c in CLASSES:
                lnc_all |= cls_idx[c]
            expressed = lnc_all & (sub >= EXPR).any(axis=1)
            n_exp = int(expressed.sum())
            for c in CLASSES:
                row[f"share_{c}"] = float((expressed & cls_idx[c]).sum()
                                          / max(n_exp, 1))
            row["share_LTR"] = float((expressed & ltr_idx).sum() / max(n_exp, 1))
            row["n_lnc_expressed"] = n_exp
            rows.append(row)
        t1 = pd.DataFrame(rows)
        t1.to_csv(safe_write_path(os.path.join(OUT, "gtex_class_by_tissue.tsv")),
                  sep="\t", index=False)
        st.note(f"tissue rows={len(t1)}")

        # ---- class-level GTEx tissue specificity (tau over tissue means)
        gcols = {t: gtex[gtex.tissue == t].index.values for t in tissues}
        mean_by_t = np.vstack([d[:, c].mean(axis=1) for t, c in gcols.items()])
        mx = mean_by_t.max(axis=0)
        expressed = mx >= 1.0          # max tissue mean TPM >= 1
        tau = np.full(len(gi), np.nan)
        tau[expressed] = ((1 - mean_by_t[:, expressed] / mx[expressed]).sum(0)
                          / (len(tissues) - 1))
        spec = gi[["symbol", "cls"]].copy()
        spec["tau_gtex"] = tau
        spec_rows = []
        for c in CLASSES + ["PCG"]:
            idx = pcg_idx if c == "PCG" else cls_idx[c]
            v = tau[idx]
            spec_rows.append(dict(cls=c, n=int(idx.sum()), tau_mean=float(v.mean()),
                                  tau_median=float(np.median(v))))
        t3 = pd.DataFrame(spec_rows)
        t3.to_csv(safe_write_path(os.path.join(OUT, "gtex_class_specificity.tsv")),
                  sep="\t", index=False)
        spec.to_csv(safe_write_path(os.path.join(OUT, "gtex_gene_tau.tsv")),
                    sep="\t")
        st.note("tau: " + "; ".join(f"{r.cls}={r.tau_median:.3f}"
                                    for r in t3.itertuples()))

        # ---- top examples per tissue (highest tissue-specific lncRNA)
        tops = []
        for t in tissues:
            cols = gcols[t]
            here = d[:, cols].mean(axis=1)
            elsewhere = np.mean([d[:, c].mean(axis=1) for t2, c in gcols.items()
                                 if t2 != t], axis=0)
            fc = here - elsewhere
            lnc_mask = lnc_all
            order = np.argsort(-fc * lnc_mask)[:5]
            for i in order:
                tops.append(dict(tissue=t, symbol=gi.symbol.iloc[i],
                                 cls=gi.cls.iloc[i], tss_in_ltr=int(ltr_idx[i]),
                                 delta_log2tpm=float(fc[i])))
        t4 = pd.DataFrame(tops)
        t4.to_csv(safe_write_path(os.path.join(OUT, "gtex_top_per_tissue.tsv")),
                  sep="\t", index=False)

        body = ("## Class expression by tissue (n>=20)\n\n"
                + df_to_md(t1.head(40), ".3f", max_rows=40)
                + "\n\n## GTEx tissue specificity (tau) by class\n\n"
                + df_to_md(t3, ".3f"))
        st.summary_md("results/tissue_disease/step_summaries/13c_gtex_tissue.md",
                      "13c GTEx tissue dimension", body)


if __name__ == "__main__":
    main()
