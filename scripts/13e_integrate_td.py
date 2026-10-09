#!/usr/bin/env python
"""13e_integrate_td.py — integrate tissue (GTEx) x disease (TCGA) x promoter
classes; cross-check FANTOM promoter activity vs GTEx body expression.

Outputs -> results/tissue_disease/tables/summary_tissue_disease.tsv,
           fantom_vs_gtex.tsv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

TD = os.path.join(WORKDIR, "data", "tissue_disease")
OUT = os.path.join(WORKDIR, "results", "tissue_disease", "tables")
RES = os.path.join(WORKDIR, "results")


def main():
    with step("13e_integrate_td") as st:
        os.makedirs(OUT, exist_ok=True)
        gt = pd.read_csv(os.path.join(OUT, "gtex_class_by_tissue.tsv"), sep="\t")
        de = pd.read_csv(os.path.join(OUT, "tcga_class_de_per_cancer.tsv"),
                         sep="\t")
        rows = []
        for c in ["divergent", "antisense", "intronic", "sense_overlapping",
                  "intergenic", "PCG"]:
            sub = de[de.cls == c]
            gt_sub = gt[f"share_{c}"].mean() if f"share_{c}" in gt else np.nan
            rows.append(dict(class_or_group=c,
                             gtex_mean_share_of_expressed_lnc=gt_sub,
                             gtex_median_frac_expr=np.nan,
                             n_cancers_de=sub.cancer.nunique(),
                             median_frac_up=sub.frac_up.median(),
                             median_delta=sub.median_delta.median()))
        # fill frac expr from gt columns
        for r in rows:
            col = f"{r['class_or_group']}_frac_expr"
            if col in gt.columns:
                r["gtex_median_frac_expr"] = gt[col].median()
        t = pd.DataFrame(rows)
        t.to_csv(safe_write_path(os.path.join(OUT, "summary_tissue_disease.tsv")),
                 sep="\t", index=False)
        st.note("summary: " + "; ".join(
            f"{r.class_or_group}: share={r.gtex_mean_share_of_expressed_lnc:.2f}"
            for r in t.itertuples()))

        # FANTOM tau vs GTEx tau consistency (gene level)
        tau_g = pd.read_csv(os.path.join(OUT, "gtex_gene_tau.tsv"), sep="\t",
                            index_col=0)
        from scipy.stats import spearmanr
        fant = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                           index_col=0, low_memory=False)[["tau"]]
        fant["base_id"] = [g.split(".")[0] for g in fant.index]
        f2 = fant.reset_index(drop=True).set_index("base_id")["tau"]
        common = tau_g.index.intersection(f2.index)
        if len(common) > 100:
            rho, p = spearmanr(tau_g.loc[common, "tau_gtex"], f2.loc[common])
            fv = pd.DataFrame({"rho": [rho], "p": [p], "n": [len(common)]})
            fv.to_csv(safe_write_path(os.path.join(OUT, "fantom_vs_gtex.tsv")),
                      sep="\t", index=False)
            st.note(f"FANTOM tau vs GTEx tau: rho={rho:.2f} (n={len(common)})")

        st.summary_md("results/tissue_disease/step_summaries/13e_integration.md",
                      "13e tissue x disease integration",
                      df_to_md(t, ".3f"))


if __name__ == "__main__":
    main()
