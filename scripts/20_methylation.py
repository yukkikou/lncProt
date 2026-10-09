#!/usr/bin/env python
"""20_methylation.py — CGGA IDH_A paired EPIC methylation x RNA (n=28).

Global methylation (median beta per sample) vs relative activation of
LTR-derived lncRNAs (LTR minus other lncRNAs), with IDH status and grade
from the clinical xlsx. Spearman correlations + partial (controlling IDH).
Outputs -> results/ccle_cgga/tables/methylation_vs_ltr.tsv + summary.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

TD = os.path.join(WORKDIR, "data", "ccle_cgga", "cgga", "CGGA_CCell_4083")
OUT = os.path.join(WORKDIR, "results", "ccle_cgga", "tables")
RES = os.path.join(WORKDIR, "results")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def partial_spearman(x, y, z):
    rx, ry, rz = pd.Series(x).rank(), pd.Series(y).rank(), pd.Series(z).rank()
    A = np.c_[rz, np.ones(len(rz))]
    bx, *_ = np.linalg.lstsq(A, rx.values, rcond=None)
    by, *_ = np.linalg.lstsq(A, ry.values, rcond=None)
    ex = rx.values - A @ bx
    ey = ry.values - A @ by
    return spearmanr(ex, ey)


def main():
    with step("20_methylation") as st:
        os.makedirs(OUT, exist_ok=True)
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", low_memory=False)
        sym2cls = dict(zip(meta.gene_name.astype(str), meta.cls))
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["tss_in_ltr"]].join(
            meta.set_index("gene_id")["gene_name"])
        sym2ltr = {str(s): int(v) for s, v in zip(te.gene_name, te.tss_in_ltr)
                   if pd.notna(s)}

        beta = pd.read_csv(os.path.join(
            TD, "CGGA_IDH_A_Methylation_EPIC_Array_20250915.txt"),
            sep="\t", index_col=0, nrows=None)
        global_beta = beta.median(axis=0)
        st.note(f"methylation: {beta.shape[0]:,} probes x {beta.shape[1]} samples")

        rna = pd.read_csv(os.path.join(TD, "CGGA_IDH_A_RNAseq_RSEM_20250915.txt"),
                          sep="\t", index_col=0)
        common = global_beta.index.intersection(rna.columns)
        st.note(f"paired samples: {len(common)}")
        rna = rna[common]
        rl = np.log2(rna.apply(pd.to_numeric, errors="coerce") + 1)
        cls = np.array([sym2cls.get(g) for g in rl.index])
        ltr = np.array([sym2ltr.get(g, 0) for g in rl.index])
        lnc = np.isin(cls, CLASSES)
        rel = (rl[lnc & (ltr == 1)].mean(axis=0)
               - rl[lnc & (ltr == 0)].mean(axis=0))
        ltr_abs = rl[lnc & (ltr == 1)].mean(axis=0)
        glob = global_beta.loc[common]

        # clinical (IDH/grade if present)
        clin = pd.read_excel(os.path.join(TD, "CGGA_IDH_A_Clinical_Information.xlsx"))
        clin.iloc[:, 0] = clin.iloc[:, 0].astype(str)
        clin = clin.set_index(clin.columns[0])
        idh_col = next((c for c in clin.columns if "idh" in str(c).lower()), None)
        idh = clin[idh_col].reindex(common).astype(str) if idh_col else None

        ok = rel.notna() & glob.notna()
        rho_g, p_g = spearmanr(glob[ok], rel[ok])
        rho_a, p_a = spearmanr(glob[ok], ltr_abs[ok])
        rows = [dict(contrast="global beta vs LTR-relative activation",
                     rho=rho_g, p=p_g, n=int(ok.sum())),
                dict(contrast="global beta vs LTR-lnc absolute mean",
                     rho=rho_a, p=p_a, n=int(ok.sum()))]
        # partial controlling IDH (encode WT=0, Mut=1)
        if idh is not None and idh.notna().sum() > 10:
            z = idh.str.contains("Wild|WT|wilde", regex=True).astype(float)
            z = z.fillna(1.0 - z).reindex(common).fillna(0.5)
            ok2 = ok & z.notna()
            rp, pp = partial_spearman(glob[ok2], rel[ok2], z[ok2])
            rows.append(dict(contrast="global beta vs LTR-relative (IDH-controlled)",
                             rho=rp, p=pp, n=int(ok2.sum())))
        t = pd.DataFrame(rows)
        t.to_csv(safe_write_path(os.path.join(OUT, "methylation_vs_ltr.tsv")),
                 sep="\t", index=False)
        st.note("; ".join(f"{r.contrast}: rho={r.rho:.2f} p={r.p:.1e}"
                          for r in t.itertuples()))
        st.summary_md("results/ccle_cgga/step_summaries/20_methylation.md",
                      "20 Paired methylation x LTR activation (n=%d)" % len(common),
                      df_to_md(t, ".3g"))


if __name__ == "__main__":
    main()
