#!/usr/bin/env python
"""83_ifn_ltr_tumour.py — is tumour LTR de-repression an interferon response?

DUX4 is absent from TCGA tumours, so the tumour-context LTR signal must have a
DUX4-independent explanation. Global demethylation and ERV induction are
frequently interferon-linked. This module asks whether the per-sample
LTR-relative lncRNA score tracks an interferon/STAT1 signature across TCGA
tumours, within and across cancer types.

Outputs
  results/revision_checks/G_ifn_ltr_per_cancer.tsv
  results/revision_checks/G_ifn_ltr_summary.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

TOIL = os.path.join(WORKDIR, "data", "tissue_disease", "toil_class_expr.npz")
TGEN = os.path.join(WORKDIR, "data", "tissue_disease", "toil_genes.tsv")
TSAM = os.path.join(WORKDIR, "data", "tissue_disease", "toil_samples.tsv")
RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "revision_checks")
IFN = ["STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP4",
       "HLA-DRA", "IDO1", "ISG15"]


def main():
    with step("83_ifn_ltr_tumour") as st:
        os.makedirs(OUT, exist_ok=True)
        g = pd.read_csv(TGEN, sep="\t")
        s = pd.read_csv(TSAM, sep="\t")
        z = np.load(TOIL, mmap_mode="r")
        M = np.asarray(z["M"])
        keep = ((s.cohort == "TCGA") & (s.condition == "tumor")).values
        cols = np.where(keep)[0]
        cancer = np.array([x if isinstance(x, str) else ""
                           for x in s.cancer.values[cols]])
        st.note("TCGA tumour samples={}".format(len(cols)))

        sym = g.symbol.astype(str).values
        ltr_rows = np.where((g.cls.notna())
                            & g.cls.isin(["divergent", "antisense", "intronic",
                                          "sense_overlapping", "intergenic"])
                            & (g.tss_in_ltr == 1))[0]
        other_lnc = np.where((g.cls.notna())
                             & g.cls.isin(["divergent", "antisense",
                                           "intronic", "sense_overlapping",
                                           "intergenic"])
                             & (g.tss_in_ltr != 1))[0]
        st.note("LTR lncRNAs={} other lncRNAs={}".format(
            len(ltr_rows), len(other_lnc)))
        X = M[:, cols]
        ltr = np.nanmean(X[ltr_rows, :], axis=0)
        oth = np.nanmean(X[other_lnc, :], axis=0)
        score = ltr - oth

        sig_rows = [np.where(sym == x)[0][0] for x in IFN if (sym == x).any()]
        st.note("IFN signature genes found: {}".format(len(sig_rows)))
        Z = X[sig_rows, :].astype(np.float64)
        Z = (Z - np.nanmean(Z, axis=1, keepdims=True)) / np.nanstd(
            Z, axis=1, keepdims=True)
        ifn = np.nanmean(Z, axis=0)

        rows = []
        for c in sorted(set(cancer)):
            idx = np.where(cancer == c)[0]
            if len(idx) < 30:
                continue
            r, p = spearmanr(score[idx], ifn[idx])
            rows.append(dict(cancer=c, n=len(idx), rho=float(r), pval=float(p)))
        T = pd.DataFrame(rows)
        pv = T.pval.values; o = np.argsort(pv)
        q = pv[o] * len(pv) / (np.arange(len(pv)) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        T["fdr"] = np.clip(q, 0, 1)[np.argsort(o)]
        T.to_csv(safe_write_path(os.path.join(
            OUT, "G_ifn_ltr_per_cancer.tsv")), sep="\t", index=False)

        # pooled within-cancer rank residual correlation
        from scipy.stats import rankdata
        df = pd.DataFrame(dict(score=score, ifn=ifn, cancer=cancer))
        df = df.dropna()
        rz = df.score.rank().values
        ri = df.ifn.rank().values
        dm = pd.get_dummies(df.cancer, drop_first=True).astype(float).values
        Xd = np.column_stack([dm, np.ones(len(df))])
        ez = rz - Xd @ np.linalg.lstsq(Xd, rz, rcond=None)[0]
        ei = ri - Xd @ np.linalg.lstsq(Xd, ri, rcond=None)[0]
        rho, p = spearmanr(ez, ei)
        S = pd.DataFrame([dict(n=len(df), n_cancers=len(T),
                               median_rho=float(T.rho.median()),
                               frac_positive=float((T.rho > 0).mean()),
                               pooled_partial_rho=float(rho),
                               pooled_pval=float(p))])
        S.to_csv(safe_write_path(os.path.join(
            OUT, "G_ifn_ltr_summary.tsv")), sep="\t", index=False)
        st.note("per-cancer median rho={:.3f}, positive in {:.0%}; pooled "
                "partial rho={:.3f} P={:.2g}".format(
                    T.rho.median(), (T.rho > 0).mean(), rho, p))
        st.summary_md(os.path.join(
            RES, "step_summaries", "83_ifn_ltr_tumour.md"), "83 IFN signature versus LTR score",
            df_to_md(T, ".3f") + "\n\n" + df_to_md(S, ".4f"))


if __name__ == "__main__":
    main()
