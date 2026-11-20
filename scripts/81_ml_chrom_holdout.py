#!/usr/bin/env python
"""81_ml_chrom_holdout.py — birth-class classifier with chromosome-held-out CV.

The main classifier used random stratified folds. Loci are spatially clustered,
so a random split can place a divergent lncRNA and its partner, or two nearby
loci, on opposite sides of the fold boundary. This module repeats the same
classifier and features with GroupKFold by chromosome and reports the drop.

Outputs
  results/revision_checks/E_ml_chrom_holdout.tsv
  results/revision_checks/E_ml_chrom_per_class.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import f1_score, confusion_matrix, cohen_kappa_score, \
    balanced_accuracy_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "revision_checks")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
SEED = 42


def build():
    m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                    index_col=0, low_memory=False)
    parts = [m[["cls", "group", "chrom"]]]
    for f, cols in [("seq_features.tsv", ["gc", "cpg_oe", "cpg_density",
                                          "TATA", "Inr", "CCAAT", "GCbox",
                                          "cpg_island", "prom_class"]),
                    ("te_promoter_stats.tsv", ["te_cov", "te_frac_LTR",
                                               "te_frac_LINE", "te_frac_SINE",
                                               "te_frac_DNA", "tss_in_te",
                                               "tss_in_ltr"]),
                    ("conservation.tsv", ["phast_mean"])]:
        d = pd.read_csv(os.path.join(RES, f), sep="\t", index_col=0,
                        low_memory=False)
        parts.append(d[[c for c in cols if c in d.columns]])
    ch = pd.read_csv(os.path.join(RES, "chromhmm_states.tsv"), sep="\t",
                     index_col=0, low_memory=False)
    parts.append(pd.concat([
        ch.eq("TssA").mean(axis=1).rename("prom_like"),
        ch.isin(["EnhG", "Enh"]).mean(axis=1).rename("enh_like")], axis=1))
    X = pd.concat(parts, axis=1)
    X["prom_class"] = X["prom_class"].map({"HCP": 2, "ICP": 1, "LCP": 0})
    return X


def main():
    with step("81_ml_chrom_holdout") as st:
        os.makedirs(OUT, exist_ok=True)
        X = build()
        sub = X[X.cls.isin(CLASSES)].copy()
        feats = [c for c in sub.columns if c not in ("cls", "group", "chrom")]
        feats = [c for c in feats if sub[c].notna().mean() > 0.8]
        Xs = sub[feats].fillna(sub[feats].median(numeric_only=True))
        y = sub.cls.values
        groups = sub.chrom.astype(str).values
        clf = HistGradientBoostingClassifier(max_depth=6, learning_rate=0.08,
                                             max_iter=400, random_state=SEED)
        gkf = GroupKFold(n_splits=5)
        yhat = np.empty(len(y), dtype=object)
        for tr, te in gkf.split(Xs, y, groups):
            clf.fit(Xs.iloc[tr], y[tr])
            yhat[te] = clf.predict(Xs.iloc[te])
        acc = float((y == yhat).mean())
        macro = f1_score(y, yhat, average="macro")
        bal = balanced_accuracy_score(y, yhat)
        kappa = cohen_kappa_score(y, yhat)
        rows = [dict(metric="top1_accuracy", value=acc),
                dict(metric="macro_f1", value=macro),
                dict(metric="balanced_accuracy", value=bal),
                dict(metric="cohens_kappa", value=kappa),
                dict(metric="n_loci", value=len(y)),
                dict(metric="folds", value=5)]
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "E_ml_chrom_holdout.tsv")), sep="\t", index=False)
        f1s = f1_score(y, yhat, average=None, labels=CLASSES)
        P = pd.DataFrame(dict(cls=CLASSES, f1=f1s))
        P.to_csv(safe_write_path(os.path.join(
            OUT, "E_ml_chrom_per_class.tsv")), sep="\t", index=False)
        st.note("chromosome-held-out: acc={:.3f} macro-F1={:.3f} "
                "balanced={:.3f} kappa={:.3f} n={}".format(
                    acc, macro, bal, kappa, len(y)))
        st.summary_md(os.path.join(
            RES, "step_summaries", "81_ml_chrom_holdout.md"),
            "81 chromosome-held-out classifier",
            df_to_md(T, ".3f") + "\n\n" + df_to_md(P, ".3f"))


if __name__ == "__main__":
    main()
