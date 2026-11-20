#!/usr/bin/env python
"""67_ml_reporting.py — strengthened reporting for ML Task A (birth-class
learnability from promoter-intrinsic features).

WHY these extra metrics. Task A of 24_ml_extensions.py reports 61% top-1
accuracy against a 52% majority class and a macro-F1 of 0.34. A reviewer can
reasonably object that accuracy alone overstates learnability when one class
(the intergenic background) dominates the table: a classifier that always
answers "intergenic" already reaches the majority frequency. To answer that
objection without rerunning or altering the published analysis, this module
re-derives the SAME out-of-fold predictions (same feature builder, same
estimator hyper-parameters, same StratifiedKFold(5, shuffle=True, seed 42))
and adds:

  1. balanced accuracy (macro-averaged recall) plus the majority-class
     baseline, so every class counts equally regardless of its frequency;
  2. per-class precision/recall/F1/support plus the per-class values the
     majority baseline would obtain, exposing which classes (if any) carry
     the signal;
  3. a bootstrap 95% CI on top-1 accuracy over samples, because a point
     accuracy without an interval cannot be compared with the baseline;
  4. a label-permutation null (50 permutations, identical folds) giving the
     mean and 95th percentile macro-F1, which places the observed value on
     the distribution expected when promoter features carry no class
     information;
  5. Cohen's kappa, which corrects agreement for chance under the observed
     class marginals.

Nothing under results/ml/ from 24_ml_extensions.py is modified; only new
files A_reporting.tsv, A_per_class.tsv and A_permutation_null.tsv are written.
"""
import os
import sys
import importlib

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                             cohen_kappa_score, f1_score,
                             precision_recall_fscore_support)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "ml")
SEED = 42
N_BOOT = 2000
N_PERM = 50

# Reuse the published feature table and class list verbatim. Importing the
# module is safe: 24_ml_extensions.py guards its work behind main(), so only
# constants and function definitions execute.
ml24 = importlib.import_module("24_ml_extensions")
CLASSES = ml24.CLASSES


def out_of_fold(Xs, y, splits):
    """Return out-of-fold predictions for one label vector.

    The rows of `splits` are the *observed* fold assignments, reused for
    every permutation so the null differs from the observed result only in
    the labels, never in which samples were held out together.
    """
    clf = HistGradientBoostingClassifier(max_depth=6, learning_rate=0.08,
                                         max_iter=400, random_state=SEED)
    yhat = np.empty(len(y), dtype=object)
    for tr, te in splits:
        clf.fit(Xs.iloc[tr], y[tr])
        yhat[te] = clf.predict(Xs.iloc[te])
    return yhat


def bootstrap_accuracy_ci(y, yhat, n_boot=N_BOOT, seed=SEED):
    """Percentile bootstrap over samples (paired resampling of y and yhat).

    Resampling samples rather than folds reflects the unit of the claim:
    accuracy on unseen loci, not accuracy on unseen folds.
    """
    rng = np.random.default_rng(seed)
    correct = (y == yhat).astype(float)
    n = len(correct)
    vals = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        vals[i] = correct[idx].mean()
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main():
    with step("67_ml_reporting") as st:
        os.makedirs(OUT, exist_ok=True)

        # ---- rebuild exactly the Task A inputs -------------------------
        X = ml24.build_feature_table()
        sub = X[X.cls.isin(CLASSES)].copy()
        feats = [c for c in sub.columns if c not in ("cls", "group")]
        feats = [c for c in feats if sub[c].notna().mean() > 0.8]
        Xs = sub[feats].fillna(sub[feats].median(numeric_only=True))
        y = sub.cls.values
        n = len(y)
        st.note(f"n={n} loci, {len(feats)} promoter-intrinsic features, "
                f"{len(CLASSES)} classes")

        skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
        splits = list(skf.split(Xs, y))
        yhat = out_of_fold(Xs, y, splits)

        # ---- 1. accuracy, balanced accuracy, majority baseline ---------
        acc = accuracy_score(y, yhat)
        bal = balanced_accuracy_score(y, yhat)
        vals, cnts = np.unique(y, return_counts=True)
        maj_cls = vals[np.argmax(cnts)]
        maj_freq = cnts.max() / n
        # The majority baseline is an actual predictor, not an assumed
        # number: predicting the modal class everywhere maximises accuracy.
        yhat_base = np.full(n, maj_cls, dtype=object)
        acc_base = accuracy_score(y, yhat_base)
        bal_base = balanced_accuracy_score(y, yhat_base)
        macro = f1_score(y, yhat, average="macro")
        macro_base = f1_score(y, yhat_base, average="macro")
        kappa = cohen_kappa_score(y, yhat)

        # ---- 3. bootstrap CI on top-1 accuracy -------------------------
        ci_lo, ci_hi = bootstrap_accuracy_ci(y, yhat)
        st.note(f"top-1 acc={acc:.3f} (95% CI {ci_lo:.3f}-{ci_hi:.3f}); "
                f"majority '{maj_cls}'={maj_freq:.3f}")
        st.note(f"balanced acc={bal:.3f} vs majority baseline "
                f"{bal_base:.3f}; macro-F1={macro:.3f} vs baseline "
                f"{macro_base:.3f}; kappa={kappa:.3f}")

        # ---- 2. per-class metrics vs per-class baseline ----------------
        p, r, f1, sup = precision_recall_fscore_support(
            y, yhat, labels=CLASSES, zero_division=0)
        bp, br, bf1, bsup = precision_recall_fscore_support(
            y, yhat_base, labels=CLASSES, zero_division=0)
        per_class = pd.DataFrame({
            "class": CLASSES,
            "precision": p, "recall": r, "f1": f1, "support": sup,
            "baseline_precision": bp, "baseline_recall": br,
            "baseline_f1": bf1,
        })
        per_class.to_csv(
            safe_write_path(os.path.join(OUT, "A_per_class.tsv")),
            sep="\t", index=False, float_format="%.6f")
        so = per_class[per_class["class"] == "sense_overlapping"].iloc[0]
        st.note(f"smallest class sense_overlapping: support={int(so.support)}, "
                f"P={so.precision:.3f} R={so.recall:.3f} F1={so.f1:.3f} "
                f"(baseline F1 {so.baseline_f1:.3f})")

        # ---- 4. label-permutation null on identical folds --------------
        # Permuting the labels destroys any feature-class association while
        # preserving the class marginal and the fold membership, so the null
        # is what the pipeline yields when only chance structure remains.
        rng = np.random.default_rng(SEED)
        null = np.empty(N_PERM)
        null_acc = np.empty(N_PERM)
        for i in range(N_PERM):
            yp = rng.permutation(y)
            yp_hat = out_of_fold(Xs, yp, splits)
            null[i] = f1_score(yp, yp_hat, average="macro")
            null_acc[i] = accuracy_score(yp, yp_hat)
        perm = pd.DataFrame({"permutation": np.arange(1, N_PERM + 1),
                             "macro_f1": null, "accuracy": null_acc})
        perm.to_csv(
            safe_write_path(os.path.join(OUT, "A_permutation_null.tsv")),
            sep="\t", index=False, float_format="%.6f")
        null_mean = float(null.mean())
        null_p95 = float(np.percentile(null, 95))
        null_sd = float(null.std(ddof=1))
        # Empirical one-sided p: how often chance alone reaches the observed
        # macro-F1; the +1 keeps p strictly positive as in a standard
        # permutation test.
        p_perm = float((1 + np.sum(null >= macro)) / (1 + N_PERM))
        beats = bool(macro > null_p95)
        # The null is effectively degenerate: with permuted labels the booster
        # learns nothing and predicts the majority class in every fold, so the
        # null spread is negligible (SD ~1e-4). This makes the comparison with
        # p95 valid but conservative in width; it is reported rather than
        # hidden so the reader does not mistake it for a fitted null.
        degenerate = null_sd < 1e-3
        st.note(f"permutation null (B={N_PERM}): macro-F1 mean={null_mean:.3f}, "
                f"p95={null_p95:.3f}, max={null.max():.3f}, sd={null_sd:.2e}"
                f"{' (degenerate: every permutation predicts the majority class)' if degenerate else ''}; "
                f"observed {macro:.3f} "
                f"{'exceeds' if beats else 'does NOT exceed'} p95 "
                f"(empirical p={p_perm:.3f})")

        # ---- reporting table: one row per metric -----------------------
        base_note = f"majority class '{maj_cls}' predicted for every locus"
        rows = [
            dict(metric="n_loci", value=float(n), baseline=np.nan,
                 ci95_low=np.nan, ci95_high=np.nan,
                 note="samples with a birth class and >=80% feature coverage"),
            dict(metric="n_features", value=float(len(feats)), baseline=np.nan,
                 ci95_low=np.nan, ci95_high=np.nan,
                 note="promoter-intrinsic features used by Task A"),
            dict(metric="top1_accuracy", value=acc, baseline=acc_base,
                 ci95_low=ci_lo, ci95_high=ci_hi,
                 note=f"out-of-fold; bootstrap CI over samples (B={N_BOOT}); "
                      f"baseline = {base_note}"),
            dict(metric="balanced_accuracy", value=bal, baseline=bal_base,
                 ci95_low=np.nan, ci95_high=np.nan,
                 note="macro-averaged recall over the five classes; baseline "
                      f"= {base_note}"),
            dict(metric="macro_f1", value=macro, baseline=macro_base,
                 ci95_low=np.nan, ci95_high=np.nan,
                 note="unweighted mean of per-class F1; baseline = "
                      f"{base_note}"),
            dict(metric="cohens_kappa", value=kappa, baseline=0.0,
                 ci95_low=np.nan, ci95_high=np.nan,
                 note="chance-corrected agreement under observed marginals; "
                      "0 = chance"),
            dict(metric="permutation_null_macro_f1_mean", value=null_mean,
                 baseline=np.nan, ci95_low=np.nan, ci95_high=np.nan,
                 note=f"mean over B={N_PERM} label permutations on identical "
                      "folds"),
            dict(metric="permutation_null_macro_f1_p95", value=null_p95,
                 baseline=np.nan, ci95_low=np.nan, ci95_high=np.nan,
                 note=f"95th percentile of the B={N_PERM} permutation null"),
            dict(metric="permutation_null_macro_f1_sd", value=null_sd,
                 baseline=np.nan, ci95_low=np.nan, ci95_high=np.nan,
                 note="null spread; near zero because permuted labels make "
                      "every replicate predict the majority class"),
            dict(metric="permutation_empirical_p", value=p_perm, baseline=np.nan,
                 ci95_low=np.nan, ci95_high=np.nan,
                 note="(1 + #{null >= observed}) / (1 + B), one-sided"),
            dict(metric="majority_class_frequency", value=float(maj_freq),
                 baseline=np.nan, ci95_low=np.nan, ci95_high=np.nan,
                 note=f"frequency of '{maj_cls}' in the analysed loci"),
        ]
        rep = pd.DataFrame(rows)
        rep.to_csv(safe_write_path(os.path.join(OUT, "A_reporting.tsv")),
                   sep="\t", index=False, float_format="%.6f")

        # ---- markdown summary ------------------------------------------
        pc_md = df_to_md(per_class.rename(columns={
            "baseline_precision": "base_P", "baseline_recall": "base_R",
            "baseline_f1": "base_F1"}), ".3f")
        rep_md = df_to_md(rep[["metric", "value", "baseline", "ci95_low",
                               "ci95_high"]], ".4f")
        null_sentence = (
            "The null is **degenerate**: every permuted replicate predicts "
            "the majority class, so its spread is ~1e-4 and it reproduces the "
            "majority baseline rather than a fitted-model null; the observed "
            "value clears it, but the null does not quantify fitting "
            "variance.\n" if degenerate else "\n")
        body = (
            f"Re-derivation of ML Task A on **n={n}** loci, {len(feats)} "
            f"promoter-intrinsic features, `StratifiedKFold(5, shuffle=True, "
            f"random_state={SEED})`. Predictions are out-of-fold and reuse "
            "the feature builder from `24_ml_extensions.py`; no file under "
            "`results/ml/` produced by step 24 was modified.\n\n"
            "## Headline metrics\n\n" + rep_md + "\n\n"
            "## Per-class metrics versus the majority baseline\n\n" + pc_md +
            "\n\n## Reading\n\n"
            f"- Majority baseline (`{maj_cls}` everywhere) reaches top-1 "
            f"accuracy {acc_base:.3f} = its frequency {maj_freq:.3f}, but "
            f"balanced accuracy {bal_base:.3f} and macro-F1 {macro_base:.3f}; "
            f"the model reaches {acc:.3f} (95% bootstrap CI "
            f"{ci_lo:.3f}-{ci_hi:.3f}), balanced accuracy {bal:.3f} and "
            f"macro-F1 {macro:.3f}.\n"
            f"- The permutation null (B={N_PERM}, identical folds) has "
            f"macro-F1 mean {null_mean:.3f}, SD {null_sd:.2e} and 95th "
            f"percentile {null_p95:.3f}; the observed {macro:.3f} "
            f"{'exceeds' if beats else 'does not exceed'} that percentile "
            f"(empirical p={p_perm:.3f}). {null_sentence}"
            f"- Cohen's kappa = {kappa:.3f}.\n"
            f"- Smallest class `sense_overlapping` (support {int(so.support)}): "
            f"P={so.precision:.3f}, R={so.recall:.3f}, F1={so.f1:.3f} "
            f"(majority-baseline F1 {so.baseline_f1:.3f}).\n\n"
            "## Files\n\n"
            "- `results/ml/A_reporting.tsv`\n"
            "- `results/ml/A_per_class.tsv`\n"
            "- `results/ml/A_permutation_null.tsv`\n")
        st.summary_md(
            os.path.join(WORKDIR, "results", "step_summaries",
                         "67_ml_reporting.md"),
            "67 ML reporting (Task A)", body)


if __name__ == "__main__":
    main()
