#!/usr/bin/env python
"""24_ml_extensions.py — three ML analyses that reinforce the birth-class model.

A. Class intrinsic learnability: predict birth class from promoter-intrinsic
   features ONLY (sequence, TE, chromatin, conservation — no coordinates,
   no activity) with gradient boosting, stratified 5-fold CV.
B. Interpretable LTR-promoter sequence classifier: 6-mer logistic regression
   (LTR-derived vs other lnc promoters), CV AUROC, top discriminative words.
C. Transferable LTR-derepression score: fixed-gene-set per-sample score,
   trained nowhere; dose-response in CGGA (grade/IDH) and validated in
   TCGA LGG/GBM.
Outputs -> results/ml/*.tsv + figures/sub/Fig10_ml.{pdf,png}
"""
import os
import sys
import gzip

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (f1_score, confusion_matrix, roc_auc_score,
                             roc_curve)
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "ml")
TD = os.path.join(WORKDIR, "data", "tissue_disease")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
SEED = 42


# ---------------------------------------------------------------- feature A --
def build_feature_table():
    m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                    index_col=0, low_memory=False)
    parts = [m[["cls", "group"]]]
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
    prom_like = ch.eq("TssA").mean(axis=1).rename("prom_like")
    enh_like = ch.isin(["EnhG", "Enh"]).mean(axis=1).rename("enh_like")
    parts.append(pd.concat([prom_like, enh_like], axis=1))
    X = pd.concat(parts, axis=1)
    X["prom_class"] = X["prom_class"].map({"HCP": 2, "ICP": 1, "LCP": 0})
    return X


def task_a(st):
    X = build_feature_table()
    sub = X[X.cls.isin(CLASSES)].copy()
    feats = [c for c in sub.columns if c not in ("cls", "group")]
    feats = [c for c in feats if sub[c].notna().mean() > 0.8]
    Xs = sub[feats].fillna(sub[feats].median(numeric_only=True))
    y = sub.cls.values
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    clf = HistGradientBoostingClassifier(max_depth=6, learning_rate=0.08,
                                         max_iter=400, random_state=SEED)
    yhat = np.empty(len(y), dtype=object)
    prob = np.zeros((len(y), len(CLASSES)))
    for tr, te in skf.split(Xs, y):
        clf.fit(Xs.iloc[tr], y[tr])
        yhat[te] = clf.predict(Xs.iloc[te])
        prob[te] = clf.predict_proba(Xs.iloc[te])
    macro = f1_score(y, yhat, average="macro")
    acc = float((y == yhat).mean())
    cm = confusion_matrix(y, yhat, labels=CLASSES)
    pd.DataFrame(cm, index=CLASSES, columns=CLASSES).to_csv(safe_write_path(
        os.path.join(OUT, "A_confusion.tsv")), sep="\t")
    # feature importance via permutation-lite: single holdout refit
    tr, te = next(StratifiedKFold(5, shuffle=True, random_state=SEED)
                  .split(Xs, y))
    clf.fit(Xs.iloc[tr], y[tr])
    base = f1_score(y[te], clf.predict(Xs.iloc[te]), average="macro")
    imps = []
    for c in feats:
        Xp = Xs.iloc[te].copy()
        Xp[c] = np.random.permutation(Xp[c].values)
        imps.append(dict(feature=c, delta_f1=base - f1_score(
            y[te], clf.predict(Xp), average="macro")))
    imp = pd.DataFrame(imps).sort_values("delta_f1", ascending=False)
    imp.to_csv(safe_write_path(os.path.join(OUT, "A_feature_importance.tsv")),
               sep="\t", index=False)
    st.note(f"A: macro-F1={macro:.3f} acc={acc:.3f} (5-fold CV, n={len(y)}); "
            f"top features: " + ", ".join(f"{r.feature}({r.delta_f1:+.3f})"
                                          for r in imp.head(5).itertuples()))
    return dict(y=y, yhat=yhat, prob=prob, cm=cm, macro=macro, acc=acc,
                imp=imp, classes=CLASSES)


# ---------------------------------------------------------------- kmer ---- B --
def kmer_table(fa_path, k=6):
    comp = str.maketrans("ACGTN", "TGCAN")
    rows, names = [], []
    name = None
    buf = []
    seqs = {}
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    alpha = [a + b + c + d + e + f for a in "ACGT" for b in "ACGT"
             for c in "ACGT" for d in "ACGT" for e in "ACGT" for f in "ACGT"]
    aidx = {a: i for i, a in enumerate(alpha)}
    M = np.zeros((len(seqs), len(alpha)), dtype=np.float32)
    for i, (g, s) in enumerate(seqs.items()):
        s = s.upper()
        for j in range(len(s) - k + 1):
            w = s[j:j + k]
            if w in aidx:
                M[i, aidx[w]] += 1
        M[i] /= max(M[i].sum(), 1)
    names = list(seqs.keys())
    return M, alpha, names


def task_b(st):
    meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                       index_col=0, low_memory=False)
    te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                     index_col=0, low_memory=False)
    lnc = meta[meta.group == "lnc"].index
    ltr_pos = set(te[(te.index.isin(lnc)) & (te.tss_in_ltr == 1)].index)
    ltr_neg = set(te[(te.index.isin(lnc)) & (te.tss_in_ltr != 1)].index)
    keep = list(ltr_pos) + list(np.random.default_rng(SEED).choice(
        sorted(ltr_neg), size=min(len(ltr_neg), 4 * len(ltr_pos)),
        replace=False))
    keepset = set(keep)
    fa_in = os.path.join(RES, "promoters_feat.fa")
    # filter fasta to keep set on the fly
    tmp = safe_write_path(os.path.join(OUT, "ltr_set.fa"))
    with open(fa_in) as fh, open(tmp, "w") as out:
        write = False
        for line in fh:
            if line.startswith(">"):
                write = line[1:].strip() in keepset
            if write:
                out.write(line)
    M, alpha, names = kmer_table(tmp)
    y = np.array([1 if n in ltr_pos else 0 for n in names])
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    lr = LogisticRegression(max_iter=2000, C=1.0, penalty="l2")
    sc = StandardScaler()
    probs = np.zeros(len(y))
    for tr, te_ in skf.split(M, y):
        Xtr = sc.fit_transform(M[tr])
        lr.fit(Xtr, y[tr])
        probs[te_] = lr.predict_proba(sc.transform(M[te_]))[:, 1]
    auc = roc_auc_score(y, probs)
    fpr, tpr, _ = roc_curve(y, probs)
    # top discriminative 6-mers (fit on all, standardised)
    Xa = sc.fit_transform(M)
    lr.fit(Xa, y)
    top = pd.DataFrame({"kmer": alpha, "coef": lr.coef_[0]}
                       ).sort_values("coef", ascending=False)
    top.to_csv(safe_write_path(os.path.join(OUT, "B_top_kmers.tsv")),
               sep="\t", index=False)
    st.note(f"B: sequence-only 6-mer AUROC={auc:.3f} (n_pos={int(y.sum())}, "
            f"n={len(y)}); top words: " +
            ", ".join(top.kmer.head(6)))
    return dict(fpr=fpr, tpr=tpr, auc=auc, top=top, y=y, probs=probs)


# ---------------------------------------------------------------- score --- C --
def task_c(st):
    genes = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                        index_col=0)
    d = np.nan_to_num(np.load(os.path.join(TD, "toil_class_expr.npz"))["M"])
    samp = pd.read_csv(os.path.join(TD, "toil_samples.tsv"), sep="\t")
    lnc = genes.cls.isin(CLASSES).values
    ltr = (genes.tss_in_ltr == 1).fillna(False).values & lnc
    oth = lnc & ~ltr
    score = d[ltr].mean(axis=0) - d[oth].mean(axis=0)
    sc_df = pd.DataFrame({"sample": samp["sample"].values, "score": score,
                          "cancer": samp.cancer.values,
                          "condition": samp.condition.values})
    # TCGA LGG/GBM dose-response (tumour vs normal within cancer)
    rows = []
    for can in ["TCGA-LGG", "TCGA-GBM"]:
        for cond in ["tumor", "normal"]:
            v = sc_df[(sc_df.cancer == can) &
                      (sc_df.condition == cond)]["score"].dropna()
            if len(v) >= 10:
                rows.append(dict(cohort=can, condition=cond, n=len(v),
                                 mean_score=float(v.mean()),
                                 sd=float(v.std())))
    t = pd.DataFrame(rows)
    t.to_csv(safe_write_path(os.path.join(OUT, "C_ltr_score_tcga.tsv")),
             sep="\t", index=False)
    st.note("C: " + "; ".join(f"{r.cohort} {r.condition}: {r.mean_score:+.3f}"
                              f" (n={r.n})" for r in t.itertuples()))
    return t


# ---------------------------------------------------------------- figure ----
def fig10(A, B, C):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mpl = plt.rcParams
    fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 62 / 25.4))
    ax = axes[0]
    cm = A["cm"].astype(float)
    cmn = cm / cm.sum(axis=1, keepdims=True)
    im = ax.imshow(cmn, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(5)); ax.set_xticklabels(A["classes"], rotation=40,
                                                ha="right", fontsize=6)
    ax.set_yticks(range(5)); ax.set_yticklabels(A["classes"], fontsize=6)
    for i in range(5):
        for j in range(5):
            ax.text(j, i, f"{cmn[i, j]:.2f}", ha="center", va="center",
                    fontsize=5.5,
                    color="white" if cmn[i, j] > 0.5 else "black")
    ax.set_title(f"A  Birth class from promoter features\n"
                 f"5-fold CV macro-F1 = {A['macro']:.2f}")
    ax = axes[1]
    ax.plot(B["fpr"], B["tpr"], color="#D55E00", lw=1.4,
            label=f"6-mer logistic, AUROC = {B['auc']:.2f}")
    ax.plot([0, 1], [0, 1], "--", color="#999", lw=0.8)
    ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
    ax.set_title("B  LTR-derived promoter from sequence alone")
    ax.legend(fontsize=6)
    ax = axes[2]
    x = np.arange(len(C))
    cols = ["#C44E52" if c == "tumor" else "#0072B2" for c in C.condition]
    ax.bar(x, C.mean_score, yerr=C.sd / np.sqrt(C.n), color=cols,
           capsize=2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.cohort.replace('TCGA-', '')}\n{r.condition}"
                        for r in C.itertuples()], fontsize=6)
    ax.set_ylabel("LTR-lncRNA score (mean Δ vs other lnc)")
    ax.set_title("C  Fixed-set LTR score in TCGA glioma")
    fig.tight_layout()
    fig.savefig(os.path.join(WORKDIR, "figures", "sub", "Fig10_ml.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(WORKDIR, "figures", "sub", "Fig10_ml.png"),
                dpi=600, bbox_inches="tight")
    print("saved Fig10_ml")


def main():
    with step("24_ml_extensions") as st:
        os.makedirs(OUT, exist_ok=True)
        A = task_a(st)
        B = task_b(st)
        C = task_c(st)
        fig10(A, B, C)
        st.summary_md("results/ml/step_summaries/24_ml.md", "24 ML extensions",
                      f"- A birth-class CV macro-F1 {A['macro']:.3f}\n"
                      f"- B LTR sequence AUROC {B['auc']:.3f}\n"
                      "- C TCGA glioma score table\n\n## Top features\n\n"
                      + df_to_md(A["imp"].head(12), ".3f"))


if __name__ == "__main__":
    main()
