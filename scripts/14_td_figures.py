#!/usr/bin/env python
"""14_td_figures.py — Fig7 (GTEx tissues), Fig8 (TCGA disease), Fig9 (integration).
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, WORKDIR
from pylib.plotting import CLASS_ORDER, CLASS_COLORS, CLASS_LABELS, save, plt

OUT = os.path.join(WORKDIR, "results", "tissue_disease", "tables")
FIGS = os.path.join(WORKDIR, "figures")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def fig7():
    p = os.path.join(OUT, "gtex_class_by_tissue.tsv")
    if not os.path.exists(p):
        return
    t = pd.read_csv(p, sep="\t").sort_values("n_lnc_expressed", ascending=False)
    top = t.head(25)
    fig, axes = plt.subplots(1, 2, figsize=(13, 8),
                             gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    y = np.arange(len(top))
    left = np.zeros(len(top))
    shares = [top[f"share_{c}"].values for c in CLASSES]
    colors = [CLASS_COLORS[c] for c in CLASSES]
    labels = [CLASS_LABELS[c] for c in CLASSES]
    for s, c, lab in zip(shares, colors, labels):
        ax.barh(y, s, left=left, color=c, label=lab, height=0.8)
        left += s
    ax.set_yticks(y)
    ax.set_yticklabels([s[:28] for s in top["tissue"]], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("share of expressed lncRNAs")
    ax.set_title("A  Birth-class composition of expressed lncRNAs (GTEx, top tissues)")
    ax.legend(fontsize=7, loc="lower right")
    ax2 = axes[1]
    ax2.scatter(top["divergent_frac_expr"], top["intergenic_frac_expr"],
                s=18, color="#888")
    for _, r in top.iterrows():
        if r["divergent_frac_expr"] > 0.5 or r["intergenic_frac_expr"] > 0.5:
            ax2.annotate(r["tissue"][:14], (r["divergent_frac_expr"],
                                            r["intergenic_frac_expr"]),
                         fontsize=6)
    ax2.plot([0, 1], [0, 1], "--", color="#ccc", lw=0.8)
    ax2.set_xlabel("divergent: fraction expressed")
    ax2.set_ylabel("intergenic: fraction expressed")
    ax2.set_title("B  Divergent vs intergenic usage per tissue")
    save(fig, "Fig7_gtex_tissues")


def fig8():
    p = os.path.join(OUT, "tcga_ltr_derepression.tsv")
    p2 = os.path.join(OUT, "tcga_class_de_per_cancer.tsv")
    if not os.path.exists(p2):
        return
    de = pd.read_csv(p2, sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(13, 8),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    # A: heatmap of frac_up (fraction of class genes significantly up in tumor)
    piv = de.pivot(index="cancer", columns="cls", values="frac_up")
    piv = piv[[c for c in CLASSES if c in piv.columns]]
    piv = piv.loc[piv.mean(axis=1).sort_values(ascending=False).index]
    im = axes[0].imshow(piv.values, aspect="auto", cmap="RdBu_r", vmin=0, vmax=0.2)
    axes[0].set_xticks(range(piv.shape[1]))
    axes[0].set_xticklabels([CLASS_LABELS[c].replace(" lncRNA", "") for c in piv.columns],
                            rotation=30, ha="right", fontsize=8)
    axes[0].set_yticks(range(len(piv)))
    axes[0].set_yticklabels([c.replace("TCGA-", "") for c in piv.index],
                            fontsize=7)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if np.isfinite(v) and (v > 0.08 or v < 0.005):
                axes[0].text(j, i, f"{v:.2f}", ha="center", va="center",
                             fontsize=5.5, color="k" if v < 0.12 else "w")
    plt.colorbar(im, ax=axes[0], shrink=0.7,
                 label="fraction of class genes up-regulated (FDR<0.05, |Δ|>1)")
    axes[0].set_title("A  Class-level tumour up-regulation per cancer")
    if os.path.exists(p):
        lt = pd.read_csv(p, sep="\t").sort_values("diff")
        axes[1].barh(np.arange(len(lt)), lt["diff"],
                     color=["#C44E52" if v > 0 else "#55A868" for v in lt["diff"]])
        axes[1].set_yticks(np.arange(len(lt)))
        axes[1].set_yticklabels([c.replace("TCGA-", "") for c in lt["cancer"]],
                                fontsize=7)
        axes[1].axvline(0, color="k", lw=0.8)
        axes[1].set_xlabel("LTR-lncRNA median Δ − other-lncRNA median Δ (log2TPM)")
        axes[1].set_title("B  LTR-derived lncRNA derepression")
    save(fig, "Fig8_tcga_disease")


def fig9():
    p = os.path.join(OUT, "summary_tissue_disease.tsv")
    if not os.path.exists(p):
        return
    t = pd.read_csv(p, sep="\t")
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    sub = t[t.class_or_group.isin(CLASSES)]
    x = sub["gtex_median_frac_expr"].values
    yv = sub["median_frac_up"].values
    ax.scatter(x, yv, s=60, c=[CLASS_COLORS[c] for c in sub.class_or_group])
    for _, r in sub.iterrows():
        ax.annotate(r.class_or_group, (r.gtex_median_frac_expr, r.median_frac_up),
                    fontsize=8, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("GTEx median fraction of tissue-expressed (normal activity)")
    ax.set_ylabel("TCGA median fraction upregulated in tumors")
    ax.set_title("Normal-tissue activity vs tumor dysregulation by birth class")
    save(fig, "Fig9_integration")


def main():
    with step("14_td_figures") as st:
        os.makedirs(FIGS, exist_ok=True)
        fig7()
        fig8()
        fig9()
        st.note("tissue/disease figures done")


if __name__ == "__main__":
    main()
