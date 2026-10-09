#!/usr/bin/env python
"""16_submission_figures.py — submission-grade Fig1-9 (Nature Communications).
Vector PDF + 600dpi PNG into figures/sub/. Okabe-Ito palette, 7pt Arial,
panel letters, significance stars from pub_stats (BH-FDR)."""
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrow, Rectangle, FancyArrowPatch
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, WORKDIR

RES = os.path.join(WORKDIR, "results")
FIGS = os.path.join(WORKDIR, "figures", "sub")

# ---- Okabe-Ito class palette (colorblind-safe) ----------------------------
C = {"PCG": "#0072B2", "divergent": "#E69F00", "antisense": "#009E73",
     "intronic": "#D55E00", "sense_overlapping": "#CC79A7",
     "intergenic": "#56B4E9", "pseudogene": "#999999",
     "transcribed_pseudogene": "#666666"}
LBL = {"PCG": "Protein-coding", "divergent": "Divergent",
       "antisense": "Antisense", "intronic": "Intronic",
       "sense_overlapping": "Sense-overlap", "intergenic": "Intergenic",
       "pseudogene": "Pseudogene", "transcribed_pseudogene": "Transcribed pseudo"}
ORDER = ["PCG", "divergent", "antisense", "intronic", "sense_overlapping",
         "intergenic", "pseudogene", "transcribed_pseudogene"]
LNC5 = ["divergent", "antisense", "intronic", "sense_overlapping", "intergenic"]

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "legend.frameon": False,
})

MM = 1 / 25.4


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(os.path.join(FIGS, f"{name}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIGS, f"{name}.png"), dpi=600,
                bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def panel(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top")


def stars(fdr):
    return "***" if fdr < 1e-10 else ("**" if fdr < 1e-3 else
                                      ("*" if fdr < 0.05 else "ns"))


def load(name, **kw):
    return pd.read_csv(os.path.join(RES, name), sep="\t", **kw)


# ============================================================ Fig 1 ==========
def fig1():
    meta = load("gene_master_table.tsv", index_col=0, low_memory=False)
    counts = meta[meta.group == "lnc"]["cls"].value_counts().reindex(LNC5)
    fig = plt.figure(figsize=(180 * MM, 62 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.18)
    ax = fig.add_subplot(gs[0])
    # ---- schematic: genome line + 5 class mini-diagrams
    ax.set_xlim(0, 10); ax.set_ylim(0, 5.6); ax.axis("off")
    rows = [
        ("Divergent", 4.7, [("pcg", 1.2, 2.2, "+"), ("lnc", 1.2, 2.0, "-")]),
        ("Antisense", 3.5, [("pcg", 1.0, 3.4, "+"), ("lnc", 1.6, 2.8, "-")]),
        ("Intronic", 2.3, [("pcg", 0.8, 4.6, "+"), ("lnc", 2.0, 2.4, "+")]),
        ("Sense-overlap", 1.1, [("pcg", 0.9, 4.4, "+"), ("lnc", 2.3, 2.4, "+")]),
        ("Intergenic", -0.1, [("pcg", 0.6, 1.8, "+"), ("pcg", 7.5, 1.8, "+"),
                              ("lnc", 4.4, 2.2, "-")]),
    ]
    for name, y0, blocks in rows:
        ax.text(0.05, y0 + 1.05, name, fontsize=7, fontweight="bold",
                color="#333")
        # genomic line
        ax.plot([0.4, 9.6], [y0 + 0.25, y0 + 0.25], color="#bbb", lw=1.2,
                zorder=1)
        for kind, x, w, strand in blocks:
            col = C["PCG"] if kind == "pcg" else C[
                {"Divergent": "divergent", "Antisense": "antisense",
                 "Intronic": "intronic", "Sense-overlap": "sense_overlapping",
                 "Intergenic": "intergenic"}[name]]
            yb = y0 + 0.25 if kind == "pcg" else y0 - 0.18
            for seg in range(3):
                sx = x + seg * (w / 3 + 0.12)
                sw = w / 3
                rect = Rectangle((sx, yb - 0.11), sw, 0.22, facecolor=col,
                                 edgecolor="none", zorder=3, alpha=0.9)
                ax.add_patch(rect)
            # arrow points downstream of the TSS
            arr_end = x + w + 0.35 if strand == "+" else x - 0.35
            ax.annotate("", xy=(arr_end, yb), xytext=(x + w / 2, yb),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0))
            # TSS dot
            ax.plot([x if strand == "+" else x + w], [yb], marker="|",
                    color="k", ms=4, mew=0.8, zorder=4)
    ax.text(0.4, 5.55, "PCG (blue)  ·  lncRNA (coloured)  ·  | = TSS",
            fontsize=6, color="#555")
    panel(ax, "a")
    ax2 = fig.add_subplot(gs[1])
    y = np.arange(len(counts))[::-1]
    ax2.barh(y, counts.values, color=[C[c] for c in counts.index],
             height=0.65)
    for yi, v in zip(y, counts.values):
        ax2.text(v + 120, yi, f"{v:,}", va="center", fontsize=6.5)
    ax2.set_yticks(y); ax2.set_yticklabels([LBL[c] for c in counts.index])
    ax2.set_xlabel("lncRNA genes (n = 19,258)")
    ax2.set_xlim(0, 11500)
    panel(ax2, "b")
    save(fig, "Fig1_classes")


# ============================================================ Fig 2 ==========
def fig2():
    by = load("seq_features_by_class.tsv").set_index("cls")
    pub = load("pub_stats.tsv")
    feats = [("gc", "GC content", "GC content"),
             ("cpg_oe", "CpG obs/exp", "CpG obs/exp"),
             ("HCP", "High-CpG promoter fraction", "HCP promoter"),
             ("TATA", "TATA-box fraction", "TATA box")]
    fig, axes = plt.subplots(1, 4, figsize=(183 * MM, 42 * MM))
    for k, (col, lab, feat) in enumerate(feats):
        ax = axes[k]
        vals = [by.loc[c, col] for c in ORDER]
        ax.bar(np.arange(len(ORDER)), vals, color=[C[c] for c in ORDER],
               width=0.72)
        ax.set_xticks(np.arange(len(ORDER)))
        ax.set_xticklabels([LBL[c].replace("Protein-coding", "PCG")
                            .replace("Transcribed pseudo", "tPseudo") for c in ORDER],
                           rotation=42, ha="right")
        ax.set_ylabel(lab)
        # stars vs PCG
        ref = by.loc["PCG", col]
        for i, c in enumerate(ORDER):
            if c == "PCG":
                continue
            r = pub[(pub.feature == feat) & (pub.cls == c)]
            if len(r):
                ax.text(i, vals[i], stars(r.fdr.iloc[0]), ha="center",
                        va="bottom", fontsize=5.5)
        ax.axhline(ref, color="#555", lw=0.5, ls=":")
        panel(ax, "abcd"[k])
    fig.tight_layout(w_pad=1.4)
    save(fig, "Fig2_sequence")


# ============================================================ Fig 3 ==========
def fig3():
    by = load("te_by_class.tsv").set_index("cls")
    age = load("te_age_by_class.tsv").set_index("cls")
    fig, axes = plt.subplots(1, 3, figsize=(160 * MM, 42 * MM))
    for k, (col, lab, feat) in enumerate([
            ("te_cov", "TE coverage of promoter", None),
            ("tss_in_ltr", "TSS inside LTR/ERV", "TSS inside LTR/ERV"),
            ("tss_in_te", "TSS inside any TE", None)]):
        ax = axes[k]
        vals = [by.loc[c, col] for c in ORDER]
        ax.bar(np.arange(len(ORDER)), vals, color=[C[c] for c in ORDER],
               width=0.72)
        ax.set_xticks(np.arange(len(ORDER)))
        ax.set_xticklabels([LBL[c].replace("Protein-coding", "PCG")
                            .replace("Transcribed pseudo", "tPseudo") for c in ORDER],
                           rotation=42, ha="right")
        ax.set_ylabel(lab)
        if feat:
            r = pub = load("pub_stats.tsv")
            for i, c in enumerate(ORDER):
                if c == "PCG":
                    continue
                rr = pub[(pub.feature == feat) & (pub.cls == c)]
                if len(rr):
                    ax.text(i, vals[i], stars(rr.fdr.iloc[0]), ha="center",
                            va="bottom", fontsize=5.5)
        panel(ax, "abc"[k])
    # age stacked (right panel replaces panel c? keep separate small)
    panel(axes[2], "c")
    save(fig, "Fig3_TE")


def fig3b():
    age = load("te_age_by_class.tsv").set_index("cls")
    order = ["PCG", "divergent", "antisense", "intronic", "intergenic",
             "pseudogene"]
    fig, ax = plt.subplots(figsize=(80 * MM, 45 * MM))
    young = [age.loc[c, "frac_young"] for c in order]
    old = [age.loc[c, "old"] / age.loc[c, "n_TSS_in_TE"] for c in order]
    anc = [age.loc[c, "frac_ancient"] for c in order]
    x = np.arange(len(order))
    ax.bar(x, young, color="#D55E00", label="primate-young")
    ax.bar(x, old, bottom=young, color="#E69F00", label="primate-old")
    ax.bar(x, anc, bottom=np.array(young) + np.array(old), color="#999999",
           label="ancient")
    ax.set_xticks(x)
    ax.set_xticklabels([LBL[c].replace("Protein-coding", "PCG") for c in order],
                       rotation=38, ha="right")
    ax.set_ylabel("fraction of TSS-donating TEs")
    ax.legend(fontsize=6)
    panel(ax, "d")
    save(fig, "Fig3d_TE_age")


# ============================================================ Fig 4 ==========
def fig4():
    by = load("cage_activity_by_class.tsv").set_index("cls")
    fig, axes = plt.subplots(1, 3, figsize=(160 * MM, 42 * MM))
    specs = [("active_frac", "CAGE-active promoters", "CAGE-active promoter"),
             ("log10_maxtpm", "Promoter output (log$_{10}$ max TPM)", None),
             ("tau", "Tissue specificity τ (CAGE)", None)]
    pub = load("pub_stats.tsv")
    for k, (col, lab, feat) in enumerate(specs):
        ax = axes[k]
        vals = [by.loc[c, col] for c in ORDER]
        ax.bar(np.arange(len(ORDER)), vals, color=[C[c] for c in ORDER],
               width=0.72)
        ax.set_xticks(np.arange(len(ORDER)))
        ax.set_xticklabels([LBL[c].replace("Protein-coding", "PCG")
                            .replace("Transcribed pseudo", "tPseudo") for c in ORDER],
                           rotation=42, ha="right")
        ax.set_ylabel(lab)
        if feat:
            for i, c in enumerate(ORDER):
                if c == "PCG":
                    continue
                rr = pub[(pub.feature == feat) & (pub.cls == c)]
                if len(rr):
                    ax.text(i, vals[i], stars(rr.fdr.iloc[0]), ha="center",
                            va="bottom", fontsize=5.5)
        panel(ax, "abc"[k])
    fig.tight_layout(w_pad=1.4)
    save(fig, "Fig4_activity")


def fig4b():
    d = load("divergent_pairs_activity.tsv")
    fig, ax = plt.subplots(figsize=(55 * MM, 45 * MM))
    data = [d["rho_divergent"].dropna(), d["rho_control_same_strand"].dropna()]
    vp = ax.violinplot(data, showmedians=True, widths=0.8)
    for b in vp["bodies"]:
        b.set_facecolor("#0072B2"); b.set_alpha(0.6)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["divergent\npairs",
                                               "same-strand\ncontrols"])
    ax.set_ylabel("Spearman ρ of CAGE activity")
    panel(ax, "d")
    save(fig, "Fig4d_pair_rho")


# ============================================================ Fig 5 ==========
def fig5():
    me = load("motif_enrichment.tsv")
    tf_ltr = load("final_gaps/tf_binding_vs_LTR.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(160 * MM, 88 * MM),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    sig = me[(me.fdr < 1e-5) & (me.log2OR > 0.3)]
    top = (sig.sort_values("log2OR", ascending=False)
           .groupby("cls").head(8))
    top = top.copy()
    top["tf"] = top["tf"].str.replace("_HUMAN", "", regex=False)
    piv = top.pivot_table(index="tf", columns="cls", values="log2OR")
    cols = [c for c in LNC5 if c in piv.columns]
    piv = piv[cols].dropna(how="all")
    piv = piv.loc[piv.max(axis=1).sort_values(ascending=False).index].head(24)
    im = ax.imshow(piv.values, aspect="auto", cmap="OrRd", vmin=0, vmax=1.6)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([LBL[c] for c in cols], rotation=35, ha="right")
    ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=6)
    plt.colorbar(im, ax=ax, shrink=0.75, label="motif log$_2$OR vs PCG",
                 pad=0.02)
    panel(ax, "a")
    ax2 = axes[1]
    t = tf_ltr.sort_values("frac_LTR_among_bound")
    y = np.arange(len(t))
    ax2.barh(y, t.frac_LTR_among_bound, color="#D55E00",
             label="LTR-derived fraction among bound lnc promoters")
    ax2.plot([t.background_LTR_frac.iloc[0]] * 2, [-0.5, len(t) - 0.5],
             ls="--", color="#0072B2", lw=1,
             label=f"background ({t.background_LTR_frac.iloc[0]:.1%})")
    for i, r in enumerate(t.itertuples()):
        ax2.text(r.frac_LTR_among_bound + 0.005, i, stars(r.p), va="center",
                 fontsize=6)
    ax2.set_yticks(y); ax2.set_yticklabels(t.tf)
    ax2.set_xlabel("fraction of TF-bound lncRNA promoters that are LTR-derived")
    ax2.legend(loc="lower right")
    ax2.set_xlim(0, 0.46)
    panel(ax2, "b")
    save(fig, "Fig5_TFs")


# ============================================================ Fig 6 ==========
def fig6():
    pc = load("robustness_partial_corr.tsv")
    ch = load("chromhmm_by_class.tsv").set_index("cls")
    pp = load("phylop/partial_corr_phylop.tsv")
    tp = load("phylop/phastcons_vs_phylop.tsv").set_index("cls")
    fig, axes = plt.subplots(1, 4,
                             figsize=(183 * MM, 50 * MM),
                             gridspec_kw={"width_ratios": [1, 1, 1, 1.15]})
    ax = axes[0]
    x = np.arange(len(pc)); w = 0.36
    ax.bar(x - w / 2, pc.rho_raw, w, color="#0072B2", label="raw")
    ax.bar(x + w / 2, pc.rho_GC_controlled, w, color="#56B4E9",
           label="GC-controlled")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels([LBL[c] for c in pc.cls],
                                         rotation=38, ha="right")
    ax.set_ylabel("Spearman ρ (class vs phastCons)")
    ax.legend()
    panel(ax, "a")
    ax2 = axes[1]
    order = [c for c in ORDER if c in ch.index]
    pc_like = ch[[c for c in ch.columns if c.endswith("_prom_like")]].mean(axis=1)
    en_like = ch[[c for c in ch.columns if c.endswith("_enh_like")]].mean(axis=1)
    x = np.arange(len(order)); w = 0.38
    ax2.bar(x - w / 2, [pc_like[c] for c in order], w, color="#0072B2",
            label="promoter-like")
    ax2.bar(x + w / 2, [en_like[c] for c in order], w, color="#E69F00",
            label="enhancer-like")
    ax2.set_xticks(x)
    ax2.set_xticklabels([LBL[c].replace("Protein-coding", "PCG")
                         .replace("Transcribed pseudo", "tPseudo") for c in order],
                        rotation=42, ha="right")
    ax2.set_ylabel("fraction of TSS (5 cell types)")
    ax2.legend()
    panel(ax2, "b")
    ax3 = axes[2]
    x = np.arange(len(pp)); w = 0.36
    ax3.bar(x - w / 2, pp.rho_raw, w, color="#0072B2", label="raw")
    ax3.bar(x + w / 2, pp.rho_GC_controlled, w, color="#56B4E9",
            label="GC-controlled")
    ax3.axhline(0, color="k", lw=0.6)
    ax3.set_xticks(x); ax3.set_xticklabels([LBL[c] for c in pp.cls],
                                           rotation=38, ha="right")
    ax3.set_ylabel("Spearman ρ (class vs phyloP)")
    ax3.legend()
    panel(ax3, "c")
    ax4 = axes[3]
    order4 = [c for c in ORDER if c in tp.index]
    x = np.arange(len(order4)); w = 0.38
    ax4.bar(x - w / 2, [tp.loc[c, "frac_accel_2"] for c in order4], w,
            color=[C[c] for c in order4], label="all promoters")
    ltr_x = [xi + w / 2 for xi, c in zip(x, order4)
             if c in LNC5 and not np.isnan(tp.loc[c, "frac_accel2_ltr"])]
    ltr_y = [tp.loc[c, "frac_accel2_ltr"] for c in order4
             if c in LNC5 and not np.isnan(tp.loc[c, "frac_accel2_ltr"])]
    ax4.bar(ltr_x, ltr_y, w, color="#BBBBBB", hatch="//",
            label="LTR-derived TSS")
    ax4.set_xticks(x)
    ax4.set_xticklabels([LBL[c].replace("Protein-coding", "PCG")
                         .replace("Transcribed pseudo", "tPseudo") for c in order4],
                        rotation=42, ha="right")
    ax4.set_ylabel("fraction of bases phyloP ≤ −2")
    ax4.legend()
    panel(ax4, "d")
    fig.tight_layout(w_pad=1.4)
    save(fig, "Fig6_conservation_chromatin")


# ============================================================ Fig 7 ==========
def fig7():
    t = load("tissue_disease/tables/gtex_class_by_tissue.tsv",
             ).sort_values("n_lnc_expressed", ascending=False).head(25)
    fig, ax = plt.subplots(figsize=(95 * MM, 95 * MM))
    y = np.arange(len(t))[::-1]
    left = np.zeros(len(t))
    for c in LNC5:
        ax.barh(y, t[f"share_{c}"], left=left, color=C[c],
                label=LBL[c], height=0.78)
        left += t[f"share_{c}"].values
    ax.set_yticks(y); ax.set_yticklabels([s[:26] for s in t.tissue],
                                         fontsize=6)
    ax.set_xlabel("share of expressed lncRNAs (GTEx)")
    ax.legend(fontsize=6, loc="lower right")
    # annotate testis LTR
    ti = [i for i, s in enumerate(t.tissue) if "Testis" in s]
    if ti:
        ax.text(0.985, y[ti[0]], "LTR-rich", ha="right", va="center",
                fontsize=6, color="#D55E00", fontweight="bold")
    panel(ax, "a")
    save(fig, "Fig7_gtex")


# ============================================================ Fig 8 ==========
def fig8():
    de = load("tissue_disease/tables/tcga_class_de_per_cancer.tsv")
    lt = load("tissue_disease/tables/tcga_ltr_derepression.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(150 * MM, 68 * MM),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    piv = de.pivot(index="cancer", columns="cls", values="frac_up")
    cols = [c for c in LNC5 if c in piv.columns]
    piv = piv[cols].loc[piv.mean(axis=1).sort_values(
        ascending=False).index]
    piv = piv.loc[[c for c in piv.index
                   if de[de.cancer == c].n_normal.max() >= 5]]
    im = ax.imshow(piv.values, aspect="auto", cmap="PuBu", vmin=0, vmax=0.2)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([LBL[c] for c in cols], rotation=35, ha="right")
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([c.replace("TCGA-", "") for c in piv.index],
                       fontsize=6)
    plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02,
                 label="fraction up-regulated\n(FDR<0.05, Δ>0.5)")
    panel(ax, "a")
    ax2 = axes[1]
    lt = lt.sort_values("diff")
    ax2.barh(np.arange(len(lt)), lt["diff"],
             color=["#D55E00" if v > 0 else "#009E73" for v in lt["diff"]])
    ax2.set_yticks(np.arange(len(lt)))
    ax2.set_yticklabels([c.replace("TCGA-", "") for c in lt.cancer],
                        fontsize=6)
    ax2.axvline(0, color="k", lw=0.7)
    ax2.set_xlabel("LTR-derived Δ − other-lncRNA Δ (log$_2$TPM)")
    panel(ax2, "b")
    # panel c (in-house GBM grade boxplot) removed from manuscript 2026-10-08;
    # cohort archived in results/gbm_own, data intact, restorable anytime
    save(fig, "Fig8_tcga")


# ============================================================ Fig 9 ==========
def fig9():
    iv = load("final_gaps/activity_vs_dysregulation.tsv")
    fg = load("final_gaps/fantom_vs_gtex_tau.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(140 * MM, 55 * MM))
    ax = axes[0]
    sub = iv[iv.cls.isin(LNC5)]
    ax.scatter(sub.gtex_median_frac_expr, sub.tcga_median_frac_up,
               s=42, c=[C[c] for c in sub.cls], zorder=3)
    for _, r in sub.iterrows():
        ax.annotate(LBL[r.cls], (r.gtex_median_frac_expr,
                                 r.tcga_median_frac_up), fontsize=6,
                    xytext=(4, 3), textcoords="offset points")
    ax.set_xlabel("GTEx median fraction of tissues expressing the class")
    ax.set_ylabel("TCGA median fraction up-regulated")
    panel(ax, "a")
    ax2 = axes[1]
    ax2.axis("off")
    ax2.text(0.02, 0.86, "Cross-platform concordance", fontsize=8,
             fontweight="bold")
    ax2.text(0.02, 0.62,
             f"Promoter activity (CAGE) vs body expression\n"
             f"(RNA-seq) tissue specificity, gene level:\n\n"
             f"Spearman ρ = {fg.rho.iloc[0]:.2f},  n = {int(fg.n.iloc[0]):,} genes",
             fontsize=8)
    ax2.text(0.02, 0.32,
             "GTEx class τ medians:\nPCG 0.59  ·  divergent 0.97\n"
             "antisense 0.99  ·  intergenic 1.00  ·  intronic 1.00",
             fontsize=7, color="#333")
    panel(ax2, "b")
    save(fig, "Fig9_integration")


def main():
    with step("16_submission_figures") as st:
        fig1(); fig2(); fig3(); fig3b(); fig4(); fig4b(); fig5(); fig6()
        fig7(); fig8(); fig9()
        st.note(f"submission figures -> {FIGS}")


if __name__ == "__main__":
    main()
