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
     # intronic and pseudogene were reassigned after a colour-vision check.
     # Okabe-Ito orange and vermillion sit at dE 12.8 under deuteranopia,
     # too close for two classes that share most panels, so intronic takes a
     # violet outside the original set and pseudogene a darker grey.
     "intronic": "#5D3A9B", "sense_overlapping": "#CC79A7",
     "intergenic": "#56B4E9", "pseudogene": "#555555",
     "transcribed_pseudogene": "#333333"}
LBL = {"PCG": "Protein-coding", "divergent": "Divergent",
       "antisense": "Antisense", "intronic": "Intronic",
       "sense_overlapping": "Sense-overlap", "intergenic": "Intergenic",
       "pseudogene": "Pseudogene", "transcribed_pseudogene": "Transcribed pseudo"}
ORDER = ["PCG", "divergent", "antisense", "intronic", "sense_overlapping",
         "intergenic", "pseudogene", "transcribed_pseudogene"]
LNC5 = ["divergent", "antisense", "intronic", "sense_overlapping", "intergenic"]

# Arial is required by Nature-family journals. Where it is absent (this
# cluster), Liberation Sans is metric-compatible with Arial, so register it
# explicitly instead of silently falling back to the wider DejaVu Sans.
from matplotlib import font_manager as _fm
for _p in ("/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
           "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf",
           "/usr/share/fonts/liberation-sans/LiberationSans-Italic.ttf",
           "/usr/share/fonts/liberation-sans/LiberationSans-BoldItalic.ttf"):
    if os.path.exists(_p):
        _fm.fontManager.addfont(_p)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "legend.frameon": False,
    # mathtext defaults to DejaVu, which reintroduces a second font in the PDF
    "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Sans",
    "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold",
    "mathtext.cal": "Liberation Sans",
    "mathtext.sf": "Liberation Sans",
})

MM = 1 / 25.4
# Full-width figures use 177 mm rather than the 180 mm journal maximum, because
# the tight bounding box adds a millimetre or two of overhang for axis labels and
# the cropped PDF must still land inside 180.


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(os.path.join(FIGS, f"{name}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIGS, f"{name}.png"), dpi=600,
                bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def panel(ax, letter):
    # Anchored just above the top-left corner of the axes. This clears both the
    # long rotated y labels and the uppermost y tick label in narrow panels.
    ax.text(-0.02, 1.15, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="right")


def stars(fdr):
    return "***" if fdr < 1e-10 else ("**" if fdr < 1e-3 else
                                      ("*" if fdr < 0.05 else "ns"))


def wilson(frac, n, z=1.96):
    """Wilson score interval for a proportion, used for the binomial panels."""
    p = np.asarray(frac, dtype=float)
    n = np.asarray(n, dtype=float)
    den = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return np.clip(centre - half, 0, 1), np.clip(centre + half, 0, 1)


def load(name, **kw):
    return pd.read_csv(os.path.join(RES, name), sep="\t", **kw)


# ============================================================ Fig 1 ==========
def fig1():
    meta = load("gene_master_table.tsv", index_col=0, low_memory=False)
    counts = meta[meta.group == "lnc"]["cls"].value_counts().reindex(LNC5)
    fig = plt.figure(figsize=(177 * MM, 62 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.32)
    ax = fig.add_subplot(gs[0])
    # ---- schematic: genome line + 5 class mini-diagrams
    ax.set_xlim(0, 10); ax.set_ylim(-0.6, 6.1); ax.axis("off")
    # Class names sit on the row baseline at the far left and every block is
    # shifted right, so a label can never fall on the row above it.
    X0 = 1.7
    rows = [
        ("Divergent", 4.7, [("pcg", 1.2, 2.2, "+"), ("lnc", 1.2, 2.0, "-")]),
        ("Antisense", 3.5, [("pcg", 1.0, 3.4, "+"), ("lnc", 1.6, 2.8, "-")]),
        ("Intronic", 2.3, [("pcg", 0.8, 4.6, "+"), ("lnc", 2.0, 2.4, "+")]),
        ("Sense-overlap", 1.1, [("pcg", 0.9, 4.4, "+"), ("lnc", 2.3, 2.4, "+")]),
        ("Intergenic", -0.1, [("pcg", 0.6, 1.8, "+"), ("pcg", 7.5, 1.8, "+"),
                              ("lnc", 4.4, 2.2, "-")]),
    ]
    for name, y0, blocks in rows:
        # centred between the PCG line and the lncRNA block of the same row
        ax.text(0.05, y0 + 0.04, name, fontsize=6.8, fontweight="bold",
                color="#333", va="center")
        # genomic line
        ax.plot([0.4 + X0, 9.1], [y0 + 0.25, y0 + 0.25], color="#bbb",
                lw=1.2, zorder=1)
        for kind, x, w, strand in blocks:
            col = C["PCG"] if kind == "pcg" else C[
                {"Divergent": "divergent", "Antisense": "antisense",
                 "Intronic": "intronic", "Sense-overlap": "sense_overlapping",
                 "Intergenic": "intergenic"}[name]]
            yb = y0 + 0.25 if kind == "pcg" else y0 - 0.18
            for seg in range(3):
                sx = x + X0 + seg * (w / 3 + 0.12)
                sw = w / 3
                rect = Rectangle((sx, yb - 0.11), sw, 0.22, facecolor=col,
                                 edgecolor="none", zorder=3, alpha=0.9)
                ax.add_patch(rect)
            # arrow points downstream of the TSS
            arr_end = (x + w + 0.35 if strand == "+" else x - 0.35) + X0
            ax.annotate("", xy=(arr_end, yb), xytext=(x + w / 2 + X0, yb),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0))
            # TSS dot
            ax.plot([(x if strand == "+" else x + w) + X0], [yb], marker="|",
                    color="k", ms=4, mew=0.8, zorder=4)
    ax.text(0.4, 5.92, "PCG (blue)  ·  lncRNA (coloured)  ·  | = TSS",
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
    sg = load("seq_features.tsv", low_memory=False)
    pub = load("pub_stats.tsv")
    # Two continuous features are shown as distributions and two binary ones as
    # proportions with Wilson intervals. A bar of means would hide the shape of
    # the GC distribution, which is the point the panel has to make.
    feats = [("gc", "GC content", "GC content", "dist"),
             ("cpg_oe", "CpG obs/exp", "CpG obs/exp", "dist"),
             ("HCP", "High-CpG fraction", "HCP promoter", "prop"),
             ("TATA", "TATA-box fraction", "TATA box", "prop")]
    fig, axes = plt.subplots(1, 4, figsize=(177 * MM, 47 * MM))
    for k, (col, lab, feat, kind) in enumerate(feats):
        ax = axes[k]
        xt = [LBL[c].replace("Protein-coding", "PCG")
              .replace("Transcribed pseudo", "tPseudo") for c in ORDER]
        if kind == "dist":
            data, top = [], []
            for c in ORDER:
                v = sg.loc[sg["cls"] == c, col].dropna().astype(float).values
                lo, hi = np.percentile(v, [1, 99])
                data.append(v[(v >= lo) & (v <= hi)])
                top.append(np.percentile(v, 97))
            parts = ax.violinplot(data, positions=np.arange(len(ORDER)),
                                  widths=0.86, showextrema=False,
                                  showmedians=False)
            for body, c in zip(parts["bodies"], ORDER):
                body.set_facecolor(C[c]); body.set_edgecolor("#FFFFFF")
                body.set_linewidth(0.3); body.set_alpha(0.95)
            ax.scatter(np.arange(len(ORDER)), [np.median(d) for d in data],
                       s=4, color="white", zorder=4, linewidths=0)
            top = np.array(top)
        else:
            vals = np.array([by.loc[c, col] for c in ORDER], float)
            ns = np.array([by.loc[c, "n"] for c in ORDER], float)
            lo, hi = wilson(vals, ns)
            for i, c in enumerate(ORDER):
                ax.errorbar(i, vals[i], yerr=[[vals[i] - lo[i]],
                                              [hi[i] - vals[i]]],
                            fmt="o", ms=3.0, color=C[c], ecolor=C[c],
                            elinewidth=0.7, capsize=1.5, zorder=3)
            top = hi
        ax.set_xticks(np.arange(len(ORDER)))
        ax.set_xticklabels(xt, rotation=45, ha="right", fontsize=5.8)
        ax.set_ylabel(lab)
        ax.set_xlim(-0.6, len(ORDER) - 0.4)
        for i, c in enumerate(ORDER):
            if c == "PCG":
                continue
            r = pub[(pub.feature == feat) & (pub.cls == c)]
            if len(r):
                ax.text(i, top[i], stars(r.fdr.iloc[0]), ha="center",
                        va="bottom", fontsize=5.5)
        if kind == "dist":
            ax.axhline(by.loc["PCG", col], color="#555", lw=0.5, ls=":")
        panel(ax, "abcd"[k])
    fig.tight_layout(w_pad=1.4)
    save(fig, "Fig2_sequence")


# ============================================================ Fig 3 ==========
def fig3():
    by = load("te_by_class.tsv").set_index("cls")
    age = load("te_age_by_class.tsv").set_index("cls")
    # wspace keeps each y label out of the panel to its left, which is where
    # it landed at the default spacing.
    fig, axes = plt.subplots(1, 3, figsize=(160 * MM, 42 * MM),
                             gridspec_kw={"wspace": 0.55})
    for k, (col, lab, feat) in enumerate([
            ("te_cov", "TE coverage of promoter", "TE coverage"),
            ("tss_in_ltr", "TSS inside LTR/ERV", "TSS inside LTR/ERV"),
            ("tss_in_te", "TSS inside any TE", "TSS inside any TE")]):
        ax = axes[k]
        vals = np.array([by.loc[c, col] for c in ORDER], float)
        kw = {}
        top = vals.copy()
        if col in ("tss_in_ltr", "tss_in_te"):
            n = np.array([by.loc[c, "n"] for c in ORDER], float)
            lo, hi = wilson(vals, n)
            kw = dict(yerr=np.vstack([vals - lo, hi - vals]), capsize=1.5,
                      error_kw=dict(elinewidth=0.6, ecolor="#333333"))
            top = hi
        ax.bar(np.arange(len(ORDER)), vals, color=[C[c] for c in ORDER],
               width=0.72, **kw)
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
                    ax.text(i, top[i], stars(rr.fdr.iloc[0]), ha="center",
                            va="bottom", fontsize=5.5)
        ax.margins(y=0.16)
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
    # The two primate bins sit next to each other in every bar, so they need to
    # separate by lightness as well as hue. Vermillion and orange differ by only
    # dE 12.8 under deuteranopia, hence the dark burnt orange for the young bin.
    ax.bar(x, young, color="#8C2D04", label="primate-young")
    ax.bar(x, old, bottom=young, color="#E69F00", label="primate-old")
    ax.bar(x, anc, bottom=np.array(young) + np.array(old), color="#999999",
           label="ancient")
    ax.set_xticks(x)
    ax.set_xticklabels([LBL[c].replace("Protein-coding", "PCG") for c in order],
                       rotation=38, ha="right")
    ax.set_ylabel("fraction of TSS-donating TEs")
    # Legend outside the axes: inside, it covered the top of the intergenic bar.
    ax.legend(fontsize=6, loc="lower center", bbox_to_anchor=(0.5, 1.02),
              ncol=3, handlelength=0.9, columnspacing=0.9, handletextpad=0.3)
    panel(ax, "d")
    save(fig, "Fig3d_TE_age")


# ============================================================ Fig 4 ==========
def fig4():
    by = load("cage_activity_by_class.tsv").set_index("cls")
    cg = load("cage_activity.tsv", low_memory=False)
    fig, axes = plt.subplots(1, 3, figsize=(160 * MM, 44 * MM))
    specs = [("active_frac", "CAGE-active promoters", "CAGE-active promoter"),
             ("log10_maxtpm", "Promoter output", "log10 max TPM (promoter output)"),
             ("tau", "Tissue specificity τ", None)]
    pub = load("pub_stats.tsv")
    for k, (col, lab, feat) in enumerate(specs):
        ax = axes[k]
        vals = np.array([by.loc[c, col] for c in ORDER], float)
        ns = np.array([by.loc[c, "n"] for c in ORDER], float)
        if col == "active_frac":
            lo, hi = wilson(vals, ns)
            err = np.vstack([vals - lo, hi - vals])
        else:
            if col == "log10_maxtpm":
                x = np.log10(cg["max_tpm"].fillna(0) + 1)
            else:
                x = cg["tau"]
            sem = np.array([x[cg["cls"] == c].dropna().sem() for c in ORDER])
            err = 1.96 * sem
        ax.bar(np.arange(len(ORDER)), vals, color=[C[c] for c in ORDER],
               width=0.72, yerr=err, capsize=1.8,
               error_kw=dict(elinewidth=0.7, ecolor="#333333"))
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
                    top = vals[i] + (err[1, i] if err.ndim == 2 else err[i])
                    ax.text(i, top, stars(rr.fdr.iloc[0]),
                            ha="center", va="bottom", fontsize=5.5)
        ax.margins(y=0.16)
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
    fig, axes = plt.subplots(1, 2, figsize=(177 * MM, 88 * MM),
                             gridspec_kw={"width_ratios": [1.35, 1],
                                          "wspace": 0.38})
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
    plt.colorbar(im, ax=ax, shrink=0.75, label="Motif log$_2$OR vs PCG",
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
    ax2.set_xlabel("Fraction of bound promoters\nthat are LTR-derived")
    # Direct label instead of a legend, which covered the MAZ bar.
    ax2.text(t.background_LTR_frac.iloc[0], -0.95,
             "background {:.1f}%".format(100 * t.background_LTR_frac.iloc[0]),
             ha="center", va="center", fontsize=5.5, color="#0072B2")
    ax2.set_ylim(-1.35, len(t) - 0.4)
    ax2.set_xlim(0, 0.46)
    panel(ax2, "b")
    save(fig, "Fig5_TFs")


# ============================================================ Fig 6 ==========
def fig6():
    pc = load("robustness_partial_corr.tsv")
    ch = load("chromhmm_by_class.tsv").set_index("cls")
    pp = load("phylop/partial_corr_phylop.tsv")
    tp = load("phylop/phastcons_vs_phylop.tsv").set_index("cls")
    ev = load("comparative/E_transcript_age.tsv", index_col=0)
    sub = load("comparative/G_promoter_substitution.tsv", index_col=0)
    meta6 = load("gene_master_table.tsv", index_col=0,
                 low_memory=False)[["cls", "group"]]
    te6 = load("te_promoter_stats.tsv", index_col=0,
               low_memory=False)[["tss_in_ltr"]]
    ev = ev.join(meta6).join(te6)
    sub = sub.join(meta6).join(te6)
    fig, axes = plt.subplots(2, 3, figsize=(177 * MM, 102 * MM))
    ax = axes[0, 0]
    # Dumbbell: the distance the second point moves is the GC effect, which a
    # pair of bars makes the reader estimate instead of see.
    y = np.arange(len(pc))[::-1]
    for yi, r in zip(y, pc.itertuples()):
        ax.plot([r.rho_raw, r.rho_GC_controlled], [yi, yi],
                color="#BBBBBB", lw=1.0, zorder=1)
    ax.scatter(pc.rho_raw, y, s=15, color="#0072B2", zorder=3, label="raw")
    ax.scatter(pc.rho_GC_controlled, y, s=15, facecolor="white",
               edgecolor="#0072B2", linewidth=0.9, zorder=3,
               label="GC-controlled")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels([LBL[c] for c in pc.cls], fontsize=6.3)
    ax.set_xlabel("Spearman ρ (class vs phastCons)")
    ax.legend(loc="lower right", fontsize=5.8)
    panel(ax, "a")
    ax2 = axes[0, 1]
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
    ax3 = axes[0, 2]
    y = np.arange(len(pp))[::-1]
    for yi, r in zip(y, pp.itertuples()):
        ax3.plot([r.rho_raw, r.rho_GC_controlled], [yi, yi],
                 color="#BBBBBB", lw=1.0, zorder=1)
    ax3.scatter(pp.rho_raw, y, s=15, color="#0072B2", zorder=3, label="raw")
    ax3.scatter(pp.rho_GC_controlled, y, s=15, facecolor="white",
                edgecolor="#0072B2", linewidth=0.9, zorder=3,
                label="GC-controlled")
    ax3.axvline(0, color="k", lw=0.6)
    ax3.set_yticks(y); ax3.set_yticklabels([LBL[c] for c in pp.cls],
                                           fontsize=6.3)
    ax3.set_xlabel("Spearman ρ (class vs phyloP)")
    ax3.legend(loc="lower right", fontsize=5.8)
    panel(ax3, "c")
    ax4 = axes[1, 0]
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
    # class-vs-PCG tests for the all-promoter bars
    ppub = load("phylop/pub_stats_phylop.tsv")
    for xi, c in zip(x, order4):
        if c == "PCG":
            continue
        rr = ppub[(ppub.feature == "frac accelerated (phylop<=-2)")
                  & (ppub.cls == c)]
        if len(rr):
            yv = tp.loc[c, "frac_accel_2"]
            ly = tp.loc[c, "frac_accel2_ltr"]
            top = max(yv, ly) if not np.isnan(ly) else yv
            ax4.text(xi - w / 2, top, stars(rr.fdr.iloc[0]),
                     ha="center", va="bottom", fontsize=5.5)
    ax4.margins(y=0.16)
    ax4.set_xticks(x)
    ax4.set_xticklabels([LBL[c].replace("Protein-coding", "PCG")
                         .replace("Transcribed pseudo", "tPseudo") for c in order4],
                        rotation=42, ha="right")
    ax4.set_ylabel("fraction of bases phyloP ≤ −2")
    # The legend covered the intronic star inside the axes. Put it above,
    # where the bars cannot reach, and give the stars headroom.
    ax4.set_ylim(0, 0.078)
    ax4.legend(fontsize=5.5, loc="lower center", bbox_to_anchor=(0.5, 1.005),
               ncol=2, handlelength=0.9, columnspacing=0.9,
               handletextpad=0.3)
    panel(ax4, "d")
    # -- e: transcript-age composition by class --------------------------
    axe = axes[1, 1]
    age_order = ["human_specific", "primate", "mammalian", "eutherian"]
    age_lab = ["Human-specific", "Primate", "Mammalian", "Eutherian"]
    age_col = ["#D9D9D9", "#56B4E9", "#0072B2", "#003B5C"]
    grp = [("PCG", ev[ev.group == "pcg"])] + \
          [(c, ev[ev.cls == c]) for c in LNC5]
    bottom = np.zeros(len(grp))
    for a, lab, col in zip(age_order, age_lab, age_col):
        vals = np.array([(g.age_class == a).mean() for _, g in grp])
        axe.bar(np.arange(len(grp)), vals, bottom=bottom, color=col,
                label=lab, width=0.76)
        bottom += vals
    axe.set_xticks(np.arange(len(grp)))
    axe.set_xticklabels(["PCG"] + [LBL[c] for c in LNC5],
                        rotation=42, ha="right", fontsize=6)
    axe.set_ylabel("fraction of loci")
    axe.set_ylim(0, 1)
    # Legend above the axes. Inside the panel it covered the stacked bars.
    axe.legend(fontsize=5, loc="upper center", bbox_to_anchor=(0.5, 1.16),
               ncol=4, handlelength=1.0, columnspacing=0.8,
               handletextpad=0.3)
    # The letter for this panel is added by the figure builder, so the panel
    # crop used in the chapter composites carries none.
    # -- f: human-macaque non-CpG substitution rate by class --------------
    axf = axes[1, 2]
    boxes = [sub[sub.group == "pcg"].rate_non.dropna().values] + \
            [sub[sub.cls == c].rate_non.dropna().values for c in LNC5]
    labs = ["PCG"] + [LBL[c] for c in LNC5]
    cols = [C["PCG"]] + [C[c] for c in LNC5]
    boxes.append(sub[(sub.cls == "intergenic") &
                     (sub.tss_in_ltr == 1)].rate_non.dropna().values)
    labs.append("Intergenic\nLTR")
    cols.append("#D55E00")
    bp = axf.boxplot(boxes, patch_artist=True, showfliers=False,
                     widths=0.6, medianprops=dict(color="k", lw=0.8))
    for patch, col in zip(bp["boxes"], cols):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
        patch.set_linewidth(0.5)
    axf.set_xticks(range(1, len(boxes) + 1))
    axf.set_xticklabels(labs, rotation=42, ha="right", fontsize=6)
    axf.set_ylabel("substitutions/site (non-CpG)")
    panel(axf, "f")
    fig.tight_layout(w_pad=1.2, h_pad=1.5)
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
    # Drop the parenthetical and the dash so GTEx names fit a single column.
    lab = [s.split(" (")[0].replace(" - ", " ") for s in t.tissue]
    ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=5.6)
    ax.set_xlabel("share of expressed lncRNAs (GTEx)")
    ax.legend(fontsize=5, loc="lower center", bbox_to_anchor=(0.5, 1.02),
              ncol=3, handlelength=1.0, columnspacing=0.8,
              handletextpad=0.3)
    # annotate testis LTR
    ti = [i for i, s in enumerate(t.tissue) if "Testis" in s]
    if ti:
        ax.text(1.02, y[ti[0]], "LTR-rich", ha="left", va="center",
                fontsize=6, color="#D55E00", fontweight="bold",
                clip_on=False)
    # Single-panel figure: the legend carries no panel letter, so neither
    # should the image.
    save(fig, "Fig7_gtex")


# ============================================================ Fig 8 ==========
def fig8():
    de = load("tissue_disease/tables/tcga_class_de_per_cancer.tsv")
    lt = load("tissue_disease/tables/tcga_ltr_derepression.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(177 * MM, 78 * MM),
                             gridspec_kw={"width_ratios": [1.15, 1, 1],
                                          "wspace": 0.45})
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
    plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02, label="Fraction up")
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
    # -- c: perturbation validation (DUX4 / DNMTi), LTR-lnc response shift ---
    pt = load("perturb/perturb_tests.tsv")
    pt = pt[pt.contrast == "LTR_lnc_vs_other"]
    order = ["GSE85461_iDUX4_fastq", "GSE85461_lentiDUX4_14h",
             "GSE56787_FSHD_vs_ctrl", "GSE337536_decitabine",
             "GSE305426_DNMTi_vs_DMSO"]
    # Single-line labels. Two-line labels collided once rotated; the accession
    # numbers are in the legend and in Methods, so they are not repeated here.
    labels = ["iDUX4", "lenti-DUX4", "FSHD", "decitabine", "DAC/AZA"]
    pt = pt.set_index("dataset").reindex(order)
    ax3 = axes[2]
    x = np.arange(len(order))
    cols = ["#D55E00", "#E69F00", "#999999", "#0072B2", "#56B4E9"]
    ax3.bar(x, pt.cliffs.values, width=0.62, color=cols)
    ax3.axhline(0, color="k", lw=0.7)
    for xi, (d, pv) in enumerate(zip(pt.cliffs.values, pt.pval.values)):
        s = "***" if pv < 1e-10 else ("**" if pv < 1e-3 else
                                      ("*" if pv < 0.05 else "ns"))
        ax3.text(xi, d + 0.012, s, ha="center", fontsize=7)
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels, fontsize=5.6, rotation=30, ha="right",
                        rotation_mode="anchor")
    ax3.set_ylabel("Cliff's δ, LTR-derived vs other lncRNA\nresponse to perturbation")
    panel(ax3, "c")
    save(fig, "Fig8_tcga")


# ============================================================ Fig 9 ==========
def fig9():
    # Panel a was a five-point class-mean scatter (GTEx activity vs TCGA
    # up-regulation, rho = 0.70, n = 5) whose message flipped between 0.4 and
    # 1.0 under leave-one-class-out; it was retired at round-7 review and
    # replaced with the per-pair matched-tau residual distributions, which
    # carry the expression-controlled specificity contrast the panel is for.
    iv = None
    fg = load("final_gaps/fantom_vs_gtex_tau.tsv")
    # wspace is set explicitly: panel b is text only, and panel c's y labels run
    # left far enough to collide with it at the default spacing.
    fig, axes = plt.subplots(1, 3, figsize=(177 * MM, 58 * MM),
                             gridspec_kw={"width_ratios": [1.0, 0.72, 1.2],
                                          "wspace": 0.42})
    ax = axes[0]
    # per-pair tau residual (lncRNA minus its expression-matched PCG), one
    # striped box per class; class medians annotated at right
    tm = load("tissue_disease/tables/tau_matched_pairs.tsv")
    tm["diff"] = tm["tau_lnc"] - tm["tau_pcg"]
    order = [c for c in ["divergent", "sense_overlapping", "antisense",
                         "intronic", "intergenic"] if c in set(tm["cls"])]
    data = [tm.loc[tm["cls"] == c, "diff"].to_numpy() for c in order]
    bp = ax.boxplot(data, vert=False, showfliers=False, widths=0.6,
                    patch_artist=True, whis=(5, 95),
                    medianprops=dict(color="black", linewidth=1))
    for patch, c in zip(bp["boxes"], order):
        patch.set_facecolor(C[c])
        patch.set_alpha(0.55)
    ax.axvline(0, color="#999", linewidth=0.8, linestyle=":", zorder=1)
    ax.set_yticklabels([LBL[c] for c in order], fontsize=6)
    ax.set_xlabel("GTEx tau residual vs matched PCG")
    for i, c in enumerate(order):
        med = float(np.median(data[i]))
        ax.annotate(f"{med:+.3f}", (1.02, 1 - (i + 0.5) / len(order)),
                    xycoords="axes fraction", fontsize=5.5, ha="left",
                    va="center", color="#333", annotation_clip=False)
    panel(ax, "a")
    ax2 = axes[1]
    ax2.axis("off")
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.text(0.02, 0.92, "Concordance and\nmatched specificity", fontsize=7.2,
             fontweight="bold", va="top")
    ax2.text(0.08, 0.70,
             "CAGE promoter vs\nRNA-seq body\ntissue specificity\n(gene level)",
             fontsize=6.5, va="top")
    ax2.text(0.08, 0.44,
             "ρ = {:.2f}\nn = {:,} genes".format(fg.rho.iloc[0],
                                                 int(fg.n.iloc[0])),
             fontsize=7, va="top")
    # The raw class tau is an expression artefact, since tau rises to 1 as
    # detection falls. Only the difference against a coding gene of matched
    # expression is interpretable, so that is what the panel shows.
    # The clipped-run value lives in the legend: printing it here pushed the
    # fourth line onto the per-class values below, which are the panel's point.
    ax2.text(0.08, 0.26, "GTEx τ excess over\nmatched PCG", fontsize=6.5,
             va="top")
    ax2.text(0.08, 0.14,
             "divergent +0.03\nantisense +0.03\nintronic +0.03\n"
             "intergenic +0.04\nsense-overlap +0.00",
             fontsize=6, va="top", color="#333")
    panel(ax2, "b")
    # -- c: class detection rank in human and mouse single cells -----------
    # Ranks, not absolute detection, because the two annotations differ in
    # completeness. Each row shows how often a class occupies each rank.
    det = ["PCG_det", "sense_overlapping_det", "divergent_det",
           "antisense_det", "intergenic_det", "intronic_det"]
    lab = {"PCG_det": "PCG", "sense_overlapping_det": "Sense-overlap",
           "divergent_det": "Divergent", "antisense_det": "Antisense",
           "intergenic_det": "Intergenic", "intronic_det": "Intronic"}
    hu = load("scbasecount/sample_class_detection.tsv")[det]
    # the 604-sample mouse set, not the 120-sample pilot: at this size every
    # class carries enough draws for the rank distribution to be stable
    mo = load("scbasecount/mouse_sample_detection_all.tsv")[det]
    ax3 = axes[2]
    shades = ["#08306B", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1",
              "#DEEBF7"]
    for i, d in enumerate(det):
        for off, df in ((0.20, hu), (-0.20, mo)):
            r = df.rank(axis=1, ascending=False)[d].value_counts(
                normalize=True).reindex(range(1, 7), fill_value=0)
            left = 0.0
            for rk, frac in r.items():
                ax3.barh(i + off, frac, left=left, height=0.34,
                         color=shades[rk - 1], edgecolor="white", lw=0.2)
                left += frac
    ax3.set_yticks(range(len(det)))
    ax3.set_yticklabels([lab[d] for d in det], fontsize=6.5)
    ax3.set_xlim(0, 1)
    ax3.set_xlabel("Fraction of samples by rank\n(top row human, bottom mouse)",
                   fontsize=6.5)
    ax3.set_title("Single-cell class hierarchy", fontsize=7)
    # One shared rank key under the whole figure. Sitting it under panel c only
    # left an empty band beneath panels a and b and forced the cropped height up.
    handles = [plt.Line2D([], [], color=shades[rk - 1], lw=4)
               for rk in range(1, 7)]
    fig.legend(handles, [str(rk) for rk in range(1, 7)], title="Rank",
               loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.14),
               fontsize=6, title_fontsize=6, handlelength=1.4,
               columnspacing=1.3, handletextpad=0.5)
    panel(ax3, "c")
    save(fig, "Fig9_integration")


def main():
    with step("16_submission_figures") as st:
        fig1(); fig2(); fig3(); fig3b(); fig4(); fig4b(); fig5(); fig6()
        fig7(); fig8(); fig9()
        st.note(f"submission figures -> {FIGS}")


if __name__ == "__main__":
    main()
