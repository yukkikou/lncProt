#!/usr/bin/env python
"""Fig. 1 v8: one style, seven panels, no duplicated trait.

A schematic, B counts, C six-species shares, D five-trait gradient,
E GC, F CpG o/e, G TATA. The old high-CpG panel is folded into D.
"""
import importlib.util
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as _fm
from matplotlib.patches import Rectangle

ROOT = os.environ.get("LNCPROT_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "figure_build", "work", "panels")
os.makedirs(OUT, exist_ok=True)
MM = 1 / 25.4

for p in ("/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
          "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf"):
    _fm.fontManager.addfont(p)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "pdf.fonttype": 42, "font.size": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6.2, "ytick.labelsize": 6.2, "legend.fontsize": 6,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
})
C = {"PCG": "#4A6FA5", "divergent": "#E69F00", "antisense": "#009E73",
     "intronic": "#5D3A9B", "sense_overlapping": "#CC79A7",
     "intergenic": "#56B4E9"}
LBL = {"PCG": "PCG", "divergent": "Divergent", "antisense": "Antisense",
       "intronic": "Intronic", "sense_overlapping": "Sense-overlap",
       "intergenic": "Intergenic"}
LNC = ["divergent", "antisense", "intronic", "sense_overlapping", "intergenic"]
ORDER = ["PCG"] + LNC


def panel(ax, s):
    ax.text(-0.02, 1.06, s, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="right")


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=600,
                bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def schematic():
    """Five-class schematic, redrawn with the same style as the rest."""
    fig, ax = plt.subplots(figsize=(90 * MM, 58 * MM))
    ax.set_xlim(0, 10); ax.set_ylim(-0.6, 5.2); ax.axis("off")
    rows = [("Divergent", "divergent", "div"), ("Antisense", "antisense", "anti"),
            ("Intronic", "intronic", "intro"), ("Sense-overlap", "sense_overlapping", "sense"),
            ("Intergenic", "intergenic", "inter")]
    for i, (name, cls, kind) in enumerate(rows):
        y = 4.4 - i * 0.95
        ax.text(0.05, y + 0.28, name, fontsize=6.8, fontweight="bold",
                color=C[cls], va="center")
        # coding gene
        ax.add_patch(Rectangle((2.2, y + 0.05), 3.2, 0.22, facecolor=C["PCG"],
                               edgecolor="none"))
        ax.plot([1.6, 6.6], [y + 0.16, y + 0.16], color="#cccccc", lw=0.8,
                zorder=0)
        lnc_y = y - 0.28
        if kind == "div":
            ax.add_patch(Rectangle((2.0, lnc_y), 3.0, 0.2, facecolor=C[cls],
                                   edgecolor="none"))
            ax.plot([2.0], [lnc_y + 0.1], marker="|", color="k", ms=5)
            ax.annotate("", xy=(1.7, lnc_y + 0.1), xytext=(2.0, lnc_y + 0.1),
                        arrowprops=dict(arrowstyle="->", lw=0.7, color="k"))
        elif kind == "anti":
            ax.add_patch(Rectangle((2.6, lnc_y), 2.6, 0.2, facecolor=C[cls],
                                   edgecolor="none"))
            ax.annotate("", xy=(2.3, lnc_y + 0.1), xytext=(2.6, lnc_y + 0.1),
                        arrowprops=dict(arrowstyle="->", lw=0.7, color="k"))
            ax.plot([5.2], [lnc_y + 0.1], marker="|", color="k", ms=5)
        elif kind == "intro":
            ax.add_patch(Rectangle((2.8, lnc_y), 1.2, 0.2, facecolor=C[cls],
                                   edgecolor="none"))
            ax.plot([2.8], [lnc_y + 0.1], marker="|", color="k", ms=5)
        elif kind == "sense":
            ax.add_patch(Rectangle((2.5, lnc_y), 2.4, 0.2, facecolor=C[cls],
                                   edgecolor="none"))
            ax.plot([2.5], [lnc_y + 0.1], marker="|", color="k", ms=5)
        else:
            ax.add_patch(Rectangle((2.0, lnc_y), 0.9, 0.2, facecolor=C[cls],
                                   edgecolor="none"))
            ax.add_patch(Rectangle((5.6, lnc_y), 0.9, 0.2, facecolor=C[cls],
                                   edgecolor="none"))
            ax.plot([2.0], [lnc_y + 0.1], marker="|", color="k", ms=5)
    ax.plot([2.2], [4.9], marker="s", color=C["PCG"], ms=3)
    ax.text(2.35, 4.9, "coding gene", fontsize=6, va="center")
    ax.plot([3.6], [4.9], marker="s", color=C["divergent"], ms=3)
    ax.text(3.75, 4.9, "lncRNA", fontsize=6, va="center")
    ax.text(4.9, 4.9, "| TSS", fontsize=6, va="center")
    panel(ax, "A")
    save(fig, "F1A_schematic")


def counts():
    gm = pd.read_csv(os.path.join(ROOT, "results", "gene_master_table.tsv"),
                     sep="\t", low_memory=False)
    n = gm[gm.group == "lnc"].cls.value_counts().reindex(LNC)
    fig, ax = plt.subplots(figsize=(52 * MM, 42 * MM))
    ax.barh(range(len(LNC)), n.values, color=[C[c] for c in LNC],
            height=0.62, edgecolor="white", linewidth=0.4)
    for i, v in enumerate(n.values):
        ax.text(v + 120, i, f"{v:,}", va="center", fontsize=6.3)
    ax.set_yticks(range(len(LNC)))
    ax.set_yticklabels([LBL[c] for c in LNC])
    ax.invert_yaxis()
    ax.set_xlabel("lncRNA genes")
    ax.set_xlim(0, 11500)
    panel(ax, "B")
    save(fig, "F1B_counts")


def species():
    tax = pd.read_csv(os.path.join(ROOT, "results", "cross_species",
                                   "taxonomy_summary.tsv"), sep="\t")
    tax = tax[tax.species != "mouse"]
    rows = [("Human", 19258, {"intergenic": 10035, "intronic": 3403,
                              "antisense": 2778, "divergent": 2692,
                              "sense_overlapping": 350}),
            ("Mouse", 9959, {"divergent": 1419, "antisense": 1342,
                             "intronic": 1318, "sense_overlapping": 92,
                             "intergenic": 5788})]
    for r in tax.itertuples():
        rows.append((r.species.capitalize(), int(r.lnc),
                     {"divergent": int(r.divergent), "antisense": int(r.antisense),
                      "intronic": int(r.intronic),
                      "sense_overlapping": int(r.sense_overlapping),
                      "intergenic": int(r.intergenic)}))
    fig, ax = plt.subplots(figsize=(84 * MM, 42 * MM))
    x = np.arange(len(rows)); bottom = np.zeros(len(rows))
    for cls in ["intergenic", "intronic", "sense_overlapping", "antisense",
                "divergent"]:
        vals = np.array([r[2].get(cls, 0) / r[1] * 100 for r in rows])
        ax.bar(x, vals, bottom=bottom, width=0.62, color=C[cls],
               edgecolor="white", linewidth=0.4, label=LBL[cls])
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([r[0] for r in rows], rotation=25, ha="right")
    ax.set_ylabel("Class share (%)")
    ax.set_ylim(0, 100)
    ax.legend(fontsize=5.6, loc="center left", bbox_to_anchor=(1.01, 0.5),
              frameon=False, handlelength=1.0, labelspacing=0.25)
    panel(ax, "C")
    save(fig, "F1C_species")


def gradient():
    seq = pd.read_csv(os.path.join(ROOT, "results", "seq_features_by_class.tsv"),
                      sep="\t", index_col=0)
    ch = pd.read_csv(os.path.join(ROOT, "results", "chromhmm_by_class.tsv"),
                     sep="\t", index_col=0)
    pl = ch[[c for c in ch.columns if c.endswith("_prom_like")]].mean(axis=1)
    te = pd.read_csv(os.path.join(ROOT, "results", "te_promoter_stats.tsv"),
                     sep="\t", index_col=0, low_memory=False)
    cov = te.groupby("cls").te_cov.mean()
    cons = pd.read_csv(os.path.join(ROOT, "results", "conservation.tsv"),
                       sep="\t", low_memory=True).groupby("cls").phast_mean.mean()
    act = pd.read_csv(os.path.join(ROOT, "results", "mouse_taxonomy",
                                   "cross_species_activity_comparison.tsv"),
                      sep="\t")
    act = act[act.species.isin(["human", ""])].set_index("cls").active_frac
    order = ["PCG"] + LNC
    lab = ["PCG", "Div", "Anti", "Intro", "Sense", "Inter"]
    traits = [("High-CpG (%)", [seq.loc[c, "HCP"] * 100 for c in order]),
              ("Promoter chromatin (%)", [pl[c] * 100 for c in order]),
              ("CAGE-active (%)", [act[c] * 100 for c in order]),
              ("phastCons", [cons[c] for c in order]),
              ("TE coverage (%)", [cov[c] * 100 for c in order])]
    fig, axes = plt.subplots(1, 5, figsize=(170 * MM, 40 * MM),
                             gridspec_kw={"wspace": 0.5})
    for ax, (title, vals) in zip(axes, traits):
        ax.plot(range(6), vals, "-o", color="#333333", ms=2.6, lw=1.0,
                mfc="white", mew=1.0)
        for i, v in enumerate(vals):
            ax.plot([i], [v], "o", ms=2.6, color=C[order[i]], zorder=3)
        ax.set_xticks(range(6)); ax.set_xticklabels(lab, rotation=45, ha="right",
                                                    fontsize=5.6)
        ax.set_title(title, fontsize=6.3, pad=3)
        lo, hi = min(vals), max(vals); pad = (hi - lo) * 0.25 + 1e-6
        ax.set_ylim(lo - pad, hi + pad)
        ax.tick_params(labelsize=5.6)
    panel(axes[0], "D")
    save(fig, "F1D_gradient")


def violins():
    sf = pd.read_csv(os.path.join(ROOT, "results", "seq_features.tsv"),
                     sep="\t", low_memory=False)
    for col, name, ylab in [("gc", "F1E_gc", "GC content"),
                            ("cpg_oe", "F1F_cpgoe", "CpG observed/expected")]:
        data = [sf.loc[sf.cls == c, col].dropna().values for c in ORDER]
        fig, ax = plt.subplots(figsize=(58 * MM, 42 * MM))
        parts = ax.violinplot(data, positions=range(len(ORDER)), widths=0.75,
                              showextrema=False, showmedians=False)
        for b, c in zip(parts["bodies"], ORDER):
            b.set_facecolor(C[c]); b.set_alpha(0.55); b.set_edgecolor("none")
        ax.scatter(range(len(ORDER)), [np.median(d) for d in data], s=6,
                   color="white", edgecolor="k", linewidth=0.4, zorder=3)
        ax.set_xticks(range(len(ORDER)))
        ax.set_xticklabels(["PCG", "Div", "Anti", "Intro", "Sense", "Inter"],
                           rotation=45, ha="right", fontsize=5.8)
        ax.set_ylabel(ylab, fontsize=6.8)
        panel(ax, name.split("_")[0][-1])
        save(fig, name)


def tata():
    sf = pd.read_csv(os.path.join(ROOT, "results", "seq_features.tsv"),
                     sep="\t", low_memory=False)
    fig, ax = plt.subplots(figsize=(58 * MM, 42 * MM))
    fr, lo, hi = [], [], []
    for c in ORDER:
        d = sf.loc[sf.cls == c, "TATA"].dropna()
        p = d.mean(); n = len(d); z = 1.96
        den = 1 + z * z / n
        centre = (p + z * z / (2 * n)) / den
        half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
        fr.append(p * 100); lo.append(max(0, (centre - half) * 100))
        hi.append((centre + half) * 100)
    ax.bar(range(len(ORDER)), fr, color=[C[c] for c in ORDER], width=0.62,
           edgecolor="white", linewidth=0.4)
    ax.errorbar(range(len(ORDER)), fr,
                yerr=[np.array(fr) - np.array(lo), np.array(hi) - np.array(fr)],
                fmt="none", ecolor="#333", lw=0.7, capsize=1.5)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(["PCG", "Div", "Anti", "Intro", "Sense", "Inter"],
                       rotation=45, ha="right", fontsize=5.8)
    ax.set_ylabel("TATA-box promoters (%)", fontsize=6.8)
    panel(ax, "G")
    save(fig, "F1G_tata")


if __name__ == "__main__":
    schematic(); counts(); species(); gradient(); violins(); tata()
