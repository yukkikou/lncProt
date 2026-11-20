#!/usr/bin/env python
"""Redraw two panels that do not survive re-layout.

Fig7c (DUX4/DNMTi response) is naturally portrait inside its 1x3 source
figure; used alone in Chapter 3 it becomes a tall narrow strip. It is redrawn
at 85 x 55 mm.

Fig9b was a text-only box (rho, n, class medians), not a plot, and duplicated
Fig9a numerically. Fig9 is redrawn with panels a and c only, keeping the shared
rank legend.
"""
import importlib.util
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = "/media/share/workdir/1610305236/lncProt"
OUT = os.path.join(ROOT, "figure_build", "work", "panels")
os.makedirs(OUT, exist_ok=True)

spec = importlib.util.spec_from_file_location(
    "s16", os.path.join(ROOT, "scripts", "16_submission_figures.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
MM = m.MM


def fig7c():
    pt = m.load("perturb/perturb_tests.tsv")
    pt = pt[pt.contrast == "LTR_lnc_vs_other"]
    order = ["GSE85461_iDUX4_fastq", "GSE85461_lentiDUX4_14h",
             "GSE56787_FSHD_vs_ctrl", "GSE337536_decitabine",
             "GSE305426_DNMTi_vs_DMSO"]
    labels = ["iDUX4", "lenti-DUX4", "FSHD", "decitabine", "DAC/AZA"]
    pt = pt.set_index("dataset").reindex(order)
    fig, ax3 = plt.subplots(figsize=(85 * MM, 55 * MM))
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
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "Fig7_c_wide.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "Fig7_c_wide.png"), dpi=600,
                bbox_inches="tight")
    plt.close(fig)
    print("saved Fig7_c_wide")


def fig4d_wide():
    d = m.load("divergent_pairs_activity.tsv")
    fig, ax = plt.subplots(figsize=(85 * MM, 55 * MM))
    data = [d["rho_divergent"].dropna(), d["rho_control_same_strand"].dropna()]
    vp = ax.violinplot(data, showmedians=True, widths=0.8)
    for b in vp["bodies"]:
        b.set_facecolor("#0072B2")
        b.set_alpha(0.6)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["divergent\npairs", "same-strand\ncontrols"])
    ax.set_ylabel("Spearman ρ of CAGE activity")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "Fig4d_wide.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "Fig4d_wide.png"), dpi=600,
                bbox_inches="tight")
    plt.close(fig)
    print("saved Fig4d_wide")


def fig5a_wide():
    pc = m.load("robustness_partial_corr.tsv")
    fig, ax = plt.subplots(figsize=(90 * MM, 55 * MM))
    y = np.arange(len(pc))[::-1]
    for yi, r in zip(y, pc.itertuples()):
        ax.plot([r.rho_raw, r.rho_GC_controlled], [yi, yi],
                color="#BBBBBB", lw=1.0, zorder=1)
    ax.scatter(pc.rho_raw, y, s=15, color="#0072B2", zorder=3, label="raw")
    ax.scatter(pc.rho_GC_controlled, y, s=15, facecolor="white",
               edgecolor="#0072B2", linewidth=0.9, zorder=3,
               label="GC-controlled")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([m.LBL[c] for c in pc.cls], fontsize=6.3)
    ax.set_xlabel("Spearman ρ (class vs phastCons)")
    ax.legend(loc="lower right", fontsize=5.8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "Fig5a_wide.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "Fig5a_wide.png"), dpi=600,
                bbox_inches="tight")
    plt.close(fig)
    print("saved Fig5a_wide")


def fig9abc():
    tm = m.load("tissue_disease/tables/tau_matched_pairs.tsv")
    tm["diff"] = tm["tau_lnc"] - tm["tau_pcg"]
    order = [c for c in ["divergent", "sense_overlapping", "antisense",
                         "intronic", "intergenic"] if c in set(tm["cls"])]
    data = [tm.loc[tm["cls"] == c, "diff"].to_numpy() for c in order]
    fig, axes = plt.subplots(1, 3, figsize=(177 * MM, 58 * MM),
                             gridspec_kw={"width_ratios": [1.0, 0.9, 1.2],
                                          "wspace": 0.46})
    ax = axes[0]
    bp = ax.boxplot(data, vert=False, showfliers=False, widths=0.6,
                    patch_artist=True, whis=(5, 95),
                    medianprops=dict(color="black", linewidth=1))
    for patch, c in zip(bp["boxes"], order):
        patch.set_facecolor(m.C[c])
        patch.set_alpha(0.55)
    ax.axvline(0, color="#999", linewidth=0.8, linestyle=":", zorder=1)
    ax.set_yticklabels([m.LBL[c] for c in order], fontsize=6)
    ax.set_xlabel("GTEx tau residual vs matched PCG")
    for i, c in enumerate(order):
        med = float(np.median(data[i]))
        ax.annotate(f"{med:+.3f}", (1.02, 1 - (i + 0.5) / len(order)),
                    xycoords="axes fraction", fontsize=5.5, ha="left",
                    va="center", color="#333", annotation_clip=False)
    m.panel(ax, "a")

    # b: per-gene concordance between the two tissue-specificity measures
    tg = m.load("final_gaps/fantom_vs_gtex_tau_genes.tsv")
    tg = tg.dropna(subset=["tau_rna", "tau_cage"])
    classes = [c for c in m.ORDER if c in set(tg.cls)]
    ax2 = axes[1]
    for c in classes:
        d = tg[tg.cls == c]
        col = "#BBBBBB" if c == "PCG" else m.C[c]
        ax2.scatter(d.tau_rna, d.tau_cage, s=0.6, alpha=0.25, color=col,
                    linewidths=0, rasterized=True,
                    label=m.LBL[c] if c != "PCG" else "PCG")
    ax2.plot([0, 1.05], [0, 1.05], ls="--", color="#666", lw=0.7)
    ax2.set_xlim(0, 1.05)
    ax2.set_ylim(0.60, 1.01)
    ax2.set_xlabel("RNA-seq τ (GTEx)", fontsize=6.5)
    ax2.set_ylabel("CAGE τ (FANTOM5)", fontsize=6.5)
    ax2.tick_params(labelsize=5.8)
    ax2.text(0.04, 0.06, "ρ = 0.72\nn = 19,730", transform=ax2.transAxes,
             fontsize=6, va="bottom")
    ax2.legend(fontsize=4.4, loc="lower right", markerscale=8,
               handletextpad=0.2, borderaxespad=0.2, labelspacing=0.2)
    m.panel(ax2, "b")

    det = ["PCG_det", "sense_overlapping_det", "divergent_det",
           "antisense_det", "intergenic_det", "intronic_det"]
    lab = {"PCG_det": "PCG", "sense_overlapping_det": "Sense-overlap",
           "divergent_det": "Divergent", "antisense_det": "Antisense",
           "intergenic_det": "Intergenic", "intronic_det": "Intronic"}
    hu = m.load("scbasecount/sample_class_detection.tsv")[det]
    mo = m.load("scbasecount/mouse_sample_detection_all.tsv")[det]
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
    m.panel(ax3, "c")
    handles = [plt.Line2D([], [], color=shades[rk - 1], lw=4)
               for rk in range(1, 7)]
    fig.legend(handles, [str(rk) for rk in range(1, 7)], title="Rank",
               loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.14),
               fontsize=6, title_fontsize=6, handlelength=1.4,
               columnspacing=1.3, handletextpad=0.5)
    fig.savefig(os.path.join(OUT, "Fig9_abc.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "Fig9_abc.png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("saved Fig9_abc")


if __name__ == "__main__":
    fig4d_wide()
    fig5a_wide()
    fig7c()
    fig9abc()
