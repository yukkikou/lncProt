#!/usr/bin/env python
"""Fig. 1B: coordinate classes across six mammals.

Bars are the five coordinate classes per species, from module 59 (macaque,
rat, rabbit, opossum plus its mouse rerun) and module 37 (mouse, the version
the manuscript cites). Overlap-defined classes are structurally absent in the
four non-mouse species because their annotations carry no lncRNA-coding exon
overlap.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager as _fm

ROOT = "/media/share/workdir/1610305236/lncProt"
OUT = os.path.join(ROOT, "figure_build", "work", "panels")
os.makedirs(OUT, exist_ok=True)

for p in ("/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
          "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf"):
    if os.path.exists(p):
        _fm.fontManager.addfont(p)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "pdf.fonttype": 42, "font.size": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6,
})
MM = 1 / 25.4

C = {"divergent": "#E69F00", "antisense": "#009E73", "intronic": "#5D3A9B",
     "sense_overlapping": "#CC79A7", "intergenic": "#56B4E9"}
LBL = {"divergent": "Divergent", "antisense": "Antisense", "intronic": "Intronic",
       "sense_overlapping": "Sense-overlap", "intergenic": "Intergenic"}
CLS = ["divergent", "antisense", "intronic", "sense_overlapping", "intergenic"]

# module 37 mouse is the version cited in the text
tax = pd.read_csv(os.path.join(ROOT, "results", "cross_species",
                               "taxonomy_summary.tsv"), sep="\t")
tax = tax[tax.species != "mouse"]
rows = [("human", 19258, dict(zip(
            ["intergenic", "intronic", "antisense", "divergent",
             "sense_overlapping"], [10035, 3403, 2778, 2692, 350]))),
        ("mouse", 9959, dict(zip(
            ["divergent", "antisense", "intronic", "sense_overlapping",
             "intergenic"], [1419, 1342, 1318, 92, 5788])))]
for r in tax.itertuples():
    rows.append((r.species, int(r.lnc),
                 {"divergent": int(r.divergent), "antisense": int(r.antisense),
                  "intronic": int(r.intronic),
                  "sense_overlapping": int(r.sense_overlapping),
                  "intergenic": int(r.intergenic)}))

fig, ax = plt.subplots(figsize=(90 * MM, 62 * MM))
n = len(rows)
x = np.arange(n)
bottom = np.zeros(n)
order = ["intergenic", "intronic", "sense_overlapping", "antisense",
         "divergent"]
for cls in order:
    vals = np.array([r[2].get(cls, 0) / r[1] * 100 for r in rows])
    ax.bar(x, vals, bottom=bottom, width=0.62, color=C[cls],
           edgecolor="white", linewidth=0.3, label=LBL[cls])
    for xi, (v, b) in enumerate(zip(vals, bottom)):
        if v >= 4:
            ax.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                    fontsize=5.2,
                    color="white" if cls in ("intergenic", "intronic")
                    else "#333")
    bottom += vals
ax.set_xticks(x)
ax.set_xticklabels([r[0].capitalize() for r in rows])
ax.set_ylabel("Class share of annotated lncRNAs (%)")
ax.set_ylim(0, 100)
ax.legend(fontsize=5.4, loc="center left", bbox_to_anchor=(1.01, 0.5),
          handlelength=1.1, labelspacing=0.3, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "Fig1B_species.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "Fig1B_species.png"), dpi=600,
            bbox_inches="tight")
print("saved Fig1B_species")
