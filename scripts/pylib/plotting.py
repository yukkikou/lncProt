#!/usr/bin/env python
"""plotting.py — shared figure style: class order/colors, save helpers.
All figure text in English. Publication-friendly defaults.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKDIR = "/media/share/workdir/1610305236/lncProt"
FIGS = os.path.join(WORKDIR, "figures")

# canonical class order + colors (colorblind-safe)
CLASS_ORDER = ["PCG", "divergent", "antisense", "intronic",
               "sense_overlapping", "intergenic", "pseudogene",
               "transcribed_pseudogene"]
CLASS_LABELS = {
    "PCG": "Protein-coding", "divergent": "Divergent lncRNA",
    "antisense": "Antisense lncRNA", "intronic": "Intronic lncRNA",
    "sense_overlapping": "Sense-overlap lncRNA", "intergenic": "Intergenic lncRNA",
    "pseudogene": "Pseudogene", "transcribed_pseudogene": "Transcribed pseudogene"}
CLASS_COLORS = {
    "PCG": "#4C72B0", "divergent": "#DD8452", "antisense": "#55A868",
    "intronic": "#C44E52", "sense_overlapping": "#8172B3",
    "intergenic": "#937860", "pseudogene": "#DA8BC3",
    "transcribed_pseudogene": "#8C8C8C", "PCG_active": "#4C72B0"}

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "font.size": 9,
    "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"saved figures/{name}.png|.pdf")
