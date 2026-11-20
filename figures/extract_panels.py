#!/usr/bin/env python
"""Split Fig. 5 (conservation) and Fig. 7 (tumour) into single panels.

Both figures are needed by more than one chapter, so the panels are rebuilt
from the plotting code with the panel letters suppressed. Each axes is saved
with its own tight bounding box, which keeps tick labels, titles and any
legend anchored to that axes.
"""
import importlib.util
import os
import sys

import matplotlib
matplotlib.use("Agg")

ROOT = "/media/share/workdir/1610305236/lncProt"
sys.path.insert(0, os.path.join(ROOT, "scripts"))
OUT = os.path.join(ROOT, "figure_build", "work", "panels")
os.makedirs(OUT, exist_ok=True)

spec = importlib.util.spec_from_file_location(
    "s16", os.path.join(ROOT, "scripts", "16_submission_figures.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def extract(fn, name, letters):
    captured = {}

    def fake_save(fig, _name):
        captured["fig"] = fig

    m.panel = lambda *a, **k: None
    m.save = fake_save
    fn()
    fig = captured["fig"]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes = fig.axes
    print(name, "axes", len(axes), "expected", len(letters))
    for ax, letter in zip(axes, letters):
        bbox = ax.get_tightbbox(renderer).transformed(
            fig.dpi_scale_trans.inverted())
        base = os.path.join(OUT, f"{name}_{letter}")
        fig.savefig(base + ".pdf", bbox_inches=bbox)
        fig.savefig(base + ".png", dpi=600, bbox_inches=bbox)
        print("  saved", base, "bbox", tuple(round(v, 2) for v in bbox.extents))


if __name__ == "__main__":
    extract(m.fig6, "Fig5", list("abcdef"))
    extract(m.fig8, "Fig7", list("abc"))
