#!/usr/bin/env python
"""Six chapter figures, one per chapter, no panel shared across chapters.

Whole figures are relettered in place (lowercase -> uppercase, continuous
across the chapter). Fig. 5 and Fig. 7 are split into single-panel crops
because their panels belong to different chapters; those crops receive their
letter when placed.
"""
import importlib.util
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = "/media/share/workdir/1610305236/lncProt"
SRC = os.path.join(ROOT, "manuscript_package", "figures")
PAN = os.path.join(ROOT, "figure_build", "work", "panels")
OUT = os.path.join(SRC, "chapters")
DPI = 600
MM = DPI / 25.4
FONT = "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf"
W = int(183 * MM)
GAP = int(5 * MM)

os.makedirs(OUT, exist_ok=True)

spec = importlib.util.spec_from_file_location(
    "bc", os.path.join(ROOT, "figure_build", "build_composites.py"))
bc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bc)


def whole(pdf, start, width_mm=183, align_top=False):
    im, mapping = bc.reletter(pdf, start, align_top=align_top)
    im = bc.crop(im)
    target = int(width_mm * MM)
    if im.width > target:
        im = im.resize((target, int(im.height * target / im.width)),
                       Image.LANCZOS)
    return im, mapping


def crop_panel(letter, name, width_mm=None):
    im = Image.open(os.path.join(PAN, f"{name}.png")).convert("RGB")
    band = int(9 * MM)
    canvas = Image.new("RGB", (im.width, im.height + band), "white")
    canvas.paste(im, (0, band))
    d = ImageDraw.Draw(canvas)
    f = ImageFont.truetype(FONT, int(9 * DPI / 72))
    d.text((int(1 * MM), int(6 * MM)), letter, font=f, fill="black",
           anchor="ls")
    if width_mm:
        t = int(width_mm * MM)
        if canvas.width > t:
            canvas = canvas.resize((t, int(canvas.height * t / canvas.width)),
                                   Image.LANCZOS)
    return canvas


def fit_row(imgs, max_up=1.25):
    """Scale one row of panels to fill 183 mm, then centre each vertically."""
    gaps = GAP * (len(imgs) - 1)
    total = sum(i.width for i in imgs)
    f = min((W - gaps) / total, max_up)
    scaled = [i.resize((int(round(i.width * f)), int(round(i.height * f))),
                       Image.LANCZOS) for i in imgs]
    h = max(i.height for i in scaled)
    w = sum(i.width for i in scaled) + gaps
    canvas = Image.new("RGB", (max(w, W), h), "white")
    x = 0
    for i in scaled:
        canvas.paste(i, (x, (h - i.height) // 2))
        x += i.width + GAP
    return canvas


def compose(rows, name, max_up=1.25):
    rows = [fit_row(r, max_up) for r in rows]
    width = max(r.width for r in rows)
    height = sum(r.height for r in rows) + GAP * (len(rows) - 1)
    canvas = Image.new("RGB", (width, height), "white")
    y = 0
    for r in rows:
        canvas.paste(r, ((width - r.width) // 2, y))
        y += r.height + GAP
    limit = int(235 * MM)
    if canvas.height > limit:
        scale = limit / canvas.height
        canvas = canvas.resize((int(canvas.width * scale), limit),
                               Image.LANCZOS)
    canvas.save(os.path.join(OUT, name + ".png"), dpi=(DPI, DPI))
    canvas.save(os.path.join(OUT, name + ".pdf"), resolution=DPI)
    print(f"{name}: {canvas.width}x{canvas.height} px "
          f"({canvas.width / MM:.0f} x {canvas.height / MM:.0f} mm)")
    return canvas


def main():
    maps = {}

    f1, m1 = whole(os.path.join(SRC, "Fig1_classes.pdf"), "A")
    f2, m2 = whole(os.path.join(SRC, "Fig2_sequence.pdf"), "C")
    compose([[f1], [f2]], "Chapter1_architecture")
    maps["Chapter1_architecture"] = [("Fig. 1", o, n) for n, o in m1] + \
                                    [("Fig. 2", o, n) for n, o in m2]

    f3, m3 = whole(os.path.join(SRC, "Fig3_TFs.pdf"), "A")
    f4, m4 = whole(os.path.join(SRC, "Fig4_activity.pdf"), "C")
    f4d, m4d = whole(os.path.join(SRC, "Fig4d_pair_rho.pdf"), "F")
    p5a = crop_panel("G", "Fig5_a", 70)
    compose([[f3, p5a], [f4, f4d]], "Chapter2_coding_adjacent")
    maps["Chapter2_coding_adjacent"] = [("Fig. 3", o, n) for n, o in m3] + \
                                       [("Fig. 4", o, n) for n, o in m4] + \
                                       [("Fig. 4", o, n) for n, o in m4d] + \
                                       [("Fig. 5", "a", "G")]

    f6, m6 = whole(os.path.join(SRC, "Fig6_TE.pdf"), "A")
    f6d, m6d = whole(os.path.join(SRC, "Fig3d_TE_age.pdf"), "D")
    p7c, m7c = whole(os.path.join(PAN, "Fig7_c_wide.pdf"), "E")
    p5b = crop_panel("F", "Fig5_b", 80)
    compose([[f6, p5b], [f6d, p7c]], "Chapter3_te_host_decayed")
    maps["Chapter3_te_host_decayed"] = [("Fig. 6", o, n) for n, o in m6] + \
        [("Fig. 6", o, n) for n, o in m6d] + \
        [("Fig. 7", "c", "E"), ("Fig. 5", "b", "F")]

    p5e = crop_panel("A", "Fig5_e", 85)
    p5c = crop_panel("B", "Fig5_c", 85)
    p5d = crop_panel("C", "Fig5_d", 85)
    p5f = crop_panel("D", "Fig5_f", 85)
    compose([[p5e, p5c], [p5d, p5f]], "Chapter4_age_constraint", max_up=1.8)
    maps["Chapter4_age_constraint"] = [("Fig. 5", "e", "A"), ("Fig. 5", "c", "B"),
                                       ("Fig. 5", "d", "C"), ("Fig. 5", "f", "D")]

    f8, m8 = whole(os.path.join(SRC, "Fig8_gtex.pdf"), "A")
    p7a = crop_panel("B", "Fig7_a", 88)
    p7b = crop_panel("C", "Fig7_b", 88)
    f9, m9 = whole(os.path.join(PAN, "Fig9_abc.pdf"), "D")
    compose([[f8, p7a, p7b], [f9]], "Chapter5_specificity")
    # redrawn Fig9: a -> D, the new scatter b -> E, detection c -> F
    maps["Chapter5_specificity"] = [("Fig. 8", "-", "A"),
                                    ("Fig. 7", "a", "B"), ("Fig. 7", "b", "C"),
                                    ("Fig. 9", "a", "D"), ("Fig. 9", "b", "E"),
                                    ("Fig. 9", "c", "F")]

    f10, m10 = whole(os.path.join(SRC, "Fig10_ml.pdf"), "A", align_top=True)
    compose([[f10]], "Chapter6_prediction")
    maps["Chapter6_prediction"] = [("Fig. 10", o, n) for n, o in m10]

    with open(os.path.join(OUT, "mapping.tsv"), "w") as fh:
        for chap, rows in maps.items():
            for src, old, new in rows:
                fh.write(f"{chap}\t{src}\t{old}\t{new}\n")


if __name__ == "__main__":
    main()
