#!/usr/bin/env python
"""Chapter composites with one uppercase letter per panel.

A composite cannot carry two lettering levels. Each panel is therefore
relettered from its source PDF text box (the lowercase letter is painted out
and the uppercase letter is drawn on the same baseline), and the letters run
A, B, C ... continuously across the whole chapter plate. Block boundaries are
whitespace only.
"""
import html
import os
import re
import subprocess

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = os.environ.get("LNCPROT_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "manuscript_package", "figures")
OUT = os.path.join(SRC, "assembled")
TMP = os.path.join(ROOT, "figure_build", "work", "reletter")
DPI = 600
MM = DPI / 25.4
FONT_BOLD = "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf"
WORD = re.compile(
    r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>')

os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)


def render(pdf):
    base = os.path.splitext(os.path.basename(pdf))[0]
    png = os.path.join(TMP, base + ".png")
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m",
                    f"-r{DPI}", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                    f"-sOutputFile={png}", pdf], check=True)
    return Image.open(png).convert("RGB")


def panel_boxes(pdf):
    """Single-glyph panel letters with their page coordinates in pixels."""
    html_path = os.path.join(TMP, "bbox.html")
    subprocess.run(["pdftotext", "-bbox", pdf, html_path], check=True)
    txt = open(html_path, encoding="utf-8").read()
    boxes = []
    for x0, y0, x1, y1, s in WORD.findall(txt):
        s = html.unescape(s)
        if len(s) != 1 or s not in "abcdefgh":
            continue
        h = float(y1) - float(y0)
        if not (8.0 <= h <= 11.5) or (float(x1) - float(x0)) < 2.5:
            continue
        boxes.append((s, float(x0), float(y0), float(x1), float(y1)))
    boxes.sort(key=lambda b: b[0])
    return boxes


def reletter(pdf, start_letter, align_top=False):
    """Return (image, mapping) with panel letters replaced from start_letter.

    align_top puts every letter on the highest baseline in the figure, which
    is used where the source figure staggered its letters to clear a two-line
    title.
    """
    im = render(pdf)
    draw = ImageDraw.Draw(im)
    font = ImageFont.truetype(FONT_BOLD, int(9 * DPI / 72))
    boxes = panel_boxes(pdf)
    found = []
    for orig, x0, y0, x1, y1 in boxes:
        px0, py0, px1, py1 = [v * DPI / 72 for v in (x0, y0, x1, y1)]
        pad = 8
        crop = np.asarray(im.crop((int(px0) - pad, int(py0) - pad,
                                   int(px1) + pad, int(py1) + pad)).convert("L"))
        ys, xs = np.where(crop < 200)
        if len(xs):
            ix0 = int(px0) - pad + xs.min()
            iy0 = int(py0) - pad + ys.min()
            ix1 = int(px0) - pad + xs.max()
            iy1 = int(py0) - pad + ys.max()
        else:
            ix0, iy0, ix1, iy1 = int(px0), int(py0), int(px1), int(py1)
        found.append((orig, ix0, iy0, ix1, iy1))
    target = min(f[4] for f in found) if (align_top and found) else None
    mapping = []
    code = ord(start_letter)
    for orig, ix0, iy0, ix1, iy1 in found:
        draw.rectangle([ix0 - 4, iy0 - 4, ix1 + 4, iy1 + 4], fill="white")
        new = chr(code)
        draw.text((ix0, target if target is not None else iy1), new,
                  font=font, fill="black", anchor="ls")
        mapping.append((new, orig))
        code += 1
    if not boxes:
        # single-panel figure: no letter to inherit. Add a white band above the
        # artwork so the new letter cannot sit on a bar or an axis label.
        band = int(9 * MM)
        canvas = Image.new("RGB", (im.width, im.height + band), "white")
        canvas.paste(im, (0, band))
        im = canvas
        draw = ImageDraw.Draw(im)
        new = chr(code)
        draw.text((int(2 * MM), int(6 * MM)), new, font=font, fill="black",
                  anchor="ls")
        mapping.append((new, None))
        code += 1
    return im, mapping


def crop(im):
    a = np.asarray(im.convert("RGB"))
    ink = (a < 250).any(axis=2)
    ys, xs = np.where(ink)
    pad = int(2 * MM)
    return im.crop((max(0, xs.min() - pad), max(0, ys.min() - pad),
                    min(im.width, xs.max() + pad), min(im.height, ys.max() + pad)))


def scale_w(im, w):
    return im.resize((w, int(round(im.height * w / im.width))), Image.LANCZOS)


def compose(blocks, width_mm, cols, name):
    """blocks: list of lists of PDFs; each inner list is stacked in one cell."""
    width = int(width_mm * MM)
    gutter = int(4 * MM)
    col_w = (width - gutter * (cols - 1)) // cols
    cells, mapping = [], []
    letter = "A"
    for bi, block in enumerate(blocks):
        parts = []
        for pdf in block:
            im, mp = reletter(pdf, letter)
            mapping += [(a, b, bi) for a, b in mp]
            letter = chr(ord(letter) + len(mp))
            parts.append(crop(im))
        # stack the parts of one cell
        h = sum(p.height for p in parts) + int(2 * MM) * (len(parts) - 1)
        cell = Image.new("RGB", (max(p.width for p in parts), h), "white")
        y = 0
        for p in parts:
            cell.paste(p, (0, y))
            y += p.height + int(2 * MM)
        cells.append(scale_w(cell, col_w))
    rows = [cells[i:i + cols] for i in range(0, len(cells), cols)]
    heights = [max(c.height for c in r) for r in rows]
    canvas = Image.new("RGB", (width, sum(heights) + gutter * (len(rows) - 1)),
                       "white")
    y = 0
    for r, h in zip(rows, heights):
        x = 0
        for c in r:
            canvas.paste(c, (x, y))
            x += col_w + gutter
        y += h + gutter
    canvas.save(os.path.join(OUT, name + ".png"), dpi=(DPI, DPI))
    canvas.save(os.path.join(OUT, name + ".pdf"), resolution=DPI)
    print(f"{name}: {canvas.width}x{canvas.height} px "
          f"({canvas.width / MM:.0f} x {canvas.height / MM:.0f} mm), "
          f"{len(mapping)} panels {mapping[0][0]}-{mapping[-1][0]}")
    return mapping


def main():
    ch1 = compose([
        [os.path.join(SRC, "Fig1_classes.pdf")],
        [os.path.join(SRC, "Fig2_sequence.pdf")],
        [os.path.join(SRC, "Fig3_TFs.pdf")],
        [os.path.join(SRC, "Fig4_activity.pdf"), os.path.join(SRC, "Fig4d_pair_rho.pdf")],
        [os.path.join(SRC, "Fig5_conservation_chromatin.pdf")],
        [os.path.join(SRC, "Fig6_TE.pdf"), os.path.join(SRC, "Fig3d_TE_age.pdf")],
    ], 183, 2, "ChapterI_origins")

    ch2 = compose([
        [os.path.join(SRC, "Fig7_tcga.pdf")],
        [os.path.join(SRC, "Fig8_gtex.pdf")],
        [os.path.join(SRC, "Fig9_integration.pdf")],
    ], 183, 2, "ChapterII_features")

    ch3 = compose([[os.path.join(SRC, "Fig10_ml.pdf")]], 183, 1,
                  "ChapterIII_function")

    fig_order = {"ChapterI_origins": ["Fig. 1", "Fig. 2", "Fig. 3", "Fig. 4",
                                      "Fig. 5", "Fig. 6"],
                 "ChapterII_features": ["Fig. 7", "Fig. 8", "Fig. 9"],
                 "ChapterIII_function": ["Fig. 10"]}
    # full legend text: each source panel keeps its description, relettered
    src_md = open(os.path.join(ROOT, "SUBMISSION_DRAFT.md"), encoding="utf-8").read()
    lstart = src_md.index("## Figure legends")
    lend = src_md.index("## Ethics statement")
    legends = {}
    for para in src_md[lstart:lend].split("\n\n"):
        mm = re.match(r"\*\*Fig\.\s*(\d+)\s*\|\s*(.*?)\*\*", para.strip(), re.S)
        if not mm:
            continue
        title = mm.group(2).strip().rstrip(".")
        rest = para.strip()[mm.end():]
        # the first mention of each letter opens a panel; later mentions are
        # cross-references inside the text and stay literal
        marks, seen = [], set()
        for m in re.finditer(r"\*\*([a-h])\*\*", rest):
            if m.group(1) in seen:
                continue
            seen.add(m.group(1))
            marks.append(m)
        panels = {}
        for k, m in enumerate(marks):
            stop = marks[k + 1].start() if k + 1 < len(marks) else len(rest)
            body = rest[m.end():stop]
            body = re.sub(r"\*\*([a-h])\*\*", r"\1", body)
            panels[m.group(1)] = body.lstrip(" ,.")
        legends[int(mm.group(1))] = (title, panels)

    title_map = {"ChapterI_origins": "Chapter I | How the promoter is supplied",
                 "ChapterII_features": "Chapter II | Features that the origins might change",
                 "ChapterIII_function": "Chapter III | What the origins predict"}
    out_lines = ["# Chapter composite legends", "",
                 "Each panel carries one uppercase letter, running A, B, C ... "
                 "across the whole plate. Panel descriptions are those of the "
                 "source figures, relettered.", ""]
    for name, maps in (("ChapterI_origins", ch1), ("ChapterII_features", ch2),
                       ("ChapterIII_function", ch3)):
        out_lines.append("## " + title_map[name])
        out_lines.append("")
        current = None
        for new, orig, bi in maps:
            fig_no = int(fig_order[name][bi].split()[-1])
            if orig is None:
                out_lines.append(f"*{legends[fig_no][0]}.* **{new}**, the whole panel.")
                current = None
                continue
            if bi != current:
                out_lines.append("")
                out_lines.append(f"*{legends[fig_no][0]}.*")
                current = bi
            text = legends[fig_no][1].get(orig, "")
            amap = {o: n for n, o, b in maps if b == bi and o}
            # remap the few within-figure cross-references; plain articles
            # such as "a TSS-placement measure" must not be touched
            for x, y in ((x, y) for x in amap for y in amap if x != y):
                text = text.replace(f"in {x} and {y}", f"in {amap[x]} and {amap[y]}")
                text = text.replace(f"Panels {x} and {y}", f"Panels {amap[x]} and {amap[y]}")
                text = text.replace(f"({x}, {y})", f"({amap[x]}, {amap[y]})")
            for x in amap:
                text = text.replace(f"in {x} show", f"in {amap[x]} show")
                text = text.replace(f"in {x} is", f"in {amap[x]} is")
                text = text.replace(f"in {x} mark", f"in {amap[x]} mark")
            out_lines.append(f"**{new}**, {text}")
        out_lines.append("")
    open(os.path.join(OUT, "Chapter_legends.md"), "w", encoding="utf-8").write(
        "\n".join(out_lines))



if __name__ == "__main__":
    main()
