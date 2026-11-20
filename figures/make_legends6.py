#!/usr/bin/env python
"""Build six chapter legends from the ten source legends and the panel map."""
import csv
import os
import re

ROOT = os.environ.get("LNCPROT_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAP = os.path.join(ROOT, "manuscript_package", "figures", "chapters", "mapping.tsv")
MD = os.path.join(ROOT, ".backup", "SUBMISSION_DRAFT.v6.0.pre-6chapter.md")
OUT = os.path.join(ROOT, "figure_build", "work", "legends6.md")

TITLES = {
 "Chapter1_architecture": ("Coordinate classes define a promoter-architecture axis", 1),
 "Chapter2_coding_adjacent": ("The coding-adjacent origin borrows the partner promoter", 2),
 "Chapter3_te_host_decayed": ("Transposon-derived, host-embedded and decayed origins", 3),
 "Chapter4_age_constraint": ("Age, constraint and shared regulation", 4),
 "Chapter5_specificity": ("Tissue, cancer and cell specificity", 5),
 "Chapter6_prediction": ("Promoter features carry the axis at its ends", 6),
}


def parse_legends():
    t = open(MD, encoding="utf-8").read()
    ls, le = t.index("## Figure legends"), t.index("## Ethics statement")
    out = {}
    for para in t[ls:le].split("\n\n"):
        m = re.match(r"\*\*Fig\.\s*(\d+)\s*\|\s*(.*?)\*\*", para.strip(), re.S)
        if not m:
            continue
        title = m.group(2).strip().rstrip(".")
        rest = para.strip()[m.end():]
        marks, seen = [], set()
        for mm in re.finditer(r"\*\*([a-h])\*\*", rest):
            if mm.group(1) in seen:
                continue
            seen.add(mm.group(1))
            marks.append(mm)
        panels = {}
        for k, mm in enumerate(marks):
            stop = marks[k + 1].start() if k + 1 < len(marks) else len(rest)
            body = re.sub(r"\*\*([a-h])\*\*", r"\1", rest[mm.end():stop]).lstrip(" ,.")
            panels[mm.group(1)] = body
        out[int(m.group(1))] = (title, panels, rest.strip())
    return out


def main():
    legends = parse_legends()
    rows = list(csv.DictReader(open(MAP), delimiter="\t",
                               fieldnames=["chapter", "fig", "old", "new"]))
    for r in rows:
        if r["fig"] == "Fig. 1B species":
            r["fig"], r["old"] = "Fig. 1", "c"
        elif r["fig"] == "Fig. 1C gradient":
            r["fig"], r["old"] = "Fig. 1", "d"
    # legend paragraphs follow the figure's own letter order (A, B, C ...)
    rows.sort(key=lambda r: (r["chapter"], r["new"]))
    # (source figure, source letter) -> (chapter number, new letter)
    lookup = {}
    for r in rows:
        n = int(r["fig"].split()[-1])
        lookup[(n, r["old"])] = (TITLES[r["chapter"]][1], r["new"])
    # Fig9b is now a real scatter, not the old text box
    lookup[(9, "b")] = (5, "E")
    # source figure -> chapter number (for whole-figure cross references)
    fig2chap = {}
    for (n, _), (c, _) in lookup.items():
        fig2chap.setdefault(n, c)

    def remap(text, src_fig):
        # one pass, with or without a panel letter: a second pass would
        # reinterpret the figure number of an already-remapped panel
        def ref_fig(m):
            n = int(m.group(1))
            x = m.group(2)
            if x:
                c, new = lookup.get((n, x), (fig2chap.get(n, n), x))
                return f"Fig. {c}{new}"
            return f"Fig. {fig2chap.get(n, n)}"
        text = re.sub(r"Fig\.\s*(\d+)([a-h])?", ref_fig, text)
        # within-figure cross references: "in a and c", "(a, b)", "Panels b and e"
        letters = [k[1] for k in lookup if k[0] == src_fig]

        def ref(letter):
            c, new = lookup.get((src_fig, letter), (None, None))
            return f"Fig. {c}{new}" if c else letter
        for x in letters:
            for y in letters:
                if x == y:
                    continue
                rx, ry = ref(x), ref(y)
                if rx[:-1] == ry[:-1]:
                    pair = rx[:-1] + rx[-1] + ", " + ry[-1]
                else:
                    pair = f"{rx} and {ry}"
                text = text.replace(f"in {x} and {y}", f"in {pair}")
                text = text.replace(f"Panels {x} and {y}",
                                    f"Panels {ref(x)} and {ref(y)}")
                rx, ry = ref(x), ref(y)
                if rx[:-1] == ry[:-1]:
                    text = text.replace(f"({x}, {y})",
                                        f"({rx[:-1]}{rx[-1]}, {ry[-1]})")
                else:
                    text = text.replace(f"({x}, {y})", f"({rx}, {ry})")
        for x in letters:
            text = text.replace(f"in {x} show", f"in {ref(x)} show")
            text = text.replace(f"in {x} is", f"in {ref(x)} is")
            text = text.replace(f"in {x} mark", f"in {ref(x)} mark")
        return text

    ch1_src = legends[1][1]
    CH1_PANELS = {
        "a": ch1_src.get("a", ""),
        "b": ch1_src.get("b", ""),
        "c": ("Class shares across six mammals. Stacked bars give the "
              "coordinate-class composition of annotated lncRNAs in human, "
              "mouse and four additional mammals. Mouse reproduces all five "
              "classes; macaque, rat, rabbit and opossum recover only the "
              "divergent and intergenic classes, because their annotations "
              "carry no lncRNA-coding exon overlap."),
        "d": ("Continuous gradient across the axis. Five promoter traits "
              "plotted along the coding-adjacent to self-assembled ordering "
              "(PCG, divergent, antisense, intronic, sense-overlapping, "
              "intergenic): high-CpG fraction, promoter-like chromatin, "
              "CAGE-active fraction, mean phastCons and TE coverage. The "
              "traits fall or rise monotonically along the axis rather than "
              "partitioning the classes into disjoint groups."),
    }

    out = ["## Figure legends", ""]
    for chap in sorted(TITLES, key=lambda k: TITLES[k][1]):
        title, num = TITLES[chap]
        parts = []
        pending = None
        for r in [r for r in rows if r["chapter"] == chap]:
            src_n = int(r["fig"].split()[-1])
            old = r["old"]
            if old == "-":
                txt = legends[src_n][2]
            elif (src_n, old) == (9, "b"):
                txt = ("Scatter of CAGE tissue specificity against RNA-seq "
                       "tissue specificity for the 19,730 genes with both "
                       "measurements, coloured by class; the dashed line is "
                       "y = x. Spearman ρ = 0.72.")
            elif src_n == 1:
                txt = CH1_PANELS.get(old, "")
            else:
                txt = legends[src_n][1].get(old, "")
            body = remap(txt, src_n).strip()
            if pending is not None:
                letter, prev = pending
                parts.append(f"**{letter}, {r['new']}**, {prev} {body}")
                pending = None
                continue
            if body.endswith(" and"):
                pending = (r["new"], body)
                continue
            parts.append(f"**{r['new']}**, {body}")
        if pending is not None:
            letter, prev = pending
            parts.append(f"**{letter}**, {prev}")
        out.append(f"**Fig. {num} | {title}.** " + " ".join(parts))
        out.append("")
    open(OUT, "w", encoding="utf-8").write("\n".join(out))
    print("wrote", OUT, len("\n".join(out)), "chars")


if __name__ == "__main__":
    main()
