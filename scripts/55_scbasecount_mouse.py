#!/usr/bin/env python
"""55_scbasecount_mouse.py — birth-class detection in mouse single cells.

Reuses the per-sample scan of module 52 with the mouse class annotation from
module 37 (GRCm38.102, same coordinate rules as the human table). The mouse
atlas shares 98.4% of its gene identifiers with that annotation, so the class
definitions are directly comparable.

The question is whether the class detection hierarchy is preserved across
species. Only ranks and relative ordering are compared, never absolute levels,
because the two annotations differ in completeness.

Outputs
  results/scbasecount/mouse_sample_detection.tsv
  results/scbasecount/mouse_hierarchy.tsv
  results/scbasecount/step_summaries/55.md
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

_spec = importlib.util.spec_from_file_location(
    "m52", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "52_scbasecount_classes.py"))
m52 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m52)
# register under its own name so worker processes can unpickle one_sample
sys.modules.setdefault("m52", m52)

OUT = os.path.join(WORKDIR, "results", "scbasecount")
ROOT = "/media/share/exchange/1610305236/scBaseCount"
MAN = os.path.join(ROOT, "manifest_mouse.tsv")


def mouse_codes():
    g = pd.read_csv(os.path.join(WORKDIR, "results", "mouse_taxonomy",
                                 "gene_master_table.tsv"),
                    sep="\t", low_memory=False).drop_duplicates("gene_id")
    g = g.set_index("gene_id")
    code = pd.Series(0, index=g.index, dtype=np.int8)
    cls = g["cls"].astype(str)
    code[cls.isin(m52.LNC5)] = cls[cls.isin(m52.LNC5)].map(m52.CODES)
    code[cls == "pseudogene"] = m52.CODES["pseudogene"]
    code[g["group"] == "pcg"] = m52.CODES["PCG"]
    return code.astype(np.int8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    with step("55_scbasecount_mouse") as st:
        os.makedirs(OUT, exist_ok=True)
        man = pd.read_csv(MAN, sep="\t")
        code_map = mouse_codes()
        rows = [dict(r) for _, r in man.iterrows()]
        if args.workers > 1:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(args.workers) as ex:
                res = list(ex.map(m52.one_sample, rows,
                                  [code_map] * len(rows), chunksize=4))
        else:
            res = [m52.one_sample(r, code_map) for r in rows]
        S = pd.DataFrame(res)
        S["tissue"] = S.group.str.replace("mouse_", "", regex=False)
        S.to_csv(safe_write_path(os.path.join(
            OUT, "mouse_sample_detection.tsv")), sep="\t", index=False)

        det = [c for c in S.columns
               if c.endswith("_det") and c.split("_det")[0] in
               ("PCG", "divergent", "antisense", "intronic",
                "sense_overlapping", "intergenic")]
        ranks = S[det].rank(axis=1, ascending=False)
        order = S[det].apply(lambda r: ">".join(
            r.sort_values(ascending=False).index.str.replace("_det", "")),
            axis=1)
        H = pd.DataFrame({
            "class": [d.replace("_det", "") for d in det],
            "mean_rank": [ranks[d].mean() for d in det],
            "n_samples": [int(S[d].notna().sum()) for d in det]})
        H.to_csv(safe_write_path(os.path.join(
            OUT, "mouse_hierarchy.tsv")), sep="\t", index=False)
        st.note("samples={} cells={}".format(len(S), int(S.n_cells_file.sum())))
        st.note("modal order: " + order.value_counts().index[0]
                + " ({}/{})".format(order.value_counts().iloc[0], len(order)))
        st.summary_md(os.path.join(OUT, "step_summaries", "55.md"),
                      "55 mouse class hierarchy",
                      df_to_md(H, ".3f") + "\n\n"
                      + df_to_md(order.value_counts().head(6)
                                 .rename("n").reset_index(), "4"))


if __name__ == "__main__":
    main()
