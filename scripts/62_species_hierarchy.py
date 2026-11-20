#!/usr/bin/env python
"""62_species_hierarchy.py — class detection in any species' single-cell set.

Generic over species. The class map comes from module 59, which builds the
five-class taxonomy from a local Ensembl GTF, and the matrices come from the
scBaseCount manifest of module 60. Because the non-model Ensembl annotations
contain no overlapping loci (see the fourth limitation), only the classes that
exist in that species' annotation are reported; nothing is imputed.

Outputs
  results/cross_species/<species>/sample_detection.tsv
  results/cross_species/<species>/hierarchy.tsv
  results/cross_species/<species>/step_summary.md
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
sys.modules.setdefault("m52", m52)

ROOT = "/media/share/exchange/1610305236/scBaseCount"
OUT = os.path.join(WORKDIR, "results", "cross_species")
PRESENT = ["PCG", "divergent", "intergenic"]


def species_codes(species):
    g = pd.read_csv(os.path.join(OUT, species, "gene_master_table.tsv"),
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
    ap.add_argument("species")
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    sp = args.species
    man_path = args.manifest or os.path.join(ROOT, "manifest_{}.tsv".format(sp))
    out = os.path.join(OUT, sp)
    os.makedirs(out, exist_ok=True)
    with step("62_{}_hierarchy".format(sp)) as st:
        man = pd.read_csv(man_path, sep="\t")
        man = man[[os.path.exists(p) for p in man.dest]]
        st.note("samples present locally={} cells={}".format(
            len(man), int(man.cells.sum())))
        if not len(man):
            st.note("no local matrices, nothing to do")
            return
        code_map = species_codes(sp)
        st.note("class map genes={} PCG={} divergent={} intergenic={}".format(
            len(code_map), int((code_map == m52.CODES["PCG"]).sum()),
            int((code_map == m52.CODES["divergent"]).sum()),
            int((code_map == m52.CODES["intergenic"]).sum())))
        rows = [dict(r) for _, r in man.iterrows()]
        if args.workers > 1:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(args.workers) as ex:
                res = list(ex.map(m52.one_sample, rows,
                                  [code_map] * len(rows), chunksize=2))
        else:
            res = [m52.one_sample(r, code_map) for r in rows]
        S = pd.DataFrame(res)
        S.to_csv(safe_write_path(os.path.join(out, "sample_detection.tsv")),
                 sep="\t", index=False)
        cols = ["{}_det".format(c) for c in PRESENT]
        R = S[cols].rank(axis=1, ascending=False)
        order = R.apply(lambda r: ">".join(
            r.sort_values().index.str.replace("_det", "")), axis=1)
        H = pd.DataFrame({
            "class": PRESENT,
            "mean_rank": [R[c].mean() for c in cols],
            "median_det": [S[c].median() for c in cols],
            "n_samples": [int(S[c].notna().sum()) for c in cols]})
        H.to_csv(safe_write_path(os.path.join(out, "hierarchy.tsv")),
                 sep="\t", index=False)
        st.note("samples={} cells={}".format(len(S),
                                             int(S.n_cells_file.sum())))
        st.note("modal order: {} ({}/{})".format(
            order.value_counts().index[0], order.value_counts().iloc[0],
            len(order)))
        st.note("PCG rank 1 in {:.1%}, divergent above intergenic in {:.1%}"
                .format((R.PCG_det == 1).mean(),
                        (R.divergent_det < R.intergenic_det).mean()))
        st.summary_md(os.path.join(out, "step_summary.md"),
                      "62 {} class hierarchy".format(sp),
                      df_to_md(H, ".4f") + "\n\n"
                      + df_to_md(order.value_counts().head(5)
                                 .rename("n").reset_index(), "4"))


if __name__ == "__main__":
    main()
