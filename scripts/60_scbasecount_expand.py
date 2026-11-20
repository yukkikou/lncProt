#!/usr/bin/env python
"""60_scbasecount_expand.py — design files for the cross-species extension.

Two designs are produced and written in the same column layout that module 51
consumes, so the verified downloader is reused unchanged.

mouse     every tissue with at least MIN_SAMPLES samples contributes up to
          PER_TISSUE samples spread evenly over the cell-count range. Spreading
          rather than taking the largest avoids a depth bias, because the class
          hierarchy is depth sensitive.
macaca    the only non-human primate with usable sample numbers. Its Ensembl
          annotation covers intergenic lncRNAs only, so this design supports a
          divergent-versus-intergenic comparison and nothing wider.

Run inside the scimilarity container, which is where pyarrow lives.

Outputs
  results/scbasecount/design_mouse_expand.tsv
  results/scbasecount/design_macaque.tsv
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR

ROOT = "/media/share/exchange/1610305236/scBaseCount"
META = os.path.join(ROOT, "metadata", "GeneFull_Ex50pAS")
OUT = os.path.join(WORKDIR, "results", "scbasecount")
MIN_CELLS = 500
MIN_SAMPLES = 20
PER_TISSUE = 4


def norm_tissue(s):
    """Collapse the 1,083 raw tissue labels into comparable names.

    Parenthetical qualifiers, the word mouse and punctuation are removed so that
    'Brain (adult)' and 'brain' land in one group.
    """
    s = str(s).strip().lower()
    s = re.sub(r"\(.*?\)", "", s)
    s = re.sub(r"\bmouse\b", "", s)
    s = re.sub(r"[^a-z ]+", " ", s)
    return " ".join(s.split())


def pick(d, per):
    """Evenly spaced samples over the cell-count range of one tissue."""
    d = d.sort_values("obs_count")
    k = min(per, len(d))
    if len(d) > k:
        idx = np.linspace(0, len(d) - 1, k).astype(int)
        d = d.iloc[idx]
    return d


def design(species, per, thr, label):
    d = pd.read_parquet(os.path.join(META, species,
                                     "sample_metadata.parquet"))
    d = d.rename(columns={"srx_accession": "srx"})
    d["tissue_n"] = d.tissue.map(norm_tissue)
    d = d[d.obs_count >= MIN_CELLS]
    d = d[~d.tissue_n.isin(["", "unsure", "not specified"])]
    vc = d.tissue_n.value_counts()
    keep = set(vc[vc >= thr].index)
    d = d[d.tissue_n.isin(keep)]
    rows = [pick(g, per) for _, g in d.groupby("tissue_n")]
    sel = pd.concat(rows).reset_index(drop=True)
    sel["group"] = sel.tissue_n.map(
        lambda t: "{}_{}".format(label, re.sub(r"[^a-z0-9]+", "_", t).strip("_")))
    sel["arm"] = "all"
    return sel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-tissue", type=int, default=PER_TISSUE)
    ap.add_argument("--min-samples", type=int, default=MIN_SAMPLES)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    with step("60_scbasecount_expand") as st:
        m = design("Mus_musculus", args.per_tissue, args.min_samples, "mouse")
        m = m[["group", "arm", "srx", "file_path", "obs_count"]]
        m.to_csv(safe_write_path(os.path.join(
            OUT, "design_mouse_expand.tsv")), sep="\t", index=False)
        st.note("mouse: {} samples, {} cells, {} tissues".format(
            len(m), int(m.obs_count.sum()), m.group.nunique()))
        st.note("mouse estimated download {:.1f} GB".format(
            m.obs_count.sum() * 15240 / 1e9))

        mm = design("Macaca_mulatta", 2, 10, "macaque")
        mm = mm[["group", "arm", "srx", "file_path", "obs_count"]]
        if len(mm) > 40:
            mm = mm.sample(40, random_state=42)
        mm.to_csv(safe_write_path(os.path.join(
            OUT, "design_macaque.tsv")), sep="\t", index=False)
        st.note("macaque: {} samples, {} cells, {} tissues".format(
            len(mm), int(mm.obs_count.sum()), mm.group.nunique()))
        st.note("macaque estimated download {:.1f} GB".format(
            mm.obs_count.sum() * 15240 / 1e9))


if __name__ == "__main__":
    main()
