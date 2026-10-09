#!/usr/bin/env python
"""04b_lnc_antisense_flag.py — add `antisense_to_lnc` flag (exon overlap with
another lncRNA on the opposite strand, e.g. TSIX vs XIST). Classes stay
defined relative to PCG; this is an interpretation flag only.
Atomic replace of gene_master_table.tsv (safe while other steps read it).
"""
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"


def bt(*args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}", SIF_BT, "bedtools", *args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[:1000])
    return r.stdout


def main():
    with step("04b_antisense_lnc_flag") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        lnc_ids = set(meta[meta.group == "lnc"].index)

        ex = pd.read_csv(os.path.join(DATA, "genes", "exons.bed"), sep="\t",
                         names=["chrom", "start", "end", "gid_tx", "score", "strand"])
        ex["gid"] = ex["gid_tx"].str.split("|").str[0]
        lex = ex[ex["gid"].isin(lnc_ids)].copy()
        lex["sc"] = 0
        p = os.path.join(DATA, "tmp_lnc_ex2.bed")
        lex[["chrom", "start", "end", "gid", "sc", "strand"]].sort_values(
            ["chrom", "start"]).to_csv(safe_write_path(p), sep="\t",
                                       index=False, header=False)
        ov = bt("intersect", "-a", p, "-b", p, "-wa", "-wb", "-S")
        anti = {}
        for ln in ov.splitlines():
            f = ln.split("\t")
            if len(f) < 12:
                continue
            a, b = f[3], f[9]
            if a != b:
                anti.setdefault(a, b)
        meta["antisense_to_lnc"] = [1 if g in anti else 0 for g in meta.index]
        st.note(f"antisense_to_lnc flagged={sum(meta['antisense_to_lnc'])}")
        tsix = meta[meta.gene_name == "TSIX"]
        st.note("TSIX flag=" + str(tsix["antisense_to_lnc"].tolist()))
        tmp = safe_write_path(os.path.join(RES, ".gene_master_table.tmp"))
        meta.to_csv(tmp, sep="\t")
        os.replace(tmp, safe_write_path(os.path.join(RES, "gene_master_table.tsv")))
        st.note("master table updated atomically")


if __name__ == "__main__":
    main()
