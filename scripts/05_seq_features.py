#!/usr/bin/env python
"""05_seq_features.py — promoter sequence features per gene/class.

Input : results/promoters_feat.fa (-500/+200), data/cpgIslandExt.txt.gz
Output: results/seq_features.tsv, results/step_summaries/05_seq_features.md
"""
import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR, df_to_md
from pylib.seqfeatures import (gc_content, cpg_oe, cpg_density, hcp_icp_lcp,
                               scan_core_elements)

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"


def bt(*args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}", SIF_BT, "bedtools", *args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bedtools failed: {' '.join(args)}\n{r.stderr[:1000]}")
    return r.stdout


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main():
    with step("05_seq_features") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                           index_col=0)

        # ---- sequence features from promoter FASTA
        rows = []
        for name, seq in read_fasta(os.path.join(RES, "promoters_feat.fa")):
            gc = gc_content(seq)
            oe = cpg_oe(seq)
            cpg_n = seq.count("CG")
            rec = dict(gene_id=name, gc=gc, cpg_oe=oe, cpg_density=cpg_density(seq),
                       cpg_count=cpg_n, prom_class=hcp_icp_lcp(gc, oe, cpg_n),
                       length=len(seq))
            rec.update(scan_core_elements(seq))
            rows.append(rec)
        feat = pd.DataFrame(rows).set_index("gene_id")
        st.note(f"promoters scored={len(feat)}")

        # ---- CpG island overlap/distance via bedtools
        rows = []
        with gzip.open(os.path.join(DATA, "cpgIslandExt.txt.gz"), "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    rows.append((f[1], int(f[2]), int(f[3]), f[4]))
        rows.sort(key=lambda r: (r[0], r[1]))   # lexicographic for bedtools
        with open(safe_write_path(os.path.join(DATA, "cpg.bed")), "w") as out:
            for c, s, e, n in rows:
                out.write(f"{c}\t{s}\t{e}\t{n}\n")
        cl = bt("closest", "-a", os.path.join(RES, "promoters_feat.bed"),
                "-b", os.path.join(DATA, "cpg.bed"), "-d", "-t", "first")
        cpg_hit, cpg_dist = {}, {}
        for ln in cl.splitlines():
            f = ln.split("\t")
            if len(f) < 11:
                continue
            gid = f[3]
            d = int(f[10]) if f[6] != "." else 10**9
            cpg_dist[gid] = d
            cpg_hit[gid] = 1 if d == 0 else 0
        feat["cpg_island"] = [cpg_hit.get(g, np.nan) for g in feat.index]
        feat["cpg_island_dist"] = [cpg_dist.get(g, np.nan) for g in feat.index]

        out = meta.join(feat, how="left")
        out.to_csv(safe_write_path(os.path.join(RES, "seq_features.tsv")), sep="\t")

        # ---- summary by class
        lnc_classes = ["divergent", "antisense", "intronic", "sense_overlapping",
                       "intergenic"]
        groups = ["PCG"] + lnc_classes + ["pseudogene", "transcribed_pseudogene"]
        rows = []
        for g in groups:
            sub = out[out["cls"] == g] if g != "PCG" else out[out["group"] == "pcg"]
            if not len(sub):
                continue
            rows.append(dict(cls=g, n=len(sub),
                             gc=sub["gc"].mean(), cpg_oe=sub["cpg_oe"].mean(),
                             HCP=(sub["prom_class"] == "HCP").mean(),
                             LCP=(sub["prom_class"] == "LCP").mean(),
                             TATA=sub["TATA"].mean(), Inr=sub["Inr"].mean(),
                             CCAAT=sub["CCAAT"].mean(), GCbox=sub["GCbox"].mean(),
                             cpg_island=sub["cpg_island"].mean()))
        sm = pd.DataFrame(rows)
        sm.to_csv(safe_write_path(os.path.join(RES, "seq_features_by_class.tsv")),
                  sep="\t", index=False)
        body = df_to_md(sm, ".3f")
        st.note("class medians: " + "; ".join(
            f"{r.cls}: GC={r.gc:.2f} o/e={r.cpg_oe:.2f} TATA={r.TATA:.2f} HCP={r.HCP:.2f}"
            for r in sm.itertuples()))
        st.summary_md("results/step_summaries/05_seq_features.md",
                      "05 promoter sequence features", body)


if __name__ == "__main__":
    main()
