#!/usr/bin/env python
"""49_cage_version_control.py — is the CAGE plus-strand bias pipeline-specific?

Compares plus fraction of:
  FANTOM5 hg38 fair+new peaks
  FANTOM5 hg38 liftover+new peaks (alternative processing)
  ENCODE K562 CAGE IDR peaks (independent pipeline)
and plus fraction in shared vs unique peaks (bedtools intersect).
Outputs -> results/strand_controls/*.tsv + step summary.
"""
import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "strand_controls")
SC = os.path.join(WORKDIR, "data", "strand_controls")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"
SETS = {
    "FANTOM_fair": os.path.join(WORKDIR, "data", "fantom",
                                "hg38_cage_peaks.bed.gz"),
    "FANTOM_liftover": os.path.join(SC,
                                    "hg38_liftover+new_CAGE_peaks.bed.gz"),
    "ENCODE_K562": os.path.join(SC, "ENCODE_K562_CAGE_idr.bed.gz"),
}


def read_bed(path):
    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rows.append((f[0], int(f[1]), int(f[2]), f[4], f[5]))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "score",
                                       "strand"])


def count(p):
    pos = int((p.strand == "+").sum())
    neg = int((p.strand == "-").sum())
    return pos, neg


def bt(args):
    r = subprocess.run(["singularity", "exec", "--no-mount", "bind-paths",
                        "--bind", "/media/share/resource:/media/share/resource:ro",
                        "--bind", f"{WORKDIR}:{WORKDIR}", SIF_BT,
                        "bedtools"] + args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[:500])
    return r.stdout


def main():
    with step("49_cage_versions") as st:
        os.makedirs(OUT, exist_ok=True)
        frames = {k: read_bed(v) for k, v in SETS.items()}
        rows = []
        for k, p in frames.items():
            pos, neg = count(p)
            rows.append(dict(set=k, n=len(p), plus=pos, minus=neg,
                             plus_frac=pos / len(p),
                             p_vs_half=binomtest(pos, pos + neg,
                                                 0.5).pvalue,
                             median_len=float((p.end - p.start).median())))
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "cage_version_plus.tsv")), sep="\t", index=False)

        # pairwise shared vs unique
        tmp = os.path.join(SC, "tmp")
        os.makedirs(tmp, exist_ok=True)

        def write_bed(name, p):
            path = os.path.join(tmp, name + ".bed")
            with open(safe_write_path(path), "w") as fh:
                for r in p.itertuples():
                    fh.write("{}\t{}\t{}\t.\t0\t{}\n".format(
                        r.chrom, r.start, r.end, r.strand))
            return path

        paths = {k: write_bed(k, p) for k, p in frames.items()}
        rows = []
        for a, b in [("FANTOM_fair", "FANTOM_liftover"),
                     ("FANTOM_fair", "ENCODE_K562"),
                     ("FANTOM_liftover", "ENCODE_K562")]:
            shared = bt(["intersect", "-a", paths[a], "-b", paths[b], "-u"])
            only = bt(["intersect", "-a", paths[a], "-b", paths[b], "-v"])
            for lab, txt in [("shared", shared), ("unique_" + a, only)]:
                arr = [x.split("\t")[5] for x in txt.splitlines() if x]
                pos = sum(1 for x in arr if x == "+")
                neg = len(arr) - pos
                if neg + pos == 0:
                    continue
                rows.append(dict(pair=a + "_vs_" + b, subset=lab,
                                 n=pos + neg, plus_frac=pos / (pos + neg),
                                 p_vs_half=binomtest(pos, pos + neg,
                                                     0.5).pvalue))
        U = pd.DataFrame(rows)
        U.to_csv(safe_write_path(os.path.join(
            OUT, "cage_overlap_plus.tsv")), sep="\t", index=False)

        st.note("version plus: " + "; ".join(
            "{} {:.3f} (n={})".format(r.set, r.plus_frac, r.n)
            for r in T.itertuples()))
        st.note("overlap plus: " + "; ".join(
            "{} {} {:.3f}".format(r.pair, r.subset, r.plus_frac)
            for r in U.itertuples() if r.subset == "shared"))
        st.summary_md("results/strand_controls/step_summaries/49.md",
                      "49 CAGE version controls",
                      "## Overall plus fraction\n\n" + df_to_md(T, ".3f")
                      + "\n\n## Shared vs unique peaks\n\n"
                      + df_to_md(U, ".3f", max_rows=30))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
