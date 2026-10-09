#!/usr/bin/env python
"""pybw_stats.py — run INSIDE pyBigWig.sif (stdlib + pyBigWig only).
Reads promoters_feat.bed + phastCons470way.bw; writes per-promoter mean/max
phastCons over the -500/+200 window. Handles chr-prefix normalization.
Output written to a path under WORKDIR only.
"""
import sys

import pyBigWig

WORKDIR = "/media/share/workdir/1610305236/lncProt"
BW = "/media/share/resource/database/PhyloP_score/hg38.phastCons470way.bw"
BED = f"{WORKDIR}/results/promoters_feat.bed"
OUT = f"{WORKDIR}/results/conservation_raw.tsv"


def main():
    bw = pyBigWig.open(BW)
    names = set(bw.chroms().keys())
    strip = not any(n.startswith("chr") for n in names)
    out = open(OUT, "w")
    out.write("gene_id\tphast_mean\tphast_max\n")
    n_missing = 0
    with open(BED) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, s, e, gid = f[0], int(f[1]), int(f[2]), f[3]
            c = chrom[3:] if (strip and chrom.startswith("chr")) else chrom
            if c not in bw.chroms():
                n_missing += 1
                out.write(f"{gid}\tnan\tnan\n")
                continue
            stats_m = bw.stats(c, s, e, type="mean")
            stats_x = bw.stats(c, s, e, type="max")
            m = stats_m[0] if stats_m and stats_m[0] is not None else float("nan")
            x = stats_x[0] if stats_x and stats_x[0] is not None else float("nan")
            out.write(f"{gid}\t{m}\t{x}\n")
    out.close()
    print(f"phastCons stats done; missing chroms={n_missing}")


if __name__ == "__main__":
    main()
