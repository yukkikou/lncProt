#!/usr/bin/env python
"""pybw_stats_phylop.py — run INSIDE pyBigWig.sif (stdlib + pyBigWig only).
Reads promoters_feat.bed + 241-mammalian-2020v2.bigWig (Cactus/Zoonomia
phyloP); writes per-promoter base-level phyloP stats over the same
-500/+200 window used for phastCons (pybw_stats.py): mean/max/min,
covered-base count, and fraction of covered bases at >=+2 / <=-2 (plus
1.3 sensitivity). Denominators are covered (non-None) bases only.
chr-prefix normalization mirrors pybw_stats.py.
"""
import math
import os

import pyBigWig

WORKDIR = "/media/share/workdir/1610305236/lncProt"
BW = "/media/share/resource/database/PhyloP_score/241-mammalian-2020v2.bigWig"
BED = f"{WORKDIR}/results/promoters_feat.bed"
OUTDIR = f"{WORKDIR}/results/phylop"
OUT = f"{OUTDIR}/phylop_raw.tsv"


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    bw = pyBigWig.open(BW)
    names = set(bw.chroms().keys())
    strip = not any(n.startswith("chr") for n in names)
    out = open(OUT, "w")
    out.write("gene_id\tphylop_mean\tphylop_max\tphylop_min\tn_cov"
              "\tfrac_cons_2\tfrac_accel_2\tfrac_cons_13\tfrac_accel_13\n")
    n_missing = 0
    n_zero = 0
    with open(BED) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, s, e, gid = f[0], int(f[1]), int(f[2]), f[3]
            c = chrom[3:] if (strip and chrom.startswith("chr")) else chrom
            if c not in names:
                n_missing += 1
                out.write(f"{gid}\tnan\tnan\tnan\t0\tnan\tnan\tnan\tnan\n")
                continue
            try:
                raw = bw.values(c, s, e)
            except RuntimeError:
                n_missing += 1
                out.write(f"{gid}\tnan\tnan\tnan\t0\tnan\tnan\tnan\tnan\n")
                continue
            vals = [v for v in raw
                    if v is not None and not math.isnan(v)]
            if not vals:
                n_zero += 1
                out.write(f"{gid}\tnan\tnan\tnan\t0\tnan\tnan\tnan\tnan\n")
                continue
            n = len(vals)
            mean = sum(vals) / n
            vmax = max(vals)
            vmin = min(vals)
            fc2 = sum(1 for v in vals if v >= 2.0) / n
            fa2 = sum(1 for v in vals if v <= -2.0) / n
            fc13 = sum(1 for v in vals if v >= 1.3) / n
            fa13 = sum(1 for v in vals if v <= -1.3) / n
            out.write(f"{gid}\t{mean}\t{vmax}\t{vmin}\t{n}"
                      f"\t{fc2}\t{fa2}\t{fc13}\t{fa13}\n")
    out.close()
    bw.close()
    print(f"phyloP stats done; missing chroms={n_missing}; "
          f"zero-coverage windows={n_zero}")


if __name__ == "__main__":
    main()
