#!/usr/bin/env python
"""pybw_patch_thresholds.py — run INSIDE pyBigWig.sif.

Promoter conserved/accelerated patch metrics at multiple phyloP/phastCons
thresholds (sensitivity for module 43/45).
"""
import os

import numpy as np
import pyBigWig

W = "/media/share/workdir/1610305236/lncProt"
PHAST = "/media/share/resource/database/PhyloP_score/hg38.phastCons470way.bw"
PHYLOP = "/media/share/resource/database/PhyloP_score/241-mammalian-2020v2.bigWig"
PROM = f"{W}/results/promoters_feat.bed"
OUT = f"{W}/results/constraint_structure/patch_thresholds.tsv"

PHYLOP_CONS = [1.3, 2.0]
PHYLOP_ACC = [1.3, 2.0]
PHAST_CONS = [0.3, 0.5, 0.7]


def runs(v, thr, ge):
    m = (~np.isnan(v)) & ((v >= thr) if ge else (v <= thr))
    if not m.any():
        return 0, 0, 0.0
    z = np.concatenate(([0], m.astype(np.int8), [0]))
    d = np.diff(z)
    st = np.where(d == 1)[0]
    en = np.where(d == -1)[0]
    lens = en - st
    return int(len(lens)), int(lens.max()), float(m.sum()) / len(m)


def main():
    bw_p = pyBigWig.open(PHAST)
    bw_l = pyBigWig.open(PHYLOP)
    hdr = ["gene_id"]
    for t in PHYLOP_CONS:
        hdr += [f"phylop_ge{t}_n", f"phylop_ge{t}_max",
                f"phylop_ge{t}_frac"]
    for t in PHYLOP_ACC:
        hdr += [f"phylop_le-{t}_n", f"phylop_le-{t}_max",
                f"phylop_le-{t}_frac"]
    for t in PHAST_CONS:
        hdr += [f"phast_ge{t}_n", f"phast_ge{t}_max", f"phast_ge{t}_frac"]
    out = open(OUT, "w")
    out.write("\t".join(hdr) + "\n")
    n = 0
    with open(PROM) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, s, e, gid = f[0], int(f[1]), int(f[2]), f[3]
            row = [gid]
            try:
                vp = np.asarray(bw_l.values(chrom, s, e), dtype=float)
            except Exception:
                vp = np.array([np.nan])
            try:
                vh = np.asarray(bw_p.values(chrom, s, e), dtype=float)
            except Exception:
                vh = np.array([np.nan])
            for t in PHYLOP_CONS:
                row += list(runs(vp, t, True))
            for t in PHYLOP_ACC:
                row += list(runs(vp, -t, False))
            for t in PHAST_CONS:
                row += list(runs(vh, t, True))
            out.write("\t".join(str(x) for x in row) + "\n")
            n += 1
    out.close()
    bw_p.close()
    bw_l.close()
    print(f"patch thresholds done: {n}")


if __name__ == "__main__":
    main()
