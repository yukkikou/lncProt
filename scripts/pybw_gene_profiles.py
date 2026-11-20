#!/usr/bin/env python
"""pybw_gene_profiles.py — run INSIDE pyBigWig.sif.

Per-gene promoter / gene-body / exon / intron / 5'-half / 3'-half mean
phastCons and phyloP, plus conserved-patch metrics in the promoter and in
the 5'/3' 2-kb body windows.
Reads genes.tsv, exons.bed, promoters_feat.bed; writes
results/constraint_structure/gene_profiles_raw.tsv.
"""
import math
import os

import numpy as np
import pyBigWig

W = "/media/share/workdir/1610305236/lncProt"
PHAST = "/media/share/resource/database/PhyloP_score/hg38.phastCons470way.bw"
PHYLOP = "/media/share/resource/database/PhyloP_score/241-mammalian-2020v2.bigWig"
GENES = f"{W}/data/genes/genes.tsv"
EXONS = f"{W}/data/genes/exons.bed"
PROM = f"{W}/results/promoters_feat.bed"
OUT = f"{W}/results/constraint_structure/gene_profiles_raw.tsv"


def stat(bw, chrom, s, e, typ):
    if e <= s:
        return float("nan")
    try:
        r = bw.stats(chrom, s, e, type=typ, exact=True)
    except Exception:
        return float("nan")
    if not r or r[0] is None:
        return float("nan")
    return float(r[0])


def patch(bw, chrom, s, e, thr, ge=True):
    if e <= s:
        return 0, 0, 0.0, 0.0, 0
    try:
        v = np.asarray(bw.values(chrom, s, e), dtype=float)
    except Exception:
        return 0, 0, 0.0, 0.0, 0
    ok = ~np.isnan(v)
    m = ok & ((v >= thr) if ge else (v <= thr))
    if not m.any():
        return 0, 0, 0.0, 0.0, 0
    z = np.concatenate(([0], m.astype(np.int8), [0]))
    d = np.diff(z)
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    lens = ends - starts
    n = len(m)
    h = n // 2
    f5, f3 = (m[:h].sum(), m[h:].sum()) if strand_plus else \
        (m[h:].sum(), m[:h].sum())
    total = int(m.sum())
    bias = float(f5 - f3) / max(total, 1)
    return int(len(lens)), int(lens.max()), float(total) / n, bias, total


def merge_iv(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if not out or s > out[-1][1]:
            out.append([s, e])
        else:
            out[-1][1] = max(out[-1][1], e)
    return out


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    bw_p = pyBigWig.open(PHAST)
    bw_l = pyBigWig.open(PHYLOP)

    prom = {}
    with open(PROM) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            prom[f[3]] = (f[0], int(f[1]), int(f[2]), f[5])

    ex = {}
    with open(EXONS) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gid = f[3].split("|")[0]
            ex.setdefault(gid, []).append((int(f[1]), int(f[2])))

    header = ["gene_id", "prom_len", "body_len", "exon_len", "intron_len"]
    for t in ["phast", "phylop"]:
        header += [f"prom_{t}", f"body_{t}", f"exon_{t}", f"intron_{t}",
                   f"p5_{t}", f"p3_{t}"]
        header += [f"prom_patch_n_{t}", f"prom_patch_max_{t}",
                   f"prom_patch_frac_{t}", f"prom_patch_bias5_{t}"]
        header += [f"b5_patch_n_{t}", f"b5_patch_max_{t}",
                   f"b5_patch_frac_{t}", f"b5_patch_bias5_{t}"]
        header += [f"b3_patch_n_{t}", f"b3_patch_max_{t}",
                   f"b3_patch_frac_{t}", f"b3_patch_bias5_{t}"]
    header += ["prom_accel_n", "prom_accel_frac"]

    out = open(OUT, "w")
    out.write("\t".join(header) + "\n")
    global strand_plus
    n = 0
    with open(GENES) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gid, chrom, start, end, strand = f[0], f[4], int(f[5]), int(f[6]), f[7]
            s, e = start - 1, end
            if e <= s:
                continue
            strand_plus = (strand == "+")
            if strand_plus:
                p5s, p5e = s, min(s + 2000, e)
                p3s, p3e = max(e - 2000, s), e
            else:
                p5s, p5e = max(e - 2000, s), e
                p3s, p3e = s, min(s + 2000, e)
            mid = (s + e) // 2
            h5s, h5e, h3s, h3e = ((s, mid, mid, e) if strand_plus
                                  else (mid, e, s, mid))
            plen = (prom[gid][2] - prom[gid][1]) if gid in prom else 0
            row = [gid, plen, e - s, 0, 0]
            vals = {}
            for t, bw in [("phast", bw_p), ("phylop", bw_l)]:
                bmean = stat(bw, chrom, s, e, "mean")
                bsum = stat(bw, chrom, s, e, "sum")
                vals[f"body_{t}"] = bmean
                vals[f"p5_{t}"] = stat(bw, chrom, h5s, h5e, "mean")
                vals[f"p3_{t}"] = stat(bw, chrom, h3s, h3e, "mean")
                # exons
                eiv = merge_iv(ex.get(gid, []))
                elen = sum(x[1] - x[0] for x in eiv)
                esum = 0.0
                for a, b in eiv:
                    x = stat(bw, chrom, a, b, "sum")
                    if not math.isnan(x):
                        esum += x
                vals[f"exon_{t}"] = esum / elen if elen > 0 else float("nan")
                ilen = (e - s) - elen
                vals[f"intron_{t}"] = ((bsum - esum) / ilen
                                       if ilen > 0 and not math.isnan(bsum)
                                       else float("nan"))
                # promoter
                if gid in prom:
                    pc, ps, pe, pstrand = prom[gid]
                    pm = stat(bw, pc, ps, pe, "mean")
                    vals[f"prom_{t}"] = pm
                    pn, pmax, pfrac, pbias, _ = patch(
                        bw, pc, ps, pe, 0.5 if t == "phast" else 2.0)
                else:
                    vals[f"prom_{t}"] = float("nan")
                    pn = pmax = 0
                    pfrac = pbias = 0.0
                vals[f"prom_patch_n_{t}"] = pn
                vals[f"prom_patch_max_{t}"] = pmax
                vals[f"prom_patch_frac_{t}"] = pfrac
                vals[f"prom_patch_bias5_{t}"] = pbias
                # body 5'/3' patches
                n5, m5, f5_, b5, _ = patch(
                    bw, chrom, p5s, p5e, 0.5 if t == "phast" else 2.0)
                n3, m3, f3_, b3, _ = patch(
                    bw, chrom, p3s, p3e, 0.5 if t == "phast" else 2.0)
                for pre, nn, mm, ff, bb in [("b5", n5, m5, f5_, b5),
                                            ("b3", n3, m3, f3_, b3)]:
                    vals[f"{pre}_patch_n_{t}"] = nn
                    vals[f"{pre}_patch_max_{t}"] = mm
                    vals[f"{pre}_patch_frac_{t}"] = ff
                    vals[f"{pre}_patch_bias5_{t}"] = bb
            # accelerated patches in promoter (phyloP <= -2)
            if gid in prom:
                pc, ps, pe, _ = prom[gid]
                an, amax, afrac, _, _ = patch(bw_l, pc, ps, pe, -2.0,
                                              ge=False)
            else:
                an = amax = 0
                afrac = 0.0
            for t in ["phast", "phylop"]:
                row += [vals.get(f"prom_{t}", float("nan")),
                        vals.get(f"body_{t}", float("nan")),
                        vals.get(f"exon_{t}", float("nan")),
                        vals.get(f"intron_{t}", float("nan")),
                        vals.get(f"p5_{t}", float("nan")),
                        vals.get(f"p3_{t}", float("nan"))]
                row += [vals[f"prom_patch_n_{t}"],
                        vals[f"prom_patch_max_{t}"],
                        vals[f"prom_patch_frac_{t}"],
                        vals[f"prom_patch_bias5_{t}"],
                        vals[f"b5_patch_n_{t}"], vals[f"b5_patch_max_{t}"],
                        vals[f"b5_patch_frac_{t}"],
                        vals[f"b5_patch_bias5_{t}"],
                        vals[f"b3_patch_n_{t}"], vals[f"b3_patch_max_{t}"],
                        vals[f"b3_patch_frac_{t}"],
                        vals[f"b3_patch_bias5_{t}"]]
                row += [an, afrac] if t == "phylop" else []
            # exon/intron lengths: rebuild for row (exon_len currently 0)
            eiv = merge_iv(ex.get(gid, []))
            elen = sum(x[1] - x[0] for x in eiv)
            row[2] = e - s
            row[3] = elen
            row[4] = (e - s) - elen
            out.write("\t".join(str(x) for x in row) + "\n")
            n += 1
    out.close()
    bw_p.close()
    bw_l.close()
    print(f"gene profiles done: {n}")


if __name__ == "__main__":
    main()
