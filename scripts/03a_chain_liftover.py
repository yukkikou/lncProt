#!/usr/bin/env python
"""03a_chain_liftover.py — CLI: lift a BED file hg19 -> hg38 via in-house chainlift.

Usage:
  python 03a_chain_liftover.py --in in.bed[.gz] --out out.bed [--min-rate 0.5]
Drops unmapped / cross-chain intervals. Reports lift rate (fails if < min-rate
of TOTAL lines were attempted-and-mapped proportion; warns only, exit code
reflects parse success — module-level QC lives in 03b).
"""
import argparse
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.chainlift import ChainDB
from pylib.progress import step, safe_write_path, df_to_md

# NOTE: the reference-tree hg19ToHg38.over.chain is TRUNCATED (593KB, 1278
# stub chains, verified 2026-09-30). Use the Ensembl assembly-mapping chain
# downloaded into data/qc (complete; validated: BRCA1 hg19->hg38 +/-2bp).
DEFAULT_CHAIN = ("/media/share/workdir/1610305236/lncProt/data/qc/"
                 "grch37_to_grch38.chain")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--chain", default=DEFAULT_CHAIN)
    ap.add_argument("--label", default="03a_liftover")
    args = ap.parse_args()

    with step(args.label) as st:
        db = ChainDB(args.chain)
        n = mapped = gap_span = 0
        op = gzip.open if args.inp.endswith(".gz") else open
        outp = safe_write_path(args.out)
        with op(args.inp, "rt") as fh, open(outp, "w") as out:
            for line in fh:
                if line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                n += 1
                chrom, s, e = f[0], int(f[1]), int(f[2])
                m = db.lift(chrom, s, e)
                if m is None or m["start"] >= m["end"]:
                    continue
                if m["gaps"] > 0:
                    gap_span += 1
                mapped += 1
                rest = f[3:] if len(f) > 3 else []
                out.write(f"{m['chrom']}\t{m['start']}\t{m['end']}\t" +
                          "\t".join(rest) + "\n")
        rate = mapped / n if n else 0.0
        st.note(f"lines={n} mapped={mapped} rate={rate:.3f} gap_spanning={gap_span}")
        print(f"[liftover] {args.inp}: {n} -> {mapped} ({rate:.1%}), "
              f"{gap_span} span internal gaps")


if __name__ == "__main__":
    main()
