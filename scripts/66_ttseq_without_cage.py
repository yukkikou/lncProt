#!/usr/bin/env python
"""66_ttseq_without_cage.py — is there initiation without a CAGE peak?

The Results say divergent lncRNAs borrow a partner promoter, yet only 62% of
their promoters carry a CAGE peak while 98% of their partners do. If the
borrowing model is right, the remaining 38% should still be transcribed, just
below the detection threshold of a cap-based assay. TT-seq measures synthesis
directly and does not require a capped 5' end, so it can separate the two
readings.

TT-seq here is K562, so the comparison uses the ENCODE K562 CAGE peaks from the
same cell type rather than the multi-tissue FANTOM set. Same cell type, same
locus, two assays.

Outputs
  results/ttseq/cage_vs_synthesis.tsv
  results/step_summaries/66_ttseq_without_cage.md
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
CAGE = os.path.join(WORKDIR, "data", "strand_controls",
                    "ENCODE_K562_CAGE_idr.bed.gz")
WIN = 250
CLS = ["divergent", "antisense", "intronic", "sense_overlapping",
       "intergenic"]


def load_cage():
    rows = []
    with gzip.open(CAGE, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            strand = f[5] if len(f) >= 6 else "+"
            mid = (int(f[1]) + int(f[2])) // 2
            rows.append((f[0], mid, strand))
    return pd.DataFrame(rows, columns=["chrom", "pos", "strand"])


def norm(c):
    return c if str(c).startswith("chr") else "chr" + str(c)


def main():
    with step("66_ttseq_without_cage") as st:
        tt = pd.read_csv(os.path.join(RES, "ttseq", "TU_mapped.tsv"),
                         sep="\t", low_memory=False)
        st.note("TT-seq units mapped to genes={}".format(len(tt)))
        cage = load_cage()
        st.note("ENCODE K562 CAGE peaks={}".format(len(cage)))

        by_chr = {}
        for c, g in cage.groupby("chrom"):
            for s, gg in g.groupby("strand"):
                by_chr[(c, s)] = np.sort(gg.pos.values)

        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        tt["chrom"] = tt.chrom.map(norm)
        hit = []
        for r in tt.itertuples():
            arr = by_chr.get((r.chrom, r.strand))
            if arr is None or not np.isfinite(r.tss):
                hit.append(False)
                continue
            j = np.searchsorted(arr, r.tss)
            ok = False
            for k in (j - 1, j, j + 1):
                if 0 <= k < len(arr) and abs(int(arr[k]) - int(r.tss)) <= WIN:
                    ok = True
                    break
            hit.append(ok)
        tt["k562_cage"] = hit
        st.note("TT-seq units with a K562 CAGE peak at the TSS={:.1%}".format(
            tt.k562_cage.mean()))

        # The TU table only lists genes with a detected unit, so its rows cannot
        # serve as a denominator. The universe is every annotated locus of the
        # class, and a missing TU is counted as no synthesis, which is the
        # conservative reading.
        have_tu = set(tt.gene_id)
        rows = []
        for c in CLS + ["PCG"]:
            uni = m[(m.cls == c) & m.tss.notna()]
            if len(uni) < 20:
                continue
            tss = uni.tss.values.astype(int)
            chrom = uni.chrom.map(norm).values
            strand = uni["strand"].values
            is_cage = np.zeros(len(uni), dtype=bool)
            for i in range(len(uni)):
                arr = by_chr.get((chrom[i], strand[i]))
                if arr is None:
                    continue
                j = np.searchsorted(arr, tss[i])
                for k in (j - 1, j, j + 1):
                    if 0 <= k < len(arr) and abs(int(arr[k]) - tss[i]) <= WIN:
                        is_cage[i] = True
                        break
            has_tu = np.array([g in have_tu for g in uni.index])
            a = int((is_cage & has_tu).sum())
            b = int((is_cage & ~has_tu).sum())
            d = int((~is_cage & has_tu).sum())
            e = int((~is_cage & ~has_tu).sum())
            p = fisher_exact([[a, b], [d, e]])[1]
            rows.append(dict(cls=c, n_loci=len(uni), cage_and_tu=a,
                             cage_no_tu=b, nocage_tu=d, nocage_no_tu=e,
                             frac_cage_with_tu=float(a / max(a + b, 1)),
                             frac_nocage_with_tu=float(d / max(d + e, 1)),
                             pval=float(p)))
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            RES, "ttseq", "cage_vs_synthesis.tsv")), sep="\t", index=False)
        for r in T.itertuples():
            st.note("{}: {} loci, TU in {:.0%} of CAGE+ and {:.0%} of CAGE- "
                    "loci, P={:.2g}".format(
                        r.cls, r.n_loci, r.frac_cage_with_tu,
                        r.frac_nocage_with_tu, r.pval))

        # how strong is synthesis without a CAGE peak
        div = tt[tt.cls == "divergent"]
        a = div[div.k562_cage].syn_per_kb.dropna()
        b = div[~div.k562_cage].syn_per_kb.dropna()
        if len(a) > 20 and len(b) > 20:
            w, pw = mannwhitneyu(a, b, alternative="two-sided")
            st.note("divergent synthesis per kb: CAGE+ median {:.4f} "
                    "(n={}), CAGE- median {:.4f} (n={}), P={:.3g}".format(
                        a.median(), len(a), b.median(), len(b), pw))
            pd.DataFrame([dict(cls="divergent", group="cage_pos", n=len(a),
                               median_syn_per_kb=float(a.median())),
                          dict(cls="divergent", group="cage_neg", n=len(b),
                               median_syn_per_kb=float(b.median()))
                          ]).to_csv(safe_write_path(os.path.join(
                              RES, "ttseq", "divergent_synthesis_by_cage.tsv")),
                              sep="\t", index=False)

        st.summary_md(os.path.join(RES, "step_summaries",
                                   "66_ttseq_without_cage.md"),
                      "66 TT-seq synthesis at CAGE-negative divergent promoters",
                      df_to_md(T, ".4f"))


if __name__ == "__main__":
    main()
