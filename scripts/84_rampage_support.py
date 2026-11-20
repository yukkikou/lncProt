#!/usr/bin/env python
"""84_rampage_support.py — independent ENCODE RAMPAGE support at the fixed TSS.

The manuscript calls a promoter active when a FANTOM5 peak used for TSS
refinement carries TPM > 1, which is circular. This module leaves the TSS set
untouched and asks whether an independent cap-based dataset, ENCODE RAMPAGE,
supports the same TSS. ENCODE RAMPAGE peaks were not used to define or refine
any TSS here.

For each gene the module records whether any of the RAMPAGE peak sets has a
peak within +/-250 bp of the major TSS, and how many peak sets do so. It then
compares class support fractions and the lncRNA-versus-partner asymmetry in
divergent pairs.

Outputs
  results/revision_checks/H_rampage_support.tsv
  results/revision_checks/H_rampage_pair.tsv
"""
import glob
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "revision_checks")
BED = os.path.join(WORKDIR, "data", "encode_rampage", "bed")
HALF = 250
CLASSES = ["PCG", "divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic", "pseudogene", "transcribed_pseudogene"]


def load_bed(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            c = f[0].replace("chr", "")
            out.setdefault(c, []).append((int(f[1]), int(f[2])))
    merged = {}
    for c, v in out.items():
        v.sort()
        s, e = [], []
        for a, b in v:
            if s and a <= e[-1]:
                e[-1] = max(e[-1], b)
            else:
                s.append(a); e.append(b)
        merged[c] = (np.array(s), np.array(e))
    return merged


def main():
    with step("84_rampage_support") as st:
        os.makedirs(OUT, exist_ok=True)
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["tss"] = pd.to_numeric(m.tss, errors="coerce")
        m["chrom_c"] = m.chrom.astype(str).str.replace("chr", "", regex=False)
        m = m.dropna(subset=["tss"]).reset_index()
        by_chr = {c: np.where(m.chrom_c == c)[0] for c in
                  set(m.chrom_c)}
        files = sorted(glob.glob(os.path.join(BED, "*.bed.gz")))
        st.note("RAMPAGE peak sets={}".format(len(files)))
        n_hit = np.zeros(len(m), dtype=np.int32)
        for k, path in enumerate(files):
            iv = load_bed(path)
            for c, idx in by_chr.items():
                t = iv.get(c)
                if t is None:
                    continue
                s, e = t
                q = m.tss.values[idx].astype(np.int64)
                j = np.searchsorted(s, q + HALF, "right") - 1
                ok = np.zeros(len(idx), dtype=bool)
                valid = j >= 0
                ok[valid] = e[j[valid]] >= (q[valid] - HALF)
                n_hit[idx[ok]] += 1
            if (k + 1) % 50 == 0:
                st.note("processed {} peak sets".format(k + 1))
        m["rampage_sets"] = n_hit
        m["rampage_any"] = (n_hit > 0).astype(int)
        rows = []
        for cls in CLASSES:
            sub = m[m.cls == cls]
            if len(sub) < 20:
                continue
            rows.append(dict(cls=cls, n=len(sub),
                             frac_any=float(sub.rampage_any.mean()),
                             median_sets=float(sub.rampage_sets.median()),
                             mean_sets=float(sub.rampage_sets.mean())))
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "H_rampage_support.tsv")), sep="\t", index=False)

        # divergent pair asymmetry
        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        base_tss = m.set_index(m.gene_id.astype(str).str.split(".").str[0])
        dp["lnc_hit"] = dp.lnc_gid.astype(str).str.split(".").str[0].map(
            base_tss.rampage_any)
        dp["pcg_hit"] = dp.pcg_gid.astype(str).str.split(".").str[0].map(
            base_tss.rampage_any)
        dp = dp.dropna(subset=["lnc_hit", "pcg_hit"])
        a = int(((dp.lnc_hit == 1) & (dp.pcg_hit == 0)).sum())
        b = int(((dp.lnc_hit == 0) & (dp.pcg_hit == 1)).sum())
        p = binomtest(min(a, b), a + b, 0.5).pvalue
        P = pd.DataFrame([dict(n_pairs=len(dp),
                               lnc_support=float(dp.lnc_hit.mean()),
                               pcg_support=float(dp.pcg_hit.mean()),
                               lnc_only=a, pcg_only=b,
                               mcnemar_binom_p=float(p))])
        P.to_csv(safe_write_path(os.path.join(
            OUT, "H_rampage_pair.tsv")), sep="\t", index=False)
        for r in T.itertuples():
            st.note("{}: RAMPAGE support {:.1%} (median sets {:.0f})".format(
                r.cls, r.frac_any, r.median_sets))
        st.note("pairs: lnc {:.1%} v partner {:.1%}; lnc-only {} pcg-only {} "
                "P={:.2g}".format(P.lnc_support.iloc[0],
                                  P.pcg_support.iloc[0], a, b, p))
        st.summary_md(os.path.join(
            RES, "step_summaries", "84_rampage_support.md"),
            "84 ENCODE RAMPAGE support at fixed TSSs",
            df_to_md(T, ".3f") + "\n\n" + df_to_md(P, ".4f"))


if __name__ == "__main__":
    main()
