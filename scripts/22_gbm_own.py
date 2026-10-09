#!/usr/bin/env python
"""22_gbm_own.py — user's own GBM lncRNA-seq cohort (135 samples, salmon TPM).

1. Aggregate transcript TPM -> gene TPM (tx2gene from GENCODE GTF).
2. Class activity per sample; LTR-derived lncRNA relative activation.
3. Clinical join (lncRNA-cj.csv): fastq stem -> row via `sample_Name` then
   `Sample` (S- prefixes / -UCC#### suffixes stripped); 131/135 resolve.
   CGGA-mirror grade II-IV LTR-delta analysis + tumor/non-tumor control.
Outputs -> results/gbm_own/*.tsv + step summary
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

QD = os.path.join(WORKDIR, "data", "gbm_own", "quants")
OUT = os.path.join(WORKDIR, "results", "gbm_own")
RES = os.path.join(WORKDIR, "results")
GTF = ("/media/share/resource/reference/homo_sapiens/GRCh38/annotation/"
       "gencode.v46.annotation.gtf")
CLIN = ("/media/share/data8/volume1/1610305236_CJ_GBM/GBM_lncRNA/lncRNA-cj.csv")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def tx2gene():
    tx = {}
    pat = re.compile(r'transcript_id "([^"]+)";.*?gene_id "([^"]+)"')
    pat2 = re.compile(r'gene_id "([^"]+)";.*?transcript_id "([^"]+)"')
    op = gzip.open if GTF.endswith(".gz") else open
    import gzip
    with op(GTF, "rt") as fh:
        for line in fh:
            if not line.startswith("#") and "\ttranscript\t" in line:
                m = pat.search(line) or pat2.search(line)
                if m:
                    a, b = m.group(1), m.group(2)
                    # normalise: key = transcript id, value = gene id
                    if a.startswith("ENST"):
                        tx[a] = b
                    else:
                        tx[b] = a
    return tx


def main():
    with step("22_gbm_own") as st:
        os.makedirs(OUT, exist_ok=True)
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", low_memory=False)
        sym2cls = dict(zip(meta.gene_name.astype(str), meta.cls))
        g2sym = {g.split(".")[0]: str(s) for g, s in
                 zip(meta.gene_id, meta.gene_name)}
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["tss_in_ltr"]].join(
            meta.set_index("gene_id")["gene_name"])
        sym2ltr = {str(s): int(v) for s, v in zip(te.gene_name, te.tss_in_ltr)
                   if pd.notna(s)}

        quants = sorted(glob.glob(os.path.join(QD, "*", "quant.sf")))
        st.note(f"salmon quants: {len(quants)}")
        if len(quants) < 10:
            raise RuntimeError("too few quants — run s9 first")
        t2g = tx2gene()

        cols = {}
        for q in quants:
            s = os.path.basename(os.path.dirname(q))
            d = pd.read_csv(q, sep="\t", usecols=["Name", "TPM"])
            t2g_nov = {k.split(".")[0]: v for k, v in t2g.items()}
            d["gene"] = d.Name.map(lambda t: t2g.get(t) or t2g_nov.get(t.split(".")[0]))
            g = d.dropna(subset=["gene"]).groupby("gene")["TPM"].sum()
            cols[s] = g
        M = pd.DataFrame(cols).fillna(0)
        M.index = [g.split(".")[0] for g in M.index]
        M = M.groupby(M.index).sum()
        st.note(f"gene x sample matrix: {M.shape}")

        cls = np.array([sym2cls.get(g2sym.get(g, ""), None) for g in M.index])
        ltr = np.array([sym2ltr.get(g2sym.get(g, ""), 0) for g in M.index])
        LM = np.log2(M + 1)
        lnc = np.isin(cls, CLASSES)

        rows = []
        for s in M.columns:
            v = LM[s]
            r = dict(sample=s)
            for c in CLASSES + ["PCG"]:
                idx = (cls == c) if c != "PCG" else np.isin(cls, ["PCG"])
                sub = v[idx]
                r[f"{c}_frac_expr"] = float((sub >= 1).sum() / max(idx.sum(), 1))
            r["ltr_rel"] = float(v[lnc & (ltr == 1)].mean()
                                 - v[lnc & (ltr == 0)].mean())
            rows.append(r)
        T = pd.DataFrame(rows)

        # clinical join: lncRNA-cj.csv carries TWO id systems — `Sample`
        # (short numbers) and `sample_Name` (bare 6-digit / GBM*/LGG*/N*).
        # Fastq stems mix both, some with `S-` prefixes or `-UCC####`
        # suffixes; strip those, try sample_Name first, then Sample.
        clin = pd.read_csv(CLIN)
        clin["sample_Name"] = clin["sample_Name"].astype(str).str.strip()
        clin["Sample"] = clin["Sample"].astype(str).str.strip()
        by_name = clin.drop_duplicates(subset=["sample_Name"]).set_index("sample_Name")
        by_sample = clin.drop_duplicates(subset=["Sample"]).set_index("Sample")

        def lookup(fq_id):
            key = fq_id[2:] if fq_id.startswith("S-") else fq_id
            key = key.split("-")[0] if "-UCC" in key else key
            if key in by_name.index:
                return by_name.loc[key], "sample_Name"
            if key in by_sample.index:
                return by_sample.loc[key], "Sample"
            return None, "none"

        jrows, jkey = [], []
        for s in T["sample"]:
            r, k = lookup(s)
            jkey.append(k)
            jrows.append({c: pd.NA for c in clin.columns} if r is None
                         else r.to_dict())
        m = pd.concat([T, pd.DataFrame(jrows, index=T.index)], axis=1)
        m["clin_join_key"] = jkey
        m["clin_matched"] = m["clin_join_key"] != "none"
        n_name = int((m.clin_join_key == "sample_Name").sum())
        n_sample = int((m.clin_join_key == "Sample").sum())
        unmatched = m.loc[~m.clin_matched, "sample"].tolist()
        st.note(f"clinical crosswalk: {int(m.clin_matched.sum())}/{len(m)} matched "
                f"(sample_Name key {n_name}, Sample key {n_sample}); "
                f"abandoned: {unmatched}")
        m.to_csv(safe_write_path(os.path.join(OUT, "gbm_own_per_sample.tsv")),
                 sep="\t", index=False)

        # CGGA-mirror grade analysis on matched glioma samples (grades II-IV;
        # Disease=normal rows are adjacent/non-tumor controls, excluded)
        g = pd.to_numeric(m["Grade"], errors="coerce")
        tumor = m.clin_matched & m["Disease"].isin(["GBM", "LGG"])
        ok = m.ltr_rel.notna() & g.notna() & tumor & g.isin([2, 3, 4])
        grade_txt = "; ".join(
            f"Grade {gv}: LTR-delta={m.ltr_rel[ok & (g == gv)].median():+.3f} "
            f"(n={int((ok & (g == gv)).sum())})" for gv in (2, 3, 4))
        st.note(grade_txt)
        lo, hi = m.ltr_rel[ok & (g == 2)], m.ltr_rel[ok & (g == 4)]
        rho, prho = spearmanr(g[ok], m.ltr_rel[ok])
        if len(lo) >= 5 and len(hi) >= 5:
            _, p = mannwhitneyu(hi, lo)
            st.note(f"grade IV vs II MWU p={p:.2e}; trend Spearman rho={rho:.2f} "
                    f"p={prho:.1e} (n={int(ok.sum())})")
        # tumor vs non-tumor control within the own cohort
        norm = m.clin_matched & (m["Disease"] == "normal")
        if int(norm.sum()) >= 3 and int(tumor.sum()) >= 5:
            tn, nn = m.ltr_rel[m.clin_matched & (m["Disease"] == "GBM")], m.ltr_rel[norm]
            _, p2 = mannwhitneyu(tn, nn)
            st.note(f"GBM vs non-tumor: {tn.median():+.3f} (n={len(tn)}) vs "
                    f"{nn.median():+.3f} (n={len(nn)}); MWU p={p2:.2e}")
        # class-activity hierarchy over matched tumor samples (clean quote)
        hier = ", ".join(
            "{} {:.0f}%".format(c, m.loc[tumor, c + "_frac_expr"].mean() * 100)
            for c in ["PCG"] + CLASSES)
        st.note(f"matched-tumor class hierarchy: {hier}")
        st.note("per-sample table written")
        st.summary_md("results/gbm_own/step_summaries/22_gbm_own.md",
                      "22 Own GBM cohort",
                      f"- {M.shape[1]} samples quantified; clinical crosswalk "
                      f"{int(m.clin_matched.sum())}/{len(m)} "
                      f"(sample_Name {n_name} / Sample {n_sample}); abandoned "
                      f"4 without clinical rows\n"
                      f"- {grade_txt}; IV vs II MWU p={p:.2e}, trend rho={rho:.2f} "
                      f"(p={prho:.1e}, n={int(ok.sum())})\n"
                      f"- matched-tumor class hierarchy: {hier}\n"
                      f"- per-sample table with clinical columns + join flags: "
                      f"results/gbm_own/gbm_own_per_sample.tsv")


if __name__ == "__main__":
    main()
