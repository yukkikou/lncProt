#!/usr/bin/env python
"""43_constraint_structure.py — promoter vs gene-body/exon constraint,
conserved-patch structure, and TE promoter-vs-exon dual role.

Inputs:
  results/constraint_structure/gene_profiles_raw.tsv (pybw_gene_profiles.py)
  data/genes/exons.bed, Dfam TE BED, te_promoter_stats.tsv
Outputs -> results/constraint_structure/*.tsv + step summary.
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "constraint_structure")
RAW = os.path.join(OUT, "gene_profiles_raw.tsv")
EXONS = os.path.join(WORKDIR, "data", "genes", "exons.bed")
TE = ("/media/share/resource/reference/homo_sapiens/GRCh38/annotation/"
      "hg38.nrph.TE.anno.sorted.bed")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def mwu(y, x):
    u, p = mannwhitneyu(y, x, alternative="two-sided")
    return float(2 * u / (len(y) * len(x)) - 1), float(p)


def sing(*args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}"] + list(args)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(" ".join(cmd) + "\n" + r.stderr[:800])
    return r.stdout


def te_overlap(exons, te, out):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}", SIF_BT, "bedtools",
           "intersect", "-a", exons, "-b", te, "-wo"]
    with open(safe_write_path(out), "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)


def main():
    with step("43_constraint") as st:
        raw = pd.read_csv(RAW, sep="\t", index_col=0, low_memory=False)
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)[
            ["cls", "group"]]
        m = meta.join(raw, how="inner")
        m["cls"] = np.where(m.group == "pcg", "PCG", m.cls)
        pcg = m[m.group == "pcg"]
        st.note("genes with profiles: {} (PCG {})".format(
            len(m), len(pcg)))

        # ---------------- A. promoter vs body/exon/intron --------------
        feats = ["prom_phast", "body_phast", "exon_phast", "intron_phast",
                 "p5_phast", "p3_phast", "prom_phylop", "body_phylop",
                 "exon_phylop", "intron_phylop", "p5_phylop", "p3_phylop"]
        m["prom_minus_body_phast"] = m.prom_phast - m.body_phast
        m["prom_minus_body_phylop"] = m.prom_phylop - m.body_phylop
        m["p5_minus_p3_phast"] = m.p5_phast - m.p3_phast
        m["p5_minus_p3_phylop"] = m.p5_phylop - m.p3_phylop
        feats += ["prom_minus_body_phast", "prom_minus_body_phylop",
                  "p5_minus_p3_phast", "p5_minus_p3_phylop"]
        pcg = m[m.group == "pcg"]
        rows = []
        for f in feats:
            x = pcg[f].dropna().values
            for c in CLASSES + ["pseudogene"]:
                y = m.loc[m.cls == c, f].dropna().values
                if len(y) < 20:
                    continue
                d, p = mwu(y, x)
                rows.append(dict(feature=f, cls=c, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=d, pval=p))
        A = pd.DataFrame(rows)
        A["fdr"] = bh_fdr(A.pval.values)
        A.to_csv(safe_write_path(os.path.join(
            OUT, "A_conservation_tests.tsv")), sep="\t", index=False)

        # ---------------- B. conserved patch structure -----------------
        pfeats = [c for c in m.columns if "_patch_" in c]
        rows = []
        for f in pfeats:
            x = pcg[f].dropna().values
            for c in CLASSES + ["pseudogene"]:
                y = m.loc[m.cls == c, f].dropna().values
                if len(y) < 20:
                    continue
                d, p = mwu(y, x)
                rows.append(dict(feature=f, cls=c, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=d, pval=p))
        B = pd.DataFrame(rows)
        B["fdr"] = bh_fdr(B.pval.values)
        B.to_csv(safe_write_path(os.path.join(
            OUT, "B_patch_tests.tsv")), sep="\t", index=False)

        # ---------------- C. TE promoter vs exon ----------------------
        te_ltr = os.path.join(WORKDIR, "data", "constraint_structure_tmp",
                              "TE_LTR.bed")
        os.makedirs(os.path.dirname(te_ltr), exist_ok=True)
        with open(TE) as fi, open(safe_write_path(te_ltr), "w") as fo:
            for line in fi:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 8 and f[7] == "LTR":
                    fo.write(line)
        ov_all = os.path.join(WORKDIR, "data", "constraint_structure_tmp",
                              "exon_TE_ov.tsv")
        ov_ltr = os.path.join(WORKDIR, "data", "constraint_structure_tmp",
                              "exon_LTR_ov.tsv")
        te_overlap(EXONS, TE, ov_all)
        te_overlap(EXONS, te_ltr, ov_ltr)

        def agg_ov(path):
            d = {}
            with open(path) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 15:
                        continue
                    gid = f[3].split("|")[0]
                    d[gid] = d.get(gid, 0) + int(f[-1])
            return pd.Series(d)

        s_all = agg_ov(ov_all)
        s_ltr = agg_ov(ov_ltr)
        m["te_exon_bp"] = s_all.reindex(m.index).fillna(0).values
        m["ltr_exon_bp"] = s_ltr.reindex(m.index).fillna(0).values
        m["te_exon_frac"] = (m.te_exon_bp /
                             m.exon_len.replace(0, np.nan)).clip(upper=1)
        m["ltr_exon_frac"] = (m.ltr_exon_bp /
                              m.exon_len.replace(0, np.nan)).clip(upper=1)
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["te_cov", "te_frac_LTR", "tss_in_ltr"]]
        m = m.join(te, how="left")
        pcg = m[m.group == "pcg"]
        rows = []
        for f in ["te_exon_frac", "ltr_exon_frac", "te_cov", "te_frac_LTR"]:
            x = pcg[f].dropna().values
            for c in CLASSES + ["pseudogene"]:
                y = m.loc[m.cls == c, f].dropna().values
                if len(y) < 20:
                    continue
                d, p = mwu(y, x)
                rows.append(dict(feature=f, cls=c, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=d, pval=p))
        ig = m[m.cls == "intergenic"]
        for f in ["te_exon_frac", "ltr_exon_frac"]:
            a = ig.loc[ig.tss_in_ltr == 1, f].dropna().values
            b = ig.loc[ig.tss_in_ltr != 1, f].dropna().values
            if len(a) >= 20 and len(b) >= 20:
                d, p = mwu(a, b)
                rows.append(dict(feature=f, cls="intergenic_LTR_vs_nonLTR",
                                 n=len(a) + len(b),
                                 median_class=float(np.median(a)),
                                 median_PCG=float(np.median(b)),
                                 cliffs_delta=d, pval=p))
        C = pd.DataFrame(rows)
        C["fdr"] = bh_fdr(C.pval.values)
        C.to_csv(safe_write_path(os.path.join(
            OUT, "C_te_dual_tests.tsv")), sep="\t", index=False)
        m[["cls", "group", "exon_len", "te_exon_bp", "ltr_exon_bp",
           "te_exon_frac", "ltr_exon_frac", "te_cov", "te_frac_LTR"]].to_csv(
            safe_write_path(os.path.join(OUT, "C_te_per_gene.tsv")),
            sep="\t")

        ok = m[["te_frac_LTR", "ltr_exon_frac"]].dropna()
        r, pr = spearmanr(ok.te_frac_LTR, ok.ltr_exon_frac)
        st.note("A: prom-body phyloP divergent d={:.3f}; intergenic "
                "d={:.3f}".format(
                    A[(A.feature == "prom_minus_body_phylop") &
                      (A.cls == "divergent")].cliffs_delta.iloc[0],
                    A[(A.feature == "prom_minus_body_phylop") &
                      (A.cls == "intergenic")].cliffs_delta.iloc[0]))
        st.note("B: prom patch frac phylop intergenic d={:.3f}".format(
            B[(B.feature == "prom_patch_frac_phylop") &
              (B.cls == "intergenic")].cliffs_delta.iloc[0]))
        st.note("C: LTR promoter vs LTR exon rho={:.3f} p={:.1e}".format(
            r, pr))
        st.summary_md("results/constraint_structure/step_summaries/43.md",
                      "43 constraint structure (promoter vs body, patches, TE)",
                      "## A. Promoter vs body/exon/intron\n\n"
                      + df_to_md(A, ".3f", max_rows=60)
                      + "\n\n## B. Patch metrics\n\n"
                      + df_to_md(B, ".3f", max_rows=60)
                      + "\n\n## C. TE promoter vs exonization\n\n"
                      + df_to_md(C, ".3f", max_rows=40)
                      + "\n\nLTR promoter vs LTR exon ρ={:.3f}, "
                        "P={:.2e}".format(r, pr))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
