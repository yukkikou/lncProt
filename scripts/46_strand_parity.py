#!/usr/bin/env python
"""46_strand_parity.py — DNA strand parity and its functional breaking.

Modes:
  A balance     — class / CAGE-peak plus-minus counts
  B te          — TSS-donating TE orientation, starts, family effects
  C length      — human gene/exon/TE length parity + power-law exponents
  D skew        — promoter GC/AT skew by class (Chargaff-related composition)
  E transcription — divergent-pair and LTR sense/antisense CAGE asymmetry
  F cross       — gene/exon length parity and exponents across 7 species
Outputs -> results/strand_parity/*.tsv + step summary.
"""
import gzip
import os
import sys
from bisect import bisect_right

import numpy as np
import pandas as pd
from scipy.stats import binomtest, ks_2samp, spearmanr
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "strand_parity")
DATA = os.path.join(WORKDIR, "data")
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


def pl_alpha(x):
    """Clauset-style continuous power-law alpha with KS-chosen xmin."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x) & (x > 0)]
    if len(x) < 100:
        return np.nan, np.nan, np.nan
    best = None
    for q in np.arange(0.5, 0.96, 0.05):
        xmin = float(np.quantile(x, q))
        tail = np.sort(x[x >= xmin])
        if len(tail) < 50 or xmin <= 0:
            continue
        alpha = 1 + len(tail) / np.sum(np.log(tail / xmin))
        cdf = np.arange(1, len(tail) + 1) / len(tail)
        theo = 1 - (tail / xmin) ** (1 - alpha)
        ks = float(np.max(np.abs(cdf - theo)))
        if best is None or ks < best[0]:
            best = (ks, alpha, xmin)
    return (best[1], best[2], best[0]) if best else (np.nan, np.nan, np.nan)


def main():
    mode = (sys.argv[1] if len(sys.argv) > 1 else "A").upper()
    meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                       sep="\t", index_col=0, low_memory=False)
    with step("46_parity_" + mode) as st:
        os.makedirs(OUT, exist_ok=True)
        if mode == "A":
            rows = []
            for c in CLASSES + ["pseudogene"]:
                sub = meta[meta.cls == c]
                pos = int((sub.strand == "+").sum())
                neg = int((sub.strand == "-").sum())
                p = binomtest(pos, pos + neg, 0.5).pvalue
                rows.append(dict(class_=c, n=pos + neg, plus=pos, minus=neg,
                                 plus_frac=pos / max(pos + neg, 1), pval=p))
            pos = int((meta[meta.group == "pcg"].strand == "+").sum())
            neg = int((meta[meta.group == "pcg"].strand == "-").sum())
            rows.append(dict(class_="PCG", n=pos + neg, plus=pos, minus=neg,
                             plus_frac=pos / max(pos + neg, 1),
                             pval=binomtest(pos, pos + neg, 0.5).pvalue))
            # CAGE peaks
            pk = gzip.open(os.path.join(DATA, "fantom",
                                        "hg38_cage_peaks.bed.gz"), "rt")
            pp = pn = 0
            for line in pk:
                f = line.split("\t")
                if len(f) > 5:
                    if f[5] == "+":
                        pp += 1
                    elif f[5] == "-":
                        pn += 1
            pk.close()
            rows.append(dict(class_="CAGE_peaks", n=pp + pn, plus=pp,
                             minus=pn, plus_frac=pp / max(pp + pn, 1),
                             pval=binomtest(pp, pp + pn, 0.5).pvalue))
            T = pd.DataFrame(rows)
            T["fdr"] = bh_fdr(T.pval.values)
            T.to_csv(safe_write_path(os.path.join(OUT, "A_balance.tsv")),
                     sep="\t", index=False)
            st.note("balance: " + "; ".join(
                "{} {:.3f}".format(r.class_, r.plus_frac) for r in T.itertuples()))
            st.summary_md("results/strand_parity/step_summaries/46_a.md",
                          "46A strand balance", df_to_md(T, ".3f"))

        elif mode == "B":
            tt = pd.read_csv(os.path.join(DATA, "tss_in_te.tsv"), sep="\t",
                             header=None,
                             names=["c", "s", "e", "gene_id", "tc", "ts",
                                    "te", "tid", "t_sub", "tstrand", "fam",
                                    "order"])
            tt = tt.merge(meta[["strand", "cls", "tss"]], left_on="gene_id",
                          right_index=True, how="inner")
            tt["same"] = (tt.tstrand == tt.strand)
            rows = []
            for lab, sub in [("all_TE_TSS", tt),
                             ("LTR_TSS", tt[tt.order == "LTR"]),
                             ("nonLTR_TSS", tt[tt.order != "LTR"])]:
                n = len(sub)
                k = int(sub.same.sum())
                p = binomtest(k, n, 0.5).pvalue
                rows.append(dict(contrast=lab, n=n, same=k, opposite=n - k,
                                 same_frac=k / max(n, 1), pval=p))
                for c in CLASSES:
                    s2 = sub[sub.cls == c]
                    if len(s2) < 20:
                        continue
                    k2 = int(s2.same.sum())
                    rows.append(dict(contrast=lab + "_" + c, n=len(s2),
                                     same=k2, opposite=len(s2) - k2,
                                     same_frac=k2 / len(s2),
                                     pval=binomtest(k2, len(s2), 0.5).pvalue))
            T = pd.DataFrame(rows)
            T["fdr"] = bh_fdr(T.pval.values)
            T.to_csv(safe_write_path(os.path.join(OUT, "B_te_orientation.tsv")),
                     sep="\t", index=False)
            # TE start offsets relative to TSS, by strand
            off = np.where(tt.strand == "+", tt.ts - tt.tss,
                           tt.tss - tt.te)
            tt = tt.assign(offset=off)
            off_rows = []
            for lab, sub in [("same", tt[tt.same]), ("opposite", tt[~tt.same])]:
                v = sub.offset.dropna().values
                if len(v) > 20:
                    off_rows.append(dict(group=lab, n=len(v),
                                         median=float(np.median(v)),
                                         mean=float(np.mean(v))))
            pd.DataFrame(off_rows).to_csv(safe_write_path(os.path.join(
                OUT, "B_te_offsets.tsv")), sep="\t", index=False)
            st.note("TE orientation: " + "; ".join(
                "{} {:.3f}".format(r.contrast, r.same_frac)
                for r in T.itertuples() if r.contrast in
                ("all_TE_TSS", "LTR_TSS", "nonLTR_TSS")))
            st.summary_md("results/strand_parity/step_summaries/46_b.md",
                          "46B TE orientation", df_to_md(T, ".3f"))

        elif mode == "C":
            genes = pd.read_csv(os.path.join(DATA, "genes", "genes.tsv"),
                                sep="\t")
            genes["len"] = genes.end - genes.start + 1
            genes = genes.merge(meta[["cls"]], left_on="gene_id",
                                right_index=True, how="left")
            rows = []
            feats = [("gene", genes)]
            ex = pd.read_csv(os.path.join(DATA, "genes", "exons.bed"),
                             sep="\t", header=None,
                             names=["chrom", "start", "end", "gid", "x",
                                    "strand"])
            ex["gene_id"] = ex.gid.str.split("|").str[0]
            ex["len"] = ex.end - ex.start
            ex = ex.merge(meta[["cls", "group"]], left_on="gene_id",
                          right_index=True, how="left")
            feats.append(("exon", ex))
            # TE lengths from BED
            te = pd.read_csv("/media/share/resource/reference/homo_sapiens/"
                             "GRCh38/annotation/hg38.nrph.TE.anno.sorted.bed",
                             sep="\t", header=None,
                             names=["c", "s", "e", "tid", "t_sub", "strand",
                                    "fam", "order"])
            te["len"] = te.e - te.s
            feats.append(("TE", te))
            for name, df in feats:
                for cls in (["PCG"] + CLASSES + ["pseudogene"] if name != "TE"
                            else ["all"]):
                    sub = (df[df.group == "pcg"] if cls == "PCG" else
                           (df[df.cls == cls] if cls != "all" else df))
                    a = sub.loc[sub.strand == "+", "len"].dropna().values
                    b = sub.loc[sub.strand == "-", "len"].dropna().values
                    if len(a) < 50 or len(b) < 50:
                        continue
                    ks = ks_2samp(a, b)
                    ap = pl_alpha(a)
                    am = pl_alpha(b)
                    rows.append(dict(feature=name, cls=cls, n_plus=len(a),
                                     n_minus=len(b),
                                     median_plus=float(np.median(a)),
                                     median_minus=float(np.median(b)),
                                     ks_stat=float(ks.statistic),
                                     ks_p=float(ks.pvalue),
                                     alpha_plus=ap[0], alpha_minus=am[0]))
            T = pd.DataFrame(rows)
            T.to_csv(safe_write_path(os.path.join(OUT, "C_length_parity.tsv")),
                     sep="\t", index=False)
            st.note("length parity: " + "; ".join(
                "{} {} KS={:.3f}".format(r.feature, r.cls, r.ks_stat)
                for r in T.itertuples() if r.cls in ("PCG", "divergent",
                                                     "intergenic")))
            st.summary_md("results/strand_parity/step_summaries/46_c.md",
                          "46C length parity / power-law", df_to_md(T, ".3f"))

        elif mode == "D":
            import pysam
            fa = pysam.FastaFile("/media/share/resource/reference/"
                                 "homo_sapiens/GRCh38/genome/GRCh38.CHR.fa")
            prom = pd.read_csv(os.path.join(RES, "promoters_feat.bed"),
                               sep="\t", header=None,
                               names=["chrom", "start", "end", "gene_id",
                                      "x", "strand"])
            prom = prom.merge(meta[["cls", "group"]], left_on="gene_id",
                              right_index=True, how="inner")
            rows = []
            for r in prom.itertuples():
                try:
                    s = fa.fetch(r.chrom, r.start, r.end).upper()
                except Exception:
                    continue
                n = len(s)
                if n < 100:
                    continue
                A, T, G, C = s.count("A"), s.count("T"), s.count("G"), s.count("C")
                gc = (G - C) / max(G + C, 1)
                at = (A - T) / max(A + T, 1)
                rows.append(dict(gene_id=r.gene_id, cls=r.cls, group=r.group,
                                 n=n, gc_skew=gc, at_skew=at))
            D = pd.DataFrame(rows)
            D.to_csv(safe_write_path(os.path.join(
                OUT, "D_promoter_skew_per_gene.tsv")), sep="\t", index=False)
            tests = []
            pcg = D[D.group == "pcg"]
            for met in ["gc_skew", "at_skew"]:
                x = pcg[met].dropna().values
                for c in CLASSES:
                    y = D.loc[D.cls == c, met].dropna().values
                    if len(y) < 20:
                        continue
                    dd, p = mwu(y, x)
                    tests.append(dict(metric=met, cls=c, n=len(y),
                                      median_class=float(np.median(y)),
                                      median_PCG=float(np.median(x)),
                                      cliffs_delta=dd, pval=p))
            T = pd.DataFrame(tests)
            T["fdr"] = bh_fdr(T.pval.values)
            T.to_csv(safe_write_path(os.path.join(
                OUT, "D_skew_tests.tsv")), sep="\t", index=False)
            st.note("promoter skew: " + "; ".join(
                "{} {} d={:.3f}".format(r.metric, r.cls, r.cliffs_delta)
                for r in T.itertuples() if r.cls in ("divergent", "intergenic")))
            st.summary_md("results/strand_parity/step_summaries/46_d.md",
                          "46D promoter GC/AT skew", df_to_md(T, ".3f"))

        elif mode == "E":
            # divergent pairs: lnc cage peak vs partner PCG cage peak
            gm = meta[["cage_peak_id", "partner_pcg", "cls", "group"]].copy()
            gm = gm.join(pd.read_csv(os.path.join(RES, "cage_activity.tsv"),
                                     sep="\t", index_col=0,
                                     low_memory=False)[["max_tpm"]])
            div = gm[(gm.cls == "divergent") & gm.partner_pcg.notna()].copy()
            div["partner"] = div.partner_pcg.astype(str)
            pv = gm[["cage_peak_id", "max_tpm"]].rename(
                columns={"cage_peak_id": "pcg_peak",
                         "max_tpm": "pcg_tpm"})
            div = div.join(pv, on="partner")
            keep = set(div.cage_peak_id.dropna().astype(str)) | \
                set(div.pcg_peak.dropna().astype(str))
            # LTR genes: opposite/same peaks within TSS +/-500
            pk = []
            with gzip.open(os.path.join(DATA, "fantom",
                                        "hg38_cage_peaks.bed.gz"), "rt") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    pk.append((f[0], int(f[1]), int(f[2]), f[3], int(f[5] == "+"),
                               f[5]))
            pdf = pd.DataFrame(pk, columns=["chrom", "start", "end", "pid",
                                            "plus", "strand"])
            tel = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                              sep="\t", index_col=0, low_memory=False)[
                ["tss_in_ltr"]]
            mt = meta.join(tel)
            ltr = mt[(mt.cls == "intergenic") &
                     (mt.tss_in_ltr == 1)][["chrom", "tss", "strand"]]
            ltr_ids = {}
            for r in ltr.itertuples():
                sub = pdf[(pdf.chrom == r.chrom) &
                          (pdf.start >= r.tss - 500) &
                          (pdf.end <= r.tss + 500)]
                same = sub[sub.strand == r.strand].pid.tolist()
                opp = sub[sub.strand != r.strand].pid.tolist()
                ltr_ids[r.Index] = (same, opp)
                keep |= set(same) | set(opp)
            mean_tpm = {}
            with gzip.open(os.path.join(DATA, "fantom",
                                        "hg38_cage_tpm_ann.osc.txt.gz"),
                           "rt") as fh:
                for line in fh:
                    if line.startswith("##"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if f[0] not in keep:
                        continue
                    try:
                        v = np.array([float(x) for x in f[7:]])
                    except Exception:
                        continue
                    mean_tpm[f[0]] = float(v.mean())
            div["lnc_tpm"] = div.cage_peak_id.map(mean_tpm)
            div["partner_tpm"] = div.pcg_peak.map(mean_tpm)
            ok = div[["lnc_tpm", "partner_tpm"]].dropna()
            rho, pr = spearmanr(ok.lnc_tpm, ok.partner_tpm)
            div.to_csv(safe_write_path(os.path.join(
                OUT, "E_divergent_tpm.tsv")), sep="\t")
            # LTR asymmetry
            rows = []
            for gid, (same, opp) in ltr_ids.items():
                st_ = np.mean([mean_tpm.get(x, np.nan) for x in same]) \
                    if same else np.nan
                op_ = np.mean([mean_tpm.get(x, np.nan) for x in opp]) \
                    if opp else np.nan
                rows.append(dict(gene_id=gid, sense_tpm=st_, anti_tpm=op_,
                                 log2_ratio=np.log2((st_ + 1) / (op_ + 1))
                                 if np.isfinite(st_) and np.isfinite(op_)
                                 else np.nan))
            L = pd.DataFrame(rows)
            L.to_csv(safe_write_path(os.path.join(
                OUT, "E_ltr_asymmetry.tsv")), sep="\t", index=False)
            st.note("divergent CAGE rho={:.3f}; LTR asymmetry median="
                    "{:.3f}".format(rho, L.log2_ratio.median()))
            st.summary_md("results/strand_parity/step_summaries/46_e.md",
                          "46E transcription asymmetry",
                          "Divergent lnc/partner CAGE Spearman rho={:.3f} "
                          "(n={}).\n\nLTR sense/antisense log2 ratio median="
                          "{:.3f} (n={}).".format(rho, len(ok),
                                                  L.log2_ratio.median(),
                                                  L.log2_ratio.notna().sum()))

        elif mode == "F":
            gtfs = [("human", os.path.join(DATA, "genes", "genes.tsv")),
                    ("mouse", os.path.join(WORKDIR, "data/comparative/gtf/Mus_musculus.GRCm39.116.gtf.gz")),
                    ("rat", os.path.join(WORKDIR, "data/comparative/gtf/Rattus_norvegicus.mRatBN7.2.110.gtf.gz")),
                    ("macaque", os.path.join(WORKDIR, "data/comparative/gtf/Macaca_mulatta.Mmul_10.116.gtf.gz")),
                    ("rabbit", os.path.join(WORKDIR, "data/comparative/gtf/Oryctolagus_cuniculus.OryCun2.0.116.gtf.gz")),
                    ("chicken", os.path.join(WORKDIR, "data/comparative/gtf/Gallus_gallus.GRCg6a.104.gtf.gz")),
                    ("opossum", os.path.join(WORKDIR, "data/comparative/gtf/Monodelphis_domestica.ASM229v1.116.gtf.gz"))]
            rows = []
            for sp, path in gtfs:
                if sp == "human":
                    g = pd.read_csv(path, sep="\t")
                    genes = [(r.strand, r.end - r.start + 1)
                             for r in g.itertuples()]
                    ex = pd.read_csv(os.path.join(DATA, "genes", "exons.bed"),
                                     sep="\t", header=None,
                                     names=["c", "s", "e", "gid", "x", "strand"])
                    exons = [(r.strand, r.e - r.s) for r in ex.itertuples()]
                else:
                    genes, exons = [], []
                    with gzip.open(path, "rt") as fh:
                        for line in fh:
                            if line.startswith("#"):
                                continue
                            f = line.rstrip("\n").split("\t")
                            if len(f) < 9:
                                continue
                            if f[2] == "gene":
                                genes.append((f[6], int(f[4]) - int(f[3]) + 1))
                            elif f[2] == "exon":
                                exons.append((f[6], int(f[4]) - int(f[3]) + 1))
                for name, vals in [("gene", genes), ("exon", exons)]:
                    a = np.array([x[1] for x in vals if x[0] == "+"])
                    b = np.array([x[1] for x in vals if x[0] == "-"])
                    ks = ks_2samp(a, b)
                    ap = pl_alpha(a)
                    am = pl_alpha(b)
                    rows.append(dict(species=sp, feature=name,
                                     n_plus=len(a), n_minus=len(b),
                                     median_plus=float(np.median(a)),
                                     median_minus=float(np.median(b)),
                                     ks_stat=float(ks.statistic),
                                     ks_p=float(ks.pvalue),
                                     alpha_plus=ap[0], alpha_minus=am[0]))
            T = pd.DataFrame(rows)
            T.to_csv(safe_write_path(os.path.join(
                OUT, "F_cross_species_parity.tsv")), sep="\t", index=False)
            st.note("cross-species: " + "; ".join(
                "{} {} KS={:.3f} a+={:.2f}/a-={:.2f}".format(
                    r.species, r.feature, r.ks_stat, r.alpha_plus,
                    r.alpha_minus) for r in T.itertuples()))
            st.summary_md("results/strand_parity/step_summaries/46_f.md",
                          "46F cross-species length parity",
                          df_to_md(T, ".3f"))
        else:
            raise SystemExit("mode A-F")


if __name__ == "__main__":
    main()
