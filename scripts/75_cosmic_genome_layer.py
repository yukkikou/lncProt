#!/usr/bin/env python
"""75_cosmic_genome_layer.py — COSMIC at the locus level.

Two questions that the manuscript can answer only with a cancer-genome
resource.

  1. Gene level. Are the coding genes that lncRNAs borrow or sit inside
     enriched for Cancer Gene Census entries, and does that hold for both
     oncogenes and tumour suppressors? ONGene covers oncogenes only, so this
     replaces the previous gene list with the census.

  2. Promoter level. Do lncRNA promoters carry more or fewer recurrent somatic
     non-coding mutations than coding promoters, and does the answer depend on
     birth class? The comparison uses a flanking window as an internal control
     for local sequencing coverage, because COSMIC is dominated by targeted
     panels and the raw count mostly measures how often a locus was sequenced.

The mutation layer uses confirmed somatic non-coding variants (v104, GRCh38).
Tissue-stratified mutation analysis is not possible from the downloaded set
because the phenotype table that maps COSMIC_PHENOTYPE_ID to tumour type is
not included; this is recorded in the step summary.

Outputs
  results/cosmic/cgc_class_tests.tsv
  results/cosmic/promoter_mutation_burden.tsv
  results/cosmic/promoter_mutation_tests.tsv
  results/step_summaries/75_cosmic_genome_layer.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
COS = os.path.join(WORKDIR, "data", "cosmic", "extracted")
OUT = os.path.join(RES, "cosmic")
CGC = os.path.join(COS, "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz")
NCV = os.path.join(COS, "Cosmic_NonCodingVariants_v104_GRCh38.tsv.gz")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
PROM = 1000          # promoter half-window
FLANK_A, FLANK_B = 5000, 7000   # control window on each side
NPERM = 1000
SEED = 42


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def partial_spearman(x, y, covars):
    d = pd.DataFrame(dict(x=x, y=y, **covars)).dropna()
    if len(d) < 30:
        return np.nan, np.nan, len(d)
    rx = d.x.rank().values
    ry = d.y.rank().values
    Z = np.column_stack([d[c].rank().values for c in covars] + [np.ones(len(d))])
    bx = np.linalg.lstsq(Z, rx, rcond=None)[0]
    by = np.linalg.lstsq(Z, ry, rcond=None)[0]
    ex = rx - Z @ bx
    ey = ry - Z @ by
    r, p = spearmanr(ex, ey)
    return float(r), float(p), len(d)


def load_cgc():
    g = pd.read_csv(CGC, sep="\t", low_memory=False)
    g["ens"] = g.SYNONYMS.astype(str).str.extract(r"(ENSG\d{11})")[0]
    role = g.ROLE_IN_CANCER.fillna("").astype(str)
    g["is_onc"] = role.str.contains("oncogene", case=False)
    g["is_tsg"] = role.str.contains("TSG", case=False)
    return g


def load_noncoding():
    """Confirmed somatic non-coding positions, sorted per chromosome."""
    acc = {}
    n = 0
    for ch in pd.read_csv(NCV, sep="\t", usecols=[10, 11, 18],
                          dtype=str, chunksize=4_000_000):
        ch.columns = ["chrom", "start", "status"]
        ok = ch.status.fillna("").str.contains("somatic", case=False)
        sub = ch[ok]
        n += len(sub)
        for c, g in sub.groupby("chrom"):
            acc.setdefault(c, []).append(
                pd.to_numeric(g.start, errors="coerce").dropna()
                .astype(np.int32).values)
    return {c: np.sort(np.concatenate(v)) for c, v in acc.items()}, n


def count(intervals, chrom, a, b):
    x = intervals.get(chrom)
    if x is None or b <= a:
        return 0
    return int(np.searchsorted(x, b, "left") - np.searchsorted(x, a, "left"))


def main():
    with step("75_cosmic_genome_layer") as st:
        os.makedirs(OUT, exist_ok=True)
        rng = np.random.default_rng(SEED)

        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["base"] = m.index.astype(str).str.split(".").str[0]
        m["tss"] = pd.to_numeric(m.tss, errors="coerce")
        m["chrom_c"] = m.chrom.astype(str).str.replace("chr", "", regex=False)
        m = m.dropna(subset=["tss"])
        pcg = m[m.group == "pcg"].copy()
        lnc = m[m.group == "lnc"].copy()

        # ---------------- 1. Cancer Gene Census ------------------------------
        cgc = load_cgc()
        cgc_sym = set(cgc.GENE_SYMBOL.astype(str))
        cgc_ens = set(cgc.ens.dropna())
        pcg["cgc"] = pcg.gene_name.astype(str).isin(cgc_sym) | \
            pcg.base.isin(cgc_ens)
        onc = cgc[cgc.is_onc]
        tsg = cgc[cgc.is_tsg]
        pcg["cgc_onc"] = pcg.gene_name.astype(str).isin(set(onc.GENE_SYMBOL)) | \
            pcg.base.isin(set(onc.ens.dropna()))
        pcg["cgc_tsg"] = pcg.gene_name.astype(str).isin(set(tsg.GENE_SYMBOL)) | \
            pcg.base.isin(set(tsg.ens.dropna()))
        st.note("CGC genes={} mapped to PCG={} (oncogene {}, TSG {})".format(
            len(cgc), int(pcg.cgc.sum()), int(pcg.cgc_onc.sum()),
            int(pcg.cgc_tsg.sum())))

        rows = []

        def add(family, test, subset, stat, p, n, note):
            rows.append(dict(family=family, test=test, subset=subset,
                             statistic=stat, pval=p, n=n, note=note))

        div = lnc[lnc.cls == "divergent"].dropna(subset=["partner_pcg"])
        part = pcg.reindex(div.partner_pcg.astype(str))
        ok = part.cgc.notna()
        for col, lab in (("cgc", "CGC"), ("cgc_onc", "oncogene"),
                         ("cgc_tsg", "tumour suppressor")):
            a = int(part[col][ok].sum()); b = int(ok.sum() - a)
            c = int(pcg[col].sum()) - a
            d = int(len(pcg) - pcg[col].sum() - b)
            orr, p = fisher_exact([[a, b], [c, d]])
            add("CGC partner", "divergent partner is {}".format(lab), lab,
                orr, p, int(ok.sum()),
                "{} of {} partners, background {:.2%}".format(
                    a, int(ok.sum()), pcg[col].mean()))

        # local CGC density with a label-shuffle null (1 Mb)
        chroms = sorted(set(pcg.chrom_c) & set(lnc.chrom_c))
        pos = {g: i for i, g in enumerate(lnc.index)}
        pcg_by, query = {}, {}
        for c in chroms:
            g = pcg[pcg.chrom_c == c].sort_values("tss")
            pcg_by[c] = dict(tss=g.tss.to_numpy(float),
                             lab=g.cgc.to_numpy().astype(np.int64))
            sub = lnc[lnc.chrom_c == c]
            if len(sub):
                query[c] = (np.array([pos[x] for x in sub.index]),
                            sub.tss.to_numpy(float))

        def win_counts(win, labels=None):
            n_p = np.zeros(len(lnc), dtype=np.int64)
            n_c = np.zeros(len(lnc), dtype=np.int64)
            for c in chroms:
                if c not in query:
                    continue
                t = pcg_by[c]["tss"]
                lab = pcg_by[c]["lab"] if labels is None else \
                    np.asarray(labels[c], dtype=np.int64)
                pre = np.concatenate([[0], np.cumsum(lab)])
                ii, q = query[c]
                lo = np.searchsorted(t, q - win, "left")
                hi = np.searchsorted(t, q + win, "right")
                n_p[ii] = hi - lo
                n_c[ii] = pre[hi] - pre[lo]
            return n_p, n_c

        n_p, n_c = win_counts(1_000_000)
        lnc["n_pcg_1m"], lnc["n_cgc_1m"] = n_p, n_c
        base = {c: pcg_by[c]["lab"] for c in chroms}
        valid = n_p > 0
        for cls in CLASSES:
            sel = (lnc.cls == cls).to_numpy() & valid
            if sel.sum() < 20:
                continue
            f = n_c[sel].sum() / n_p[sel].sum()
            null = np.empty(NPERM)
            for k in range(NPERM):
                vals = rng.permutation(np.concatenate([base[c] for c in chroms]))
                perm, off = {}, 0
                for c in chroms:
                    L = len(base[c]); perm[c] = vals[off:off + L]; off += L
                _, nc = win_counts(1_000_000, perm)
                null[k] = nc[sel].sum() / n_p[sel].sum()
            p2 = min(1.0, 2 * min((1 + (null >= f).sum()) / (NPERM + 1),
                                  (1 + (null <= f).sum()) / (NPERM + 1)))
            add("CGC density", "local CGC share within 1 Mb", cls,
                float(f - null.mean()), p2, int(sel.sum()),
                "observed {:.3%} v null {:.3%}".format(f, null.mean()))

        # ---------------- 2. promoter mutation burden ------------------------
        st.note("loading non-coding variants ...")
        intervals, nvar = load_noncoding()
        st.note("confirmed somatic non-coding variants={:,}".format(nvar))

        m2 = m.copy()
        m2["n_prom"] = [count(intervals, r.chrom_c, int(r.tss) - PROM,
                              int(r.tss) + PROM) for r in m2.itertuples()]
        m2["n_flank"] = [
            count(intervals, r.chrom_c, int(r.tss) - FLANK_B,
                  int(r.tss) - FLANK_A)
            + count(intervals, r.chrom_c, int(r.tss) + FLANK_A,
                    int(r.tss) + FLANK_B) for r in m2.itertuples()]
        # Rates per bp, smoothed, then the log ratio of promoter to flank
        m2["prom_rate"] = m2.n_prom / (2.0 * PROM)
        m2["flank_rate"] = m2.n_flank / (2.0 * (FLANK_B - FLANK_A))
        m2["log2_prom_flank"] = np.log2((m2.prom_rate + 0.5) /
                                        (m2.flank_rate + 0.5))

        # covariates
        sf = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                         low_memory=False)
        sf["base"] = sf.gene_id.astype(str).str.split(".").str[0]
        sf = sf.drop_duplicates("base").set_index("base")[["gc"]]
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        ph["base"] = ph.index.astype(str).str.split(".").str[0]
        ph = ph.drop_duplicates("base").set_index("base")[
            ["phast_mean", "phylop_mean", "tss_in_ltr"]]
        gt = pd.read_csv(os.path.join(
            RES, "tissue_disease", "tables", "gtex_gene_tau.tsv"), sep="\t")
        gt = gt.rename(columns={"Unnamed: 0": "gene_id"})
        gt["base"] = gt.gene_id.astype(str).str.split(".").str[0]
        gt = gt.drop_duplicates("base").set_index("base")[
            ["tau_gtex", "max_tissue_mean"]]
        cb = pd.read_csv(os.path.join(
            RES, "cancer_bridge", "C_bridge_pergene.tsv"), sep="\t",
            index_col=0)
        cb["base"] = cb.base_id.astype(str)
        cb = cb.drop_duplicates("base").set_index("base")[["n_up", "n_de"]]
        B = m2.set_index("base").join(sf).join(ph).join(gt).join(cb)
        out_cols = ["cls", "group", "n_prom", "n_flank", "prom_rate",
                    "flank_rate", "log2_prom_flank", "gc", "phast_mean",
                    "phylop_mean", "tss_in_ltr", "tau_gtex", "n_up", "n_de"]
        B[out_cols].to_csv(safe_write_path(os.path.join(
            OUT, "promoter_mutation_burden.tsv")), sep="\t")

        # class medians and matched contrasts against coding promoters
        pcgB = B[(B.group == "pcg") & B.log2_prom_flank.notna()].copy()
        pcgB["lg"] = pcgB.log2_prom_flank
        for cls in CLASSES:
            sub = B[(B.cls == cls) & B.log2_prom_flank.notna()].copy()
            if len(sub) < 20:
                continue
            # match on GC within 0.05, nearest neighbour without replacement
            used = np.zeros(len(pcgB), bool)
            pidx = pcgB.reset_index()
            order = np.argsort(pidx.gc.values)
            gs = pidx.gc.values[order]
            mm = np.full(len(sub), -1)
            for i, (_, r) in enumerate(sub.iterrows()):
                j = np.searchsorted(gs, r.gc)
                best, bd = -1, np.inf
                for k in (j - 1, j, j + 1, j - 2, j + 2):
                    if 0 <= k < len(gs):
                        pj = order[k]
                        if used[pj]:
                            continue
                        dd = abs(gs[k] - r.gc)
                        if dd < bd and dd <= 0.05:
                            best, bd = pj, dd
                if best >= 0:
                    used[best] = True
                    mm[i] = best
            ok = mm >= 0
            if ok.sum() < 20:
                continue
            d = sub.log2_prom_flank.values[ok] - pidx.lg.values[mm[ok]]
            u, p = mannwhitneyu(sub.log2_prom_flank.values[ok],
                                pidx.lg.values[mm[ok]])
            add("promoter burden", "lncRNA v GC-matched PCG", cls,
                float(np.median(d)), p, int(ok.sum()),
                "median log2 promoter/flank difference")

        # partnered v unpartnered coding promoters
        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        partnered = set(dp.pcg_gid.astype(str).str.split(".").str[0])
        pcgB["partnered"] = pcgB.index.isin(partnered)
        pcgB = pcgB.dropna(subset=["log2_prom_flank", "gc", "max_tissue_mean"])
        a = pcgB[pcgB.partnered].log2_prom_flank
        b = pcgB[~pcgB.partnered].log2_prom_flank
        if len(a) > 50 and len(b) > 50:
            u, p = mannwhitneyu(a, b)
            add("promoter burden", "partnered v unpartnered PCG promoter",
                "PCG", float(np.median(a) - np.median(b)), p,
                int(len(a) + len(b)), "raw medians; see partial test below")
            r, pp, nn = partial_spearman(
                pcgB.log2_prom_flank, pcgB.partnered.astype(float),
                {"gc": pcgB.gc, "expr": np.log10(
                    pcgB.max_tissue_mean.clip(lower=1e-3))})
            add("promoter burden", "partial Spearman partnered",
                "PCG", r, pp, nn, "controls GC and tau (expression proxy)")

        # LTR contrast within intergenic
        ig = B[(B.cls == "intergenic") & B.log2_prom_flank.notna()
               & B.tss_in_ltr.notna()]
        al = ig[ig.tss_in_ltr == 1].log2_prom_flank
        bl = ig[ig.tss_in_ltr != 1].log2_prom_flank
        if len(al) > 20 and len(bl) > 20:
            u, p = mannwhitneyu(al, bl)
            add("promoter burden", "LTR-derived v other intergenic", "intergenic",
                float(np.median(al) - np.median(bl)), p,
                int(len(al) + len(bl)), "median log2 promoter/flank")

        # conservation and disease associations, within lncRNA and PCG
        for grp, sub in (("lncRNA", B[B.group == "lnc"]),
                         ("PCG", B[B.group == "pcg"])):
            r, p, nn = partial_spearman(
                sub.log2_prom_flank, sub.phast_mean,
                {"gc": sub.gc, "expr": np.log10(sub.max_tissue_mean.clip(lower=1e-3))})
            add("promoter burden", "partial Spearman phastCons", grp, r, p, nn,
                "controls GC and expression")
            dz = sub[(sub.n_de >= 1)]
            if len(dz) > 100:
                r, p, nn = partial_spearman(
                    dz.log2_prom_flank, dz.n_up,
                    {"gc": dz.gc, "expr": np.log10(dz.max_tissue_mean.clip(lower=1e-3))})
                add("promoter burden", "partial Spearman pan-cancer up count",
                    grp, r, p, nn, "controls GC and expression")

        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "promoter_mutation_tests.tsv")), sep="\t", index=False)
        # keep the CGC tests in their own table as well
        T[T.family.str.startswith("CGC")].to_csv(safe_write_path(
            os.path.join(OUT, "cgc_class_tests.tsv")), sep="\t", index=False)

        for r in T.itertuples():
            st.note("{} | {} | {} | stat={:.3g} n={} FDR={:.3g} | {}".format(
                r.family, r.test, r.subset, r.statistic, r.n, r.fdr, r.note))
        st.note("NOTE: tissue-stratified mutation analysis not possible, the "
                "COSMIC phenotype table is not in the downloaded set.")
        st.summary_md(os.path.join(
            RES, "step_summaries", "75_cosmic_genome_layer.md"),
            "75 COSMIC locus-level layer",
            df_to_md(T, ".4f"))


if __name__ == "__main__":
    main()
