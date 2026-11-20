#!/usr/bin/env python
"""76_cosmic_tissue_matching.py — does the tumour tissue of promoter mutations
match the normal tissue where the lncRNA is expressed?

The previous COSMIC pass (module 75) rejected promoter mutation density because
it tracks sequencing coverage. This module asks a different question that does
not depend on absolute density. For every gene, it records which primary
tumour sites contribute non-coding somatic mutations inside the promoter, then
asks whether the dominant tumour site is the tissue in which the gene is most
highly expressed in GTEx.

The comparison is restricted to primary sites that map cleanly onto a GTEx
tissue group (for example lung, liver, breast, colon, brain, skin, blood).
Ambiguous categories such as soft_tissue, bone and upper_aerodigestive_tract
are excluded. Concordance is tested against 1,000 permutations of the dominant
tumour site across genes, so tumour-type frequencies are preserved.

Outputs
  results/cosmic/promoter_site_spectrum.tsv
  results/cosmic/tissue_matching_tests.tsv
  results/step_summaries/76_cosmic_tissue_matching.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
COS = os.path.join(WORKDIR, "data", "cosmic", "extracted")
OUT = os.path.join(RES, "cosmic")
NCV = os.path.join(COS, "Cosmic_NonCodingVariants_v104_GRCh38.tsv.gz")
CLS = os.path.join(COS, "Cosmic_Classification_v104_GRCh38.tsv.gz")
TOIL = os.path.join(WORKDIR, "data", "tissue_disease", "toil_class_expr.npz")
TSAM = os.path.join(WORKDIR, "data", "tissue_disease", "toil_samples.tsv")
TGEN = os.path.join(WORKDIR, "data", "tissue_disease", "toil_genes.tsv")
PROM = 1000
NPERM = 1000
SEED = 42

# COSMIC primary site -> shared tissue group
SITE_GROUP = {
    "lung": "lung", "liver": "liver", "breast": "breast",
    "prostate": "prostate", "large_intestine": "colon",
    "kidney": "kidney", "thyroid": "thyroid", "stomach": "stomach",
    "pancreas": "pancreas", "testis": "testis", "ovary": "ovary",
    "endometrium": "uterus", "cervix": "cervix",
    "fallopian_tube": "fallopian", "vagina": "vagina",
    "haematopoietic_and_lymphoid_tissue": "blood",
    "central_nervous_system": "brain", "skin": "skin",
    "oesophagus": "esophagus", "salivary_gland": "salivary",
    "pituitary": "pituitary", "adrenal_gland": "adrenal",
    "small_intestine": "small_intestine",
}

# GTEx tissue label -> shared tissue group
GTEX_GROUP = {
    "Lung": "lung", "Liver": "liver", "Breast - Mammary Tissue": "breast",
    "Prostate": "prostate", "Colon - Transverse": "colon",
    "Colon - Sigmoid": "colon", "Kidney - Cortex": "kidney",
    "Thyroid": "thyroid", "Stomach": "stomach", "Pancreas": "pancreas",
    "Testis": "testis", "Ovary": "ovary", "Uterus": "uterus",
    "Cervix - Ectocervix": "cervix", "Cervix - Endocervix": "cervix",
    "Fallopian Tube": "fallopian", "Vagina": "vagina",
    "Whole Blood": "blood", "Cells - EBV-transformed lymphocytes": "blood",
    "Spleen": "blood", "Skin - Sun Exposed (Lower leg)": "skin",
    "Skin - Not Sun Exposed (Suprapubic)": "skin",
    "Esophagus - Mucosa": "esophagus", "Esophagus - Muscularis": "esophagus",
    "Esophagus - Gastroesophageal Junction": "esophagus",
    "Minor Salivary Gland": "salivary", "Pituitary": "pituitary",
    "Adrenal Gland": "adrenal",
    "Small Intestine - Terminal Ileum": "small_intestine",
}


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def rowwise_spearman(A, B):
    """Spearman correlation between matching rows of two matrices."""
    Ar = np.apply_along_axis(rankdata, 1, A)
    Br = np.apply_along_axis(rankdata, 1, B)
    Az = Ar - Ar.mean(axis=1, keepdims=True)
    Bz = Br - Br.mean(axis=1, keepdims=True)
    den = np.sqrt((Az * Az).sum(axis=1) * (Bz * Bz).sum(axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (Az * Bz).sum(axis=1) / den
    r[den == 0] = np.nan
    return r


def main():
    with step("76_cosmic_tissue_matching") as st:
        os.makedirs(OUT, exist_ok=True)

        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["base"] = m.index.astype(str).str.split(".").str[0]
        m["tss"] = pd.to_numeric(m.tss, errors="coerce")
        m["chrom_c"] = m.chrom.astype(str).str.replace("chr", "", regex=False)
        m = m.dropna(subset=["tss"])
        genes = m.reset_index(drop=True)
        gi = {g: i for i, g in enumerate(genes.base)}

        # promoter index per chromosome
        by_chr = {}
        for c, g in genes.groupby("chrom_c"):
            o = np.argsort(g.tss.values)
            by_chr[c] = (g.tss.values[o], g.index.values[o])

        # classification
        cl = pd.read_csv(CLS, sep="\t")
        cl["group"] = cl.PRIMARY_SITE.map(SITE_GROUP)
        cl = cl.dropna(subset=["group"])
        pheno = dict(zip(cl.COSMIC_PHENOTYPE_ID.astype(str), cl.group))
        groups = sorted(set(SITE_GROUP.values()))
        gidx = {g: i for i, g in enumerate(groups)}
        st.note("classification rows={} mapped sites={} tissue groups={}".format(
            len(cl), len(pheno), len(groups)))

        SPT = os.path.join(OUT, "promoter_site_spectrum.tsv")
        if "--reuse-spectrum" in sys.argv and os.path.exists(SPT):
            spec = pd.read_csv(SPT, sep="\t")
            st.note("reusing promoter spectrum: {}".format(SPT))
        else:
            counts = np.zeros((len(genes), len(groups)), dtype=np.int32)
            nvar = 0
            used = 0
            for ch in pd.read_csv(NCV, sep="\t", usecols=[6, 10, 11, 18],
                                  dtype=str, chunksize=4_000_000):
                ch.columns = ["pheno", "chrom", "start", "status"]
                ch = ch[ch.status.fillna("").str.contains("somatic", case=False)]
                nvar += len(ch)
                ch["site"] = ch.pheno.map(pheno)
                ch = ch.dropna(subset=["site"])
                ch["start"] = pd.to_numeric(ch.start, errors="coerce")
                ch = ch.dropna(subset=["start"])
                for c, sub in ch.groupby("chrom"):
                    if c not in by_chr:
                        continue
                    tss, gids = by_chr[c]
                    pos = sub.start.values.astype(np.int64)
                    lo = np.searchsorted(tss, pos - PROM, "left")
                    hi = np.searchsorted(tss, pos + PROM, "right")
                    h = np.where(hi > lo)[0]
                    if not len(h):
                        continue
                    si = sub.site.map(gidx).values
                    for k in h:
                        for j in range(lo[k], hi[k]):
                            counts[gids[j], si[k]] += 1
                            used += 1
            st.note("confirmed somatic variants={:,}; promoter assignments={:,}"
                    .format(nvar, used))
            spec = pd.DataFrame(counts, columns=groups)
            spec.insert(0, "gene", genes.base.values)
            tot = spec[groups].sum(axis=1).values
            with np.errstate(divide="ignore", invalid="ignore"):
                p = spec[groups].values / tot[:, None]
                ent = -(np.where(p > 0, p * np.log(p), 0)).sum(axis=1)
            spec["n_site_mut"] = tot
            spec["n_sites"] = (spec[groups].values > 0).sum(axis=1)
            spec["entropy"] = ent
            spec["dominant_site"] = np.where(
                tot > 0,
                np.array(groups)[np.argmax(spec[groups].values, axis=1)], "")
            spec.to_csv(safe_write_path(SPT), sep="\t", index=False)
        spec = spec.rename(columns={"gene": "base"})

        # ---------------- GTEx tissue means ---------------------------------
        z = np.load(TOIL, mmap_mode="r")
        M = np.asarray(z["M"])
        samples = pd.read_csv(TSAM, sep="\t")
        tgenes = pd.read_csv(TGEN, sep="\t")
        tgenes["base"] = tgenes.gene_id.astype(str).str.split(".").str[0]
        tgenes = tgenes.drop_duplicates("base")
        gexpr = pd.DataFrame({"base": tgenes.base.values})
        gcols = samples[(samples.cohort == "GTEx")
                        & samples.tissue.notna()].index.values
        st.note("GTEx samples with tissue={}".format(len(gcols)))
        tissues = samples.tissue.values[gcols]
        gmap = pd.Series(tissues).map(GTEX_GROUP)
        keep = gmap.notna().values
        tcols = gcols[keep]
        tgroups = gmap.values[keep]
        means = {}
        for t in sorted(set(tgroups)):
            idx = tcols[tgroups == t]
            means[t] = M[:, idx].mean(axis=1)
        E = pd.DataFrame(means)
        E.insert(0, "base", tgenes.base.values)
        Efull = E.copy()
        Efull.to_csv(safe_write_path(os.path.join(
            OUT, "gtex_tissue_means_by_group.tsv")), sep="\t", index=False)
        # top and top-3 groups per gene
        arr = E[list(E.columns[1:])].values
        order = np.argsort(arr, axis=1)[:, ::-1]
        names = np.array(E.columns[1:])
        top1 = names[order[:, 0]]
        top3 = [set(names[order[i, :3]]) for i in range(len(E))]
        E["gtex_top"] = top1
        E["gtex_top3"] = top3
        E = E[["base", "gtex_top", "gtex_top3"]]
        del M, z

        S = spec.merge(E, on="base", how="left")
        S = S.merge(genes[["base", "cls", "group"]], on="base", how="left")
        S = S[~S.cls.isna()]
        S["dominant_ok"] = S.dominant_site.isin(list(GTEX_GROUP.values()))
        S["concordant"] = (S.dominant_ok
                           & (S.dominant_site == S.gtex_top))
        S["in_top3"] = [d in t3 if isinstance(t3, set) else False
                        for d, t3 in zip(S.dominant_site, S.gtex_top3)]
        S.drop(columns=["gtex_top3"]).to_csv(safe_write_path(os.path.join(
            OUT, "promoter_site_spectrum_annotated.tsv")), sep="\t",
            index=False)

        rows = []

        def add(test, subset, stat, p, n, note):
            rows.append(dict(test=test, subset=subset, statistic=stat,
                             pval=p, n=n, note=note))

        # ---------------- mutation spectrum vs expression profile ---------
        mutcols = [c for c in groups if c in Efull.columns]
        Em = Efull.set_index("base").reindex(spec.base)[mutcols].to_numpy(float)
        Mm = spec[mutcols].to_numpy(float)
        gser = spec.base.map(genes.set_index("base").group).values
        elig = ((Mm.sum(1) >= 5) & (Em.sum(1) > 0)
                & (Mm.std(1) > 0) & (Em.std(1) > 0))
        obs_r = rowwise_spearman(Mm[elig], Em[elig])
        rng = np.random.default_rng(SEED)
        K = 100
        null_all = np.empty(K)
        null_lnc = np.empty(K)
        null_pcg = np.empty(K)
        is_lnc = (gser[elig] == "lnc")
        is_pcg = (gser[elig] == "pcg")
        for k in range(K):
            perm = rng.permutation(len(mutcols))
            rp = rowwise_spearman(Mm[elig], Em[elig][:, perm])
            null_all[k] = np.nanmean(rp)
            null_lnc[k] = np.nanmean(rp[is_lnc])
            null_pcg[k] = np.nanmean(rp[is_pcg])
        for lab, mask, null in (
                ("all genes", np.ones(len(obs_r), bool), null_all),
                ("lncRNA", is_lnc, null_lnc), ("PCG", is_pcg, null_pcg)):
            o = np.nanmean(obs_r[mask])
            p = (1 + (null >= o).sum()) / (K + 1)
            add("profile match",
                "mutation spectrum vs GTEx profile, {}".format(lab),
                float(o - null.mean()), p, int(mask.sum()),
                "mean row-wise Spearman {:.3f} v null {:.3f}".format(
                    o, null.mean()))


        # concordance in genes with enough mapped mutations
        for minm in (3, 5):
            d = S[S.dominant_ok & (S.n_site_mut >= minm)].copy()
            if len(d) < 50:
                continue
            obs = d.concordant.mean()
            obs3 = d.in_top3.mean()
            rng = np.random.default_rng(SEED)
            null = np.empty(NPERM); null3 = np.empty(NPERM)
            dom = d.dominant_site.values
            top1v = d.gtex_top.values
            top3v = d.gtex_top3.values
            for k in range(NPERM):
                perm = rng.permutation(dom)
                null[k] = (perm == top1v).mean()
                null3[k] = np.mean([(p in t3) if isinstance(t3, set) else False
                                    for p, t3 in zip(perm, top3v)])
            p1 = (1 + (null >= obs).sum()) / (NPERM + 1)
            p3 = (1 + (null3 >= obs3).sum()) / (NPERM + 1)
            add("tissue match",
                "dominant site = top GTEx, all genes >= {} mut".format(minm),
                float(obs - null.mean()), p1, int(len(d)),
                "observed {:.1%} v null {:.1%}".format(obs, null.mean()))
            add("tissue match",
                "dominant site in top-3 GTEx, all genes >= {} mut".format(minm),
                float(obs3 - null3.mean()), p3, int(len(d)),
                "observed {:.1%} v null {:.1%}".format(obs3, null3.mean()))

            # by group, lncRNAs only
            ln = d[d.group == "lnc"]
            if len(ln) > 50:
                obs_l = ln.concordant.mean()
                null_l = np.empty(NPERM)
                for k in range(NPERM):
                    null_l[k] = (rng.permutation(ln.dominant_site.values)
                                 == ln.gtex_top.values).mean()
                add("tissue match",
                    "dominant site = top GTEx, lncRNA >= {} mut".format(minm),
                    float(obs_l - null_l.mean()),
                    (1 + (null_l >= obs_l).sum()) / (NPERM + 1), int(len(ln)),
                    "observed {:.1%} v null {:.1%}".format(
                        obs_l, null_l.mean()))

        # entropy vs GTEx tau
        gt = pd.read_csv(os.path.join(
            RES, "tissue_disease", "tables", "gtex_gene_tau.tsv"), sep="\t")
        gt = gt.rename(columns={"Unnamed: 0": "base"})
        gt["base"] = gt.base.astype(str).str.split(".").str[0]
        gt = gt.drop_duplicates("base").set_index("base")
        S2 = S.set_index("base").join(gt[["tau_gtex", "max_tissue_mean"]])
        for grp in ("lnc", "pcg"):
            d = S2[(S2.group == grp) & (S2.n_site_mut >= 3)]
            d = d.dropna(subset=["tau_gtex"])
            if len(d) < 50:
                continue
            r, p = spearmanr(d.entropy, d.tau_gtex, nan_policy="omit")
            add("entropy",
                "Spearman entropy vs GTEx tau, {}".format(grp),
                r, p, int(len(d)), "genes with >= 3 mapped mutations")

        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "tissue_matching_tests.tsv")), sep="\t", index=False)
        for r in T.itertuples():
            st.note("{} | {} | {} | stat={:.3f} n={} FDR={:.3g} | {}".format(
                r.test, r.subset, r.test, r.statistic, r.n, r.fdr, r.note))
        st.summary_md(os.path.join(
            RES, "step_summaries", "76_cosmic_tissue_matching.md"),
            "76 COSMIC tumour site against GTEx tissue",
            df_to_md(T, ".4f"))


if __name__ == "__main__":
    main()
