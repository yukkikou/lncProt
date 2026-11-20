#!/usr/bin/env python
"""73_oncogene_proximity.py — do lncRNA birth classes sit near known oncogenes?

The manuscript has no cancer-gene layer. The obvious resource would be the
COSMIC Cancer Gene Census, which is not on this cluster. ONGene is available
(/media/share/resource/database/ONGene_list/ONGene.txt) and is a curated
oncogene list with Ensembl identifiers, so the same question can be asked from
a public resource without a licence.

Two questions, both local-density controlled.

  1. Divergent lncRNAs inherit their partner's promoter. Are their partners
     enriched for oncogenes relative to the 20,065 coding genes?
  2. For the other classes, is the local density of oncogenes around the
     lncRNA TSS higher than expected? The null shuffles oncogene labels among
     coding genes, which preserves both gene density and the class-specific
     distance structure, so a class that simply lives in gene-dense regions
     cannot pass.

Outputs
  results/cancer_bridge/oncogene_proximity_per_gene.tsv
  results/cancer_bridge/oncogene_proximity_tests.tsv
  results/step_summaries/73_oncogene_proximity.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "cancer_bridge")
ONGENE = "/media/share/resource/database/ONGene_list/ONGene.txt"
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
WINDOWS = [100_000, 1_000_000]
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


def ensembl_id(links):
    for t in str(links).split("|"):
        if t.startswith("Ensembl:ENSG"):
            return t.split(":")[1].split(".")[0]
    return None


def main():
    with step("73_oncogene_proximity") as st:
        on = pd.read_csv(ONGENE, sep="\t", low_memory=False)
        on["ens"] = on.Links.map(ensembl_id)
        onc = set(on.ens.dropna().astype(str))
        st.note("ONGene rows={} unique Ensembl={}".format(len(on), len(onc)))

        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["base"] = m.index.astype(str).str.split(".").str[0]
        m["tss"] = pd.to_numeric(m.tss, errors="coerce")
        m["chrom_c"] = m.chrom.astype(str).str.replace("chr", "", regex=False)
        m = m.dropna(subset=["tss"])
        m["is_onc"] = m.base.isin(onc)

        pcg = m[m.group == "pcg"].copy()
        lnc = m[m.group == "lnc"].copy()
        st.note("PCG={} oncogene PCG={} ({:.2%}); lncRNA={}".format(
            len(pcg), int(pcg.is_onc.sum()), pcg.is_onc.mean(), len(lnc)))

        # ---- per-chromosome sorted arrays, used by both the observed pass and
        # the label-shuffle null
        chroms = sorted(set(pcg.chrom_c) & set(lnc.chrom_c))
        pcg_by = {}
        for c in chroms:
            g = pcg[pcg.chrom_c == c].sort_values("tss")
            pcg_by[c] = dict(tss=g.tss.to_numpy(dtype=float),
                             onc=g.is_onc.to_numpy().astype(np.int64),
                             idx=g.index.to_numpy())
        lnc_by = {c: lnc[lnc.chrom_c == c] for c in chroms}

        # Precompute the query side once. Rebuilding these arrays inside every
        # permutation would dominate the 1,000-shuffle null.
        pos = {g: i for i, g in enumerate(lnc.index)}
        query = {}
        for c in chroms:
            sub = lnc_by[c]
            if not len(sub):
                continue
            query[c] = (np.array([pos[g] for g in sub.index]),
                        sub.tss.to_numpy(dtype=float))

        def window_counts(win, labels_by_chrom=None):
            """Per-lncRNA (#PCG, #oncogene) inside +/- win, for one label set."""
            n_pcg = np.zeros(len(lnc), dtype=np.int64)
            n_onc = np.zeros(len(lnc), dtype=np.int64)
            for c in chroms:
                if c not in query:
                    continue
                t = pcg_by[c]["tss"]
                lab = (pcg_by[c]["onc"] if labels_by_chrom is None
                       else np.asarray(labels_by_chrom[c], dtype=np.int64))
                pre = np.concatenate([[0], np.cumsum(lab)])
                ii, q = query[c]
                lo = np.searchsorted(t, q - win, side="left")
                hi = np.searchsorted(t, q + win, side="right")
                n_pcg[ii] = hi - lo
                n_onc[ii] = pre[hi] - pre[lo]
            return n_pcg, n_onc

        # observed
        obs = {}
        for w in WINDOWS:
            obs[w] = window_counts(w)
        for j, w in enumerate(WINDOWS):
            lnc["n_pcg_{}k".format(w // 1000)] = obs[w][0]
            lnc["n_onc_{}k".format(w // 1000)] = obs[w][1]

        # nearest PCG and its oncogene status
        near_d, near_onc = np.full(len(lnc), np.nan), np.zeros(len(lnc), bool)
        pos = {g: i for i, g in enumerate(lnc.index)}
        for c in chroms:
            t = pcg_by[c]["tss"]
            o = pcg_by[c]["onc"].astype(bool)
            sub = lnc_by[c]
            if not len(sub):
                continue
            q = sub.tss.to_numpy(dtype=float)
            j = np.clip(np.searchsorted(t, q), 1, len(t) - 1)
            left, right = j - 1, j
            dl = np.abs(q - t[left]); dr = np.abs(q - t[right])
            take = np.where(dl <= dr, left, right)
            ii = np.array([pos[g] for g in sub.index])
            near_d[ii] = np.minimum(dl, dr)
            near_onc[ii] = o[take]
        lnc["nearest_pcg_dist"] = near_d
        lnc["nearest_pcg_oncogene"] = near_onc

        rows = []

        def add(test, subset, stat, p, n, note):
            rows.append(dict(test=test, subset=subset, statistic=stat,
                             pval=p, n=n, note=note))

        # ---- 1. divergent partners -----------------------------------------
        div = lnc[lnc.cls == "divergent"].dropna(subset=["partner_pcg"])
        partner = m.reindex(div.partner_pcg.astype(str))
        ok = partner.is_onc.notna()
        a = int(partner.is_onc[ok].sum())
        b = int(ok.sum() - a)
        c = int(pcg.is_onc.sum()) - a
        d = int(len(pcg) - pcg.is_onc.sum() - b)
        orr, p = fisher_exact([[a, b], [c, d]])
        add("divergent partner", "partner PCG is an ONGene", orr, p, int(ok.sum()),
            "{} of {} partners, background {:.2%}".format(
                a, int(ok.sum()), pcg.is_onc.mean()))

        # ---- 2. local oncogene density, label-shuffle null ------------------
        rng = np.random.default_rng(SEED)
        base_labels = {c: pcg_by[c]["onc"] for c in chroms}
        for w in WINDOWS:
            n_pcg, n_onc = obs[w]
            valid = n_pcg > 0
            tags = ["{}k".format(w // 1000)]
            for cls in CLASSES:
                sel = (lnc.cls == cls).to_numpy() & valid
                if sel.sum() < 20:
                    continue
                o = n_onc[sel].sum(); n = n_pcg[sel].sum()
                f = o / n
                # null: shuffle oncogene labels among coding genes, 1,000x
                null = np.empty(NPERM)
                for k in range(NPERM):
                    perm = {}
                    all_idx = np.concatenate([np.arange(len(base_labels[c]))
                                              for c in chroms])
                    vals = np.concatenate([base_labels[c] for c in chroms])
                    vals = rng.permutation(vals)
                    off = 0
                    for c in chroms:
                        L = len(base_labels[c])
                        perm[c] = vals[off:off + L]; off += L
                    _, n_onc_p = window_counts(w, perm)
                    null[k] = n_onc_p[sel].sum() / n
                p_ge = (1 + (null >= f).sum()) / (NPERM + 1)
                p_le = (1 + (null <= f).sum()) / (NPERM + 1)
                p2 = min(1.0, 2 * min(p_ge, p_le))
                add("local density {}k".format(w // 1000), cls,
                    float(f - null.mean()), p2, int(sel.sum()),
                    "observed {:.3%} v null {:.3%}".format(f, null.mean()))

        # ---- 3. LTR-derived intergenic -------------------------------------
        ig = lnc[(lnc.cls == "intergenic")].copy()
        if "tss_in_ltr" not in ig.columns:
            te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                             sep="\t", index_col=0, low_memory=False)
            ig = ig.join(te[["tss_in_ltr"]])
        a = ig[ig.tss_in_ltr == 1]; b = ig[ig.tss_in_ltr != 1]
        if len(a) > 20 and len(b) > 20:
            fa = a.n_onc_1000k.sum() / max(a.n_pcg_1000k.sum(), 1)
            fb = b.n_onc_1000k.sum() / max(b.n_pcg_1000k.sum(), 1)
            u, p = mannwhitneyu(a.n_onc_1000k / a.n_pcg_1000k.clip(lower=1),
                                b.n_onc_1000k / b.n_pcg_1000k.clip(lower=1))
            add("LTR", "intergenic LTR-derived v other", float(fa - fb), p,
                int(len(a) + len(b)),
                "pooled oncogene share within 1 Mb {:.3%} v {:.3%}".format(fa, fb))

        # ---- 4. does oncogene proximity track tumour up-regulation? --------
        cb = pd.read_csv(os.path.join(OUT, "C_bridge_pergene.tsv"), sep="\t",
                         index_col=0)
        cb["base"] = cb.base_id.astype(str)
        cb = cb.drop_duplicates("base").set_index("base")
        lnc2 = lnc.copy()
        lnc2["base"] = lnc2.index.astype(str).str.split(".").str[0]
        j = lnc2.set_index("base").join(
            cb[["n_up", "n_de", "breadth"]], how="inner")
        j = j[(j.n_de >= 1) & (j.n_pcg_1000k > 0)]
        if len(j) > 100:
            f = j.n_onc_1000k / j.n_pcg_1000k
            r, p = spearmanr(f, j.n_up, nan_policy="omit")
            add("tumour", "oncogene density v pan-cancer up count", r, p,
                int(len(j)), "Spearman, lncRNAs with a DE call")
            j2 = j.dropna(subset=["breadth"])
            if len(j2) > 100:
                from numpy import polyfit, polyval
                x1 = f.loc[j2.index].rank().values
                x2 = pd.Series(j2.n_up.values).rank().values
                x3 = j2.breadth.rank().values
                e1 = x1 - polyval(polyfit(x3, x1, 1), x3)
                e2 = x2 - polyval(polyfit(x3, x2, 1), x3)
                r2, p2 = spearmanr(e1, e2)
                add("tumour", "partial, controlling expression breadth", r2, p2,
                    int(len(j2)), "rank residuals on breadth")

        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "oncogene_proximity_tests.tsv")), sep="\t", index=False)
        keep = ["cls", "chrom_c", "tss", "is_onc", "partner_pcg",
                "n_pcg_100k", "n_onc_100k", "n_pcg_1000k", "n_onc_1000k",
                "nearest_pcg_dist", "nearest_pcg_oncogene"]
        lnc[keep].to_csv(safe_write_path(os.path.join(
            OUT, "oncogene_proximity_per_gene.tsv")), sep="\t")

        for r in T.itertuples():
            st.note("{} | {} | stat={:.3g} n={} FDR={:.3g} | {}".format(
                r.test, r.subset, r.statistic, r.n, r.fdr, r.note))
        st.summary_md(os.path.join(
            RES, "step_summaries", "73_oncogene_proximity.md"),
            "73 lncRNA classes and oncogene proximity",
            df_to_md(T, ".4f"))


if __name__ == "__main__":
    main()
