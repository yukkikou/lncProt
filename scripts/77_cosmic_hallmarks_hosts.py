#!/usr/bin/env python
"""77_cosmic_hallmarks_hosts.py — which cancer processes do the borrowed and
host promoters serve?

Module 75 showed that divergent partners are enriched for Cancer Gene Census
genes, with the signal at tumour suppressors. Two natural follow-ups use the
same census without any new data.

  1. Host genes. Intronic and antisense lncRNAs sit inside a coding gene. Is
     that host enriched for census genes, and does the same tumour-suppressor
     bias appear?
  2. Cancer hallmarks. Among census genes, are the ones chosen as partners or
     hosts concentrated in particular hallmark categories, or do they span the
     census uniformly? The background is restricted to census genes, so the
     answer is not a restatement of census membership.
  3. Fusion role and age. Are fusion-role genes over-represented among
     partners, and are pairs at census genes older than pairs at other coding
     genes?

Outputs
  results/cosmic/cgc_hallmark_host_tests.tsv
  results/step_summaries/77_cosmic_hallmarks_hosts.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
COS = os.path.join(WORKDIR, "data", "cosmic", "extracted")
OUT = os.path.join(RES, "cosmic")
CGC = os.path.join(COS, "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz")
HALL = os.path.join(COS, "Cosmic_CancerGeneCensusHallmarksOfCancer_v104_GRCh38.tsv.gz")

HALLMARKS = [
    "proliferative signalling", "escaping programmed cell death",
    "genome instability and mutations", "cell division control",
    "angiogenesis", "change of cellular energetics",
    "escaping immune response to cancer",
    "tumour promoting inflammation", "cell replicative immortality",
    "invasion and metastasis", "suppression of growth", "senescence",
    "differentiation and development", "global regulation of gene expression",
]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def main():
    with step("77_cosmic_hallmarks_hosts") as st:
        g = pd.read_csv(CGC, sep="\t", low_memory=False)
        g["ens"] = g.SYNONYMS.astype(str).str.extract(r"(ENSG\d{11})")[0]
        role = g.ROLE_IN_CANCER.fillna("").astype(str)
        g["is_onc"] = role.str.contains("oncogene", case=False)
        g["is_tsg"] = role.str.contains("TSG", case=False)
        g["is_fusion"] = role.str.contains("fusion", case=False)

        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["base"] = m.index.astype(str).str.split(".").str[0]
        pcg = m[m.group == "pcg"].copy()

        def flag(df, col):
            sy = set(g.loc[g[col], "GENE_SYMBOL"].astype(str))
            en = set(g.loc[g[col], "ens"].dropna())
            return df.gene_name.astype(str).isin(sy) | df.base.isin(en)

        for col, lab in (("is_onc", "onc"), ("is_tsg", "tsg"),
                         ("is_fusion", "fusion")):
            pcg[lab] = flag(pcg, col)
        # census membership from the full table, independent of role string
        all_sy = set(g.GENE_SYMBOL.astype(str))
        all_en = set(g.ens.dropna())
        pcg["cgc"] = pcg.gene_name.astype(str).isin(all_sy) | pcg.base.isin(all_en)
        st.note("PCG census={} onc={} tsg={} fusion={}".format(
            int(pcg.cgc.sum()), int(pcg.onc.sum()), int(pcg.tsg.sum()),
            int(pcg.fusion.sum())))

        rows = []

        def add(family, test, subset, stat, p, n, note):
            rows.append(dict(family=family, test=test, subset=subset,
                             statistic=stat, pval=p, n=n, note=note))

        # ---------------- host genes of intronic and antisense --------------
        for cls in ("intronic", "antisense"):
            sub = m[(m.cls == cls) & m.partner_pcg.notna()].copy()
            # one row per distinct host gene. Many lncRNAs share a host, so a
            # per-lncRNA test would treat the same gene as independent evidence
            # several times over.
            hosts = sorted(set(
                sub.partner_pcg.astype(str).str.split(".").str[0]))
            host = pcg.set_index("base").reindex(hosts)
            ok = host.cgc.notna().values
            for col, lab in (("cgc", "CGC"), ("onc", "oncogene"),
                             ("tsg", "tumour suppressor"),
                             ("fusion", "fusion")):
                a = int(host[col].values[ok].sum()); b = int(ok.sum() - a)
                c = int(pcg[col].sum()) - a
                d = int(len(pcg) - pcg[col].sum() - b)
                orr, p = fisher_exact([[a, b], [c, d]])
                add("host gene", "host is {}".format(lab), cls, orr, p,
                    int(ok.sum()),
                    "{} of {} distinct hosts, background {:.2%}".format(
                        a, int(ok.sum()), pcg[col].mean()))

        # ---------------- hallmark specificity among census genes ----------
        h = pd.read_csv(HALL, sep="\t", low_memory=False)
        h = h[h.HALLMARK.isin(HALLMARKS)]
        hp = h[["GENE_SYMBOL", "HALLMARK"]].drop_duplicates()
        st.note("hallmark pairs={} genes={}".format(len(hp),
                                                    hp.GENE_SYMBOL.nunique()))
        cgc_pcg = pcg[pcg.cgc].copy()
        partner_bases = set(
            m[m.cls == "divergent"].partner_pcg.dropna().astype(str)
            .str.split(".").str[0])
        cgc_pcg["partner"] = cgc_pcg.base.isin(partner_bases)
        hosts_all = set()
        for cls in ("intronic", "antisense"):
            hosts_all |= set(m[(m.cls == cls)].partner_pcg.dropna()
                             .astype(str).str.split(".").str[0])
        cgc_pcg["host"] = cgc_pcg.base.isin(hosts_all)
        for hname in HALLMARKS:
            genes = set(hp.loc[hp.HALLMARK == hname, "GENE_SYMBOL"])
            in_h = cgc_pcg.gene_name.astype(str).isin(genes)
            for sub, lab in ((cgc_pcg.partner, "divergent partner"),
                             (cgc_pcg.host, "intronic/antisense host")):
                a = int((in_h & sub).sum()); b = int((~in_h & sub).sum())
                c = int((in_h & ~sub).sum()); d = int((~in_h & ~sub).sum())
                if a + b == 0:
                    continue
                orr, p = fisher_exact([[a, b], [c, d]])
                add("hallmark", hname, lab, orr, p, int(a + b),
                    "{} of {} census {} genes".format(a, a + b, lab))

        # ---------------- age and conservation of census pairs -------------
        ev = pd.read_csv(os.path.join(
            RES, "comparative", "E_transcript_age.tsv"), sep="\t",
            index_col=0)
        ev["base"] = ev.index.astype(str).str.split(".").str[0]
        ev = ev.drop_duplicates("base").set_index("base")[["age_level"]]
        div = m[m.cls == "divergent"].dropna(subset=["partner_pcg"]).copy()
        div["pbase"] = div.partner_pcg.astype(str).str.split(".").str[0]
        div["cgc_p"] = div.pbase.map(pcg.set_index("base").cgc).fillna(False)
        div["tsg_p"] = div.pbase.map(pcg.set_index("base").tsg).fillna(False)
        div = div.join(ev, on="base")
        for flag, lab in (("cgc_p", "CGC"), ("tsg_p", "TSG")):
            a = div[div[flag]].age_level.dropna()
            b = div[~div[flag]].age_level.dropna()
            if len(a) > 20 and len(b) > 20:
                u, p = mannwhitneyu(a, b)
                add("age", "transcript age, {} partner vs other".format(lab),
                    "divergent", float(a.median() - b.median()), p,
                    int(len(a) + len(b)),
                    "median age level {:.0f} v {:.0f}".format(
                        a.median(), b.median()))

        T = pd.DataFrame(rows)
        for fam, sub in T.groupby("family"):
            T.loc[sub.index, "fdr"] = bh_fdr(sub.pval.values)
        T = T.sort_values(["family", "pval"])
        T.to_csv(safe_write_path(os.path.join(
            OUT, "cgc_hallmark_host_tests.tsv")), sep="\t", index=False)

        for r in T.itertuples():
            st.note("{} | {} | {} | stat={:.3g} n={} FDR={:.3g}".format(
                r.family, r.test, r.subset, r.statistic, r.n, r.fdr))
        st.summary_md(os.path.join(
            RES, "step_summaries", "77_cosmic_hallmarks_hosts.md"),
            "77 census hallmarks and host genes",
            df_to_md(T, ".4f"))


if __name__ == "__main__":
    main()
