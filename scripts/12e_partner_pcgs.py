#!/usr/bin/env python
"""12e_partner_pcgs.py — what makes a PCG likely to host a divergent lncRNA?

Compare divergent-partner PCGs vs all other PCGs (activity, conservation,
GC, chromatin) + GO enrichment via g:Profiler REST (network; degrade
gracefully if unreachable).
Outputs: results/partner_pcg_traits.tsv, results/partner_pcg_go.tsv
"""
import json
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")


def gprofiler(query_genes, organism="hsapiens"):
    """POST gene symbols to g:Profiler; return tidy DataFrame or None."""
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    payload = {
        "organism": organism,
        "query": query_genes,
        "sources": ["GO:BP", "GO:MF"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "no_evidences": True,
    }
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "lncProt/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        res = json.load(r)
    rows = [dict(term=f"{r_['source']}:{r_['native']}", name=r_["name"],
                 pval=r_["p_value"], overlap=r_["query_size"],
                 term_size=r_["term_size"], effective=r_["effective_domain_size"])
            for r_ in res["result"]]
    return pd.DataFrame(rows).sort_values("pval").head(25) if rows else None


def main():
    with step("12e_partner_pcgs") as st:
        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        seq = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["gc", "cpg_oe"]]
        cons = pd.read_csv(os.path.join(RES, "conservation.tsv"), sep="\t",
                           index_col=0, low_memory=False)[["phast_mean"]]
        act = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                          index_col=0, low_memory=False)[["active", "breadth"]]
        m = m.join(seq).join(cons).join(act)

        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        partners = set(dp["pcg_gid"])
        pcg = m[m.group == "pcg"].copy()
        pcg["is_partner"] = pcg.index.isin(partners)
        st.note(f"unique partner PCGs={pcg.is_partner.sum()}")

        rows = []
        for feat in ["gc", "cpg_oe", "phast_mean", "active", "breadth"]:
            a = pcg.loc[pcg.is_partner, feat].dropna()
            b = pcg.loc[~pcg.is_partner, feat].dropna()
            rows.append(dict(feature=feat, n_partner=len(a), n_other=len(b),
                             mean_partner=float(a.mean()),
                             mean_other=float(b.mean()),
                             ratio=float(a.mean() / b.mean()) if b.mean() else np.nan))
        t = pd.DataFrame(rows)
        t.to_csv(safe_write_path(os.path.join(RES, "partner_pcg_traits.tsv")),
                 sep="\t", index=False)
        st.note("partner traits: " + "; ".join(
            f"{r.feature}: {r.mean_partner:.3f} vs {r.mean_other:.3f}"
            for r in t.itertuples()))

        # GO enrichment (top-level: functions enriched among partner PCGs)
        symbols = pcg.loc[pcg.is_partner, "gene_name"].dropna().tolist()
        go = None
        try:
            go = gprofiler(symbols)
        except Exception as e:
            st.note(f"g:Profiler unreachable ({type(e).__name__}) — GO skipped")
        if go is not None and len(go):
            go.to_csv(safe_write_path(os.path.join(RES, "partner_pcg_go.tsv")),
                      sep="\t", index=False)
            st.note("top GO: " + "; ".join(
                f"{r.name} (p={r.pval:.1e})" for r in go.head(5).itertuples()))
        else:
            st.note("no GO result")

        body = ("## Partner-PCG vs other PCG traits\n\n" + df_to_md(t, ".3f")
                + ("\n\n## Partner-PCG GO enrichment (g:Profiler, SCS)\n\n"
                   + df_to_md(go, ".2g", max_rows=25) if go is not None
                   else "\n\n_(GO unavailable)_"))
        st.summary_md("results/step_summaries/12e_partner_pcgs.md",
                      "12e partner PCG traits", body)


if __name__ == "__main__":
    main()
