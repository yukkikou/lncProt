#!/usr/bin/env python
"""68_tumour_breadth_adjusted.py — does class still predict tumour up-regulation
once expression breadth is accounted for?

The Results report the fraction of genes up-regulated per cancer as an
unadjusted comparison. Within lncRNAs, that fraction rises steeply with
expression breadth, and the classes differ in breadth by more than an order of
magnitude, so the ranking among classes may be a breadth ranking. This module
fits a logistic model with breadth as a covariate and reports the adjusted
difference between classes, plus a directly matched comparison.

Reads
  results/cancer_bridge/C_bridge_pergene.tsv   per-gene pan-cancer DE counts,
                                               breadth and class

Outputs
  results/cancer_bridge/D_breadth_adjusted.tsv
  results/step_summaries/68_tumour_breadth.md
"""
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "cancer_bridge")
CLS = ["intergenic", "intronic", "antisense", "divergent", "sense_overlapping"]


def main():
    with step("68_tumour_breadth") as st:
        d = pd.read_csv(os.path.join(OUT, "C_bridge_pergene.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        d = d[d.cls.isin(CLS)].dropna(subset=["breadth"]).copy()
        # Genes never tested in a cancer cannot be counted as not up-regulated.
        # The earlier model included them and inflated the class effects.
        d = d[d.n_de >= 1].copy()
        d["up"] = (d.n_up > 0).astype(int)
        d["log_breadth"] = np.log1p(d.breadth)
        # intergenic is the reference because it is the largest class and the
        # one with the lowest breadth, so every other class is read against it
        d["cls"] = pd.Categorical(d.cls, categories=CLS)
        st.note("genes={} up-rate={:.3f}".format(len(d), d.up.mean()))

        rows = []
        mod = smf.logit("up ~ C(cls) + log_breadth", data=d).fit(disp=0)
        for k in mod.params.index:
            if not k.startswith("C(cls)"):
                continue
            name = k.split("[T.")[-1].rstrip("]")
            rows.append(dict(term=name, coef=float(mod.params[k]),
                             odds_ratio=float(np.exp(mod.params[k])),
                             ci_low=float(np.exp(mod.conf_int().loc[k, 0])),
                             ci_high=float(np.exp(mod.conf_int().loc[k, 1])),
                             pval=float(mod.pvalues[k]),
                             model="logit up ~ class + log(breadth)"))
        rows.append(dict(term="log_breadth", coef=float(mod.params["log_breadth"]),
                         odds_ratio=float(np.exp(mod.params["log_breadth"])),
                         ci_low=float(np.exp(mod.conf_int().loc[
                             "log_breadth", 0])),
                         ci_high=float(np.exp(mod.conf_int().loc[
                             "log_breadth", 1])),
                         pval=float(mod.pvalues["log_breadth"]),
                         model="logit up ~ class + log(breadth)"))

        # unadjusted comparison, for the same contrasts, so the shrinkage is visible
        base = d[d.cls == "intergenic"]
        for c in CLS[1:]:
            s = d[d.cls == c]
            rows.append(dict(term=c, coef=np.nan,
                             odds_ratio=float((s.up.mean() / (1 - s.up.mean()))
                                              / (base.up.mean() /
                                                 (1 - base.up.mean()))),
                             ci_low=np.nan, ci_high=np.nan, pval=np.nan,
                             model="unadjusted odds ratio vs intergenic"))
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "D_breadth_adjusted.tsv")), sep="\t", index=False)

        adj = T[T.model.str.startswith("logit") & (T.term != "log_breadth")]
        st.note("breadth odds ratio per log unit {:.2f} (P={:.2g})".format(
            float(np.exp(mod.params["log_breadth"])),
            float(mod.pvalues["log_breadth"])))
        st.note("adjusted class odds ratios vs intergenic: " + ", ".join(
            "{}={:.2f}".format(r.term, r.odds_ratio) for r in adj.itertuples()))

        # matched comparison inside the middle breadth band
        lo, hi = d.breadth.quantile([0.4, 0.6])
        s = d[(d.breadth >= lo) & (d.breadth <= hi)]
        M = s.groupby("cls", observed=True).agg(
            n=("up", "size"), up_rate=("up", "mean"),
            breadth_med=("breadth", "median")).reset_index()
        M.to_csv(safe_write_path(os.path.join(
            OUT, "D_breadth_matched.tsv")), sep="\t", index=False)
        st.note("matched on breadth {:.0f}-{:.0f}: ".format(lo, hi)
                + ", ".join("{}={:.3f} (n={})".format(
                    r.cls, r.up_rate, r.n) for r in M.itertuples()))

        st.summary_md(os.path.join(RES, "step_summaries",
                                   "68_tumour_breadth.md"),
                      "68 tumour up-regulation adjusted for expression breadth",
                      df_to_md(T, ".4f") + "\n\nMatched band\n\n"
                      + df_to_md(M, ".4f"))


if __name__ == "__main__":
    main()
