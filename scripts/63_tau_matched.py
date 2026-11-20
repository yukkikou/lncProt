#!/usr/bin/env python
"""63_tau_matched.py — is the class difference in GTEx tau an expression artefact?

tau rises towards 1 whenever a gene is detected in few tissues, and lncRNAs are
detected in fewer tissues than protein-coding genes because they are expressed
more weakly. The class comparison in the Results therefore risks restating an
expression difference. This module matches every lncRNA to protein-coding genes
of the same expression level and recomputes tau.

Matching is nearest neighbour on log10(mean expression of the strongest tissue),
within a caliper of 0.1 log units, one-to-one without replacement.

Outputs
  results/tissue_disease/tables/tau_expression_matched.tsv
  results/tissue_disease/tables/tau_by_expression_bin.tsv
  results/step_summaries/63_tau_matched.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
TD = os.path.join(RES, "tissue_disease", "tables")
CLS = ["divergent", "antisense", "intronic", "sense_overlapping",
       "intergenic"]
CALIPER = 0.1


def match_one_to_one(target, pool, caliper):
    """Nearest neighbour on the matching variable, without replacement.

    Both vectors are already log10 expression. A greedy pass from the rarest
    target gene outwards keeps the matching honest at the low end, where most
    lncRNAs sit.
    """
    used = np.zeros(len(pool), dtype=bool)
    order = np.argsort(pool)
    pool_sorted = pool[order]
    idx_out = np.full(len(target), -1, dtype=int)
    for i in np.argsort(target):
        j = np.searchsorted(pool_sorted, target[i])
        best, bestd = -1, np.inf
        for k in (j - 1, j, j + 1, j - 2, j + 2):
            if 0 <= k < len(pool_sorted):
                pidx = order[k]
                if used[pidx]:
                    continue
                d = abs(pool_sorted[k] - target[i])
                if d < bestd and d <= caliper:
                    best, bestd = pidx, d
        if best >= 0:
            used[best] = True
            idx_out[i] = best
    return idx_out


def expression_bins(g):
    """tau as a function of expression, for PCGs alone.

    This curve is the null a class comparison has to beat.
    """
    pcg = g[g.cls == "PCG"].copy()
    pcg["bin"] = pd.qcut(pcg.logexpr, 10, labels=False, duplicates="drop")
    return pcg.groupby("bin").agg(
        n=("tau_gtex", "size"),
        expr_med=("logexpr", "median"),
        tau_med=("tau_gtex", "median")).reset_index()


def main():
    with step("63_tau_matched") as st:
        g = pd.read_csv(os.path.join(TD, "gtex_gene_tau.tsv"), sep="\t")
        g = g.rename(columns={"Unnamed: 0": "gene_id"})
        g = g.dropna(subset=["tau_gtex"])
        g = g[g.max_tissue_mean > 0]
        g["logexpr"] = np.log10(g.max_tissue_mean)
        st.note("genes with tau={} (PCG={})".format(
            len(g), int((g.cls == "PCG").sum())))

        pcg = g[g.cls == "PCG"].reset_index(drop=True)
        rows = []
        for c in CLS:
            sub = g[g.cls == c].reset_index(drop=True)
            if len(sub) < 20:
                continue
            m = match_one_to_one(sub.logexpr.values, pcg.logexpr.values,
                                 CALIPER)
            ok = m >= 0
            if ok.sum() < 20:
                continue
            lt = sub.tau_gtex.values[ok]
            pt = pcg.tau_gtex.values[m[ok]]
            rows.append(dict(cls=c, n_lnc=len(sub), n_matched=int(ok.sum()),
                             tau_lnc_median=float(np.median(lt)),
                             tau_pcg_matched=float(np.median(pt)),
                             tau_pcg_all=float(np.median(pcg.tau_gtex)),
                             median_diff_matched=float(np.median(lt - pt)),
                             expr_lnc_median=float(np.median(
                                 sub.logexpr.values[ok])),
                             expr_pcg_median=float(np.median(
                                 pcg.logexpr.values[m[ok]]))))
        R = pd.DataFrame(rows)
        if len(R):
            w, pw = wilcoxon(R.median_diff_matched) if len(R) > 5 else (np.nan,
                                                                       np.nan)
            st.note("per-class matched differences: " + ", ".join(
                "{}={:+.3f}".format(r.cls, r.median_diff_matched)
                for r in R.itertuples()))
        R.to_csv(safe_write_path(os.path.join(
            TD, "tau_expression_matched.tsv")), sep="\t", index=False)

        B = expression_bins(g)
        B.to_csv(safe_write_path(os.path.join(
            TD, "tau_by_expression_bin.tsv")), sep="\t", index=False)
        st.note("tau rises from {:.3f} in the lowest expression decile to "
                "{:.3f} in the highest among PCGs".format(
                    B.tau_med.iloc[0], B.tau_med.iloc[-1]))

        # pooled test: all lncRNAs with a matched PCG against that PCG
        allm = []
        for c in CLS:
            sub = g[g.cls == c].reset_index(drop=True)
            m = match_one_to_one(sub.logexpr.values, pcg.logexpr.values,
                                 CALIPER)
            ok = m >= 0
            for a, b in zip(sub.tau_gtex.values[ok], pcg.tau_gtex.values[m[ok]]):
                allm.append((c, a, b))
        A = pd.DataFrame(allm, columns=["cls", "tau_lnc", "tau_pcg"])
        if len(A) > 20:
            w, pw = wilcoxon(A.tau_lnc, A.tau_pcg)
            st.note("pooled matched comparison n={} median diff {:+.3f} "
                    "P={:.3g}".format(len(A), float(np.median(
                        A.tau_lnc - A.tau_pcg)), pw))
            A.to_csv(safe_write_path(os.path.join(
                TD, "tau_matched_pairs.tsv")), sep="\t", index=False)

        # ---- developmental breadth, same confound, same control
        dev_path = os.path.join(RES, "developmental", "per_gene_summary.tsv")
        if os.path.exists(dev_path):
            d = pd.read_csv(dev_path, sep="\t")
            h = d[d.species == "human"].copy()
            master = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                                 sep="\t", index_col=0, low_memory=False)
            master["base"] = master.index.str.split(".").str[0]
            h["base"] = h.gene.str.split(".").str[0]
            h = h.merge(master[["base", "group"]].drop_duplicates("base"),
                        on="base", how="left")
            h = h.dropna(subset=["breadth_frac", "mean_tpm", "group"])
            h["logexpr"] = np.log10(h.mean_tpm.clip(lower=1e-3))
            pcgh = h[h.group == "pcg"]
            x = pcgh.logexpr.values
            y = pcgh.breadth_frac.values
            order = np.argsort(x)
            xs, ys = x[order], y[order]
            lnch = h[h.group == "lnc"]
            j = np.clip(np.searchsorted(xs, lnch.logexpr.values), 0,
                        len(xs) - 1)
            jl = np.clip(j - 1, 0, len(xs) - 1)
            take = np.where(np.abs(xs[j] - lnch.logexpr.values) <=
                            np.abs(xs[jl] - lnch.logexpr.values), j, jl)
            okm = np.abs(xs[take] - lnch.logexpr.values) <= CALIPER
            D = pd.DataFrame([dict(
                metric="developmental breadth fraction",
                n_matched=int(okm.sum()),
                lnc_median=float(np.median(
                    lnch.breadth_frac.values[okm])),
                pcg_matched_median=float(np.median(ys[take][okm])),
                pcg_unmatched_median=float(np.median(y)),
                lnc_unmatched_median=float(lnch.breadth_frac.median()))])
            D.to_csv(safe_write_path(os.path.join(
                TD, "developmental_breadth_matched.tsv")), sep="\t",
                index=False)
            st.note("developmental breadth: lnc {:.3f} versus matched PCG "
                    "{:.3f} (unmatched PCG {:.3f})".format(
                        D.lnc_median.iloc[0], D.pcg_matched_median.iloc[0],
                        D.pcg_unmatched_median.iloc[0]))

        st.summary_md(os.path.join(RES, "step_summaries", "63_tau_matched.md"),
                      "63 expression-matched GTEx tau",
                      df_to_md(R, ".3f") + "\n\ntau by PCG expression decile\n\n"
                      + df_to_md(B, ".3f"))


if __name__ == "__main__":
    main()
