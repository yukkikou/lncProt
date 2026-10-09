#!/usr/bin/env python
"""15_final_gaps.py — publication gap-filling:
 (a) FANTOM tau vs GTEx tau concordance (fixed index alignment)
 (b) LTR-lncRNA enrichment in testis (Fisher vs all other tissues)
 (c) pan-cancer LTR derepression (permutation of cancer labels on delta)
 (d) DUX4/PHOX2B/BATF3 ChIP (ReMap) overlap with LTR-derived lnc promoters
 (e) gene-level normal-activity vs tumour-dysregulation association
Outputs -> results/final_gaps/*.tsv + step summary
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

TD = os.path.join(WORKDIR, "data", "tissue_disease")
OUT = os.path.join(WORKDIR, "results", "final_gaps")
RES = os.path.join(WORKDIR, "results")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def main():
    with step("15_final_gaps") as st:
        os.makedirs(OUT, exist_ok=True)
        genes = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                            index_col=0)

        # ---- (a) FANTOM tau vs GTEx tau (aligned by base ENSG)
        tau_g = pd.read_csv(os.path.join(OUT, "..", "tissue_disease", "tables",
                                         "gtex_gene_tau.tsv"), sep="\t",
                            index_col=0)
        fant = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                           index_col=0, low_memory=False)[["tau"]]
        fant = fant[fant.tau.notna()]
        f_idx = pd.Index([g.split(".")[0] for g in fant.index])
        f_ser = pd.Series(fant.tau.values, index=f_idx)
        f_ser = f_ser[~f_ser.index.duplicated(keep="first")]
        common = tau_g.index.intersection(f_ser.index)
        rho, p = spearmanr(tau_g.loc[common, "tau_gtex"],
                           f_ser.loc[common])
        pd.DataFrame({"rho": [rho], "p": [p], "n": [len(common)]}).to_csv(
            safe_write_path(os.path.join(OUT, "fantom_vs_gtex_tau.tsv")),
            sep="\t", index=False)
        st.note(f"(a) FANTOM~GTEx tau rho={rho:.3f} p={p:.1e} n={len(common)}")

        # ---- (b) testis LTR enrichment (per-gene Fisher)
        gt = pd.read_csv(os.path.join(OUT, "..", "tissue_disease", "tables",
                                      "gtex_class_by_tissue.tsv"), sep="\t")
        testis = gt[gt.tissue.str.contains("Testis")]
        others = gt[~gt.tissue.str.contains("Testis")]
        exp_lnc_testis = testis.share_LTR.iloc[0] * testis.n_lnc_expressed.iloc[0]
        tot_testis = testis.n_lnc_expressed.iloc[0]
        exp_lnc_other = others.share_LTR.mean() * others.n_lnc_expressed.sum()
        tot_other = others.n_lnc_expressed.sum()
        orr, p = fisher_exact([[round(exp_lnc_testis), round(tot_testis - exp_lnc_testis)],
                               [round(exp_lnc_other), round(tot_other - exp_lnc_other)]])
        st.note(f"(b) testis LTR share {testis.share_LTR.iloc[0]:.3f} vs "
                f"others {others.share_LTR.mean():.3f}; OR={orr:.2f} p={p:.1e}")

        # ---- (c) pan-cancer LTR derepression (label permutation)
        lt = pd.read_csv(os.path.join(OUT, "..", "tissue_disease", "tables",
                                      "tcga_ltr_derepression.tsv"), sep="\t")
        obs = lt["diff"].mean()
        rng = np.random.default_rng(11)
        # permute LTR labels at gene level within each cancer is heavy;
        # here bootstrap across cancers as the pan-cancer statistic
        boot = [lt["diff"].sample(frac=1.0, replace=True, random_state=rng)
                .mean() for _ in range(10000)]
        lo, hi = np.percentile(boot, [2.5, 97.5])
        p_perm = float(np.mean(np.array(boot) <= 0))
        pd.DataFrame({"mean_diff": [obs], "ci_lo": [lo], "ci_hi": [hi],
                      "p_boot_vs0": [p_perm], "n_cancers": [len(lt)]}).to_csv(
            safe_write_path(os.path.join(OUT, "pancancer_ltr_derepression.tsv")),
            sep="\t", index=False)
        st.note(f"(c) pan-cancer LTR delta={obs:.3f} 95%CI[{lo:.3f},{hi:.3f}] "
                f"P(<=0)={p_perm:.3f}")

        # ---- (d) candidate TF ChIP overlap with LTR-derived promoters
        rh = pd.read_csv(os.path.join(WORKDIR, "data", "remap_hits.tsv.gz"),
                         sep="\t", names=["gene_id", "tf", "cell", "peak"])
        rh["base"] = rh.gene_id.str.split(".").str[0]
        ltr_set = set(genes.index[(genes.tss_in_ltr == 1).fillna(False)])
        lnc_bases = set(genes.index[genes.cls.isin(CLASSES)])
        rows = []
        for tf in ["DUX4", "PHOX2B", "BATF3", "CTCF", "MAZ"]:
            bound = set(rh[rh.tf == tf]["base"])
            if not bound:
                continue
            b_ltr = len(bound & ltr_set)
            b_all = len(bound & lnc_bases)
            n_ltr = len(ltr_set)
            n_all = len(lnc_bases)
            orr, p = fisher_exact([[b_ltr, b_all - b_ltr],
                                   [n_ltr - b_ltr, (n_all - n_ltr) - (b_all - b_ltr)]])
            rows.append(dict(tf=tf, bound_lnc=b_all, bound_LTR_derived=b_ltr,
                             frac_LTR_among_bound=b_ltr / max(b_all, 1),
                             background_LTR_frac=n_ltr / n_all,
                             fisher_OR=orr, p=p))
        tfd = pd.DataFrame(rows)
        tfd.to_csv(safe_write_path(os.path.join(OUT, "tf_binding_vs_LTR.tsv")),
                   sep="\t", index=False)
        st.note("(d) " + "; ".join(f"{r.tf}: LTR-frac bound={r.frac_LTR_among_bound:.3f} "
                f"bg={r.background_LTR_frac:.3f} p={r.p:.1e}"
                for r in tfd.itertuples()))

        # ---- (e) gene-level normal activity vs tumour change
        de = pd.read_csv(os.path.join(OUT, "..", "tissue_disease", "tables",
                                      "tcga_class_de_per_cancer.tsv"), sep="\t")
        # per class: correlation between GTEx median frac_expr and TCGA median frac_up (n=5)
        rows = []
        for c in CLASSES:
            sub = de[de.cls == c]
            rows.append(dict(cls=c,
                             gtex_median_frac_expr=gt[f"{c}_frac_expr"].median(),
                             tcga_median_frac_up=sub.frac_up.median()))
        iv = pd.DataFrame(rows)
        rho2, p2 = spearmanr(iv.gtex_median_frac_expr, iv.tcga_median_frac_up)
        iv.to_csv(safe_write_path(os.path.join(OUT, "activity_vs_dysregulation.tsv")),
                  sep="\t", index=False)
        st.note(f"(e) normal-activity vs dysregulation rho={rho2:.2f} p={p2:.2f} "
                f"(class level, n={len(iv)})")

        st.summary_md("results/step_summaries/15_final_gaps.md",
                      "15 final gap-filling statistics",
                      f"- FANTOM~GTEx tau: rho={rho:.3f}, p={p:.1e}, n={len(common)}\n"
                      f"- Testis LTR enrichment: OR={orr:.2f}, p={p:.1e} (per-tissue aggregate)\n"
                      f"- Pan-cancer LTR derepression: {obs:.3f} [{lo:.3f},{hi:.3f}], P={p_perm:.3f}\n\n"
                      "## TF binding vs LTR-derived promoters\n\n"
                      + df_to_md(tfd, ".3f")
                      + f"\n\n## Class activity vs dysregulation\n\n" + df_to_md(iv, ".3f")
                      + f"\nrho={rho2:.2f} (p={p2:.2f})")


if __name__ == "__main__":
    main()
