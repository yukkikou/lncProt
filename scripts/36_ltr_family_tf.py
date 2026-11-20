#!/usr/bin/env python
"""36_ltr_family_tf.py — LTR family x TF motif modules by birth class.

Per-gene TSS-donating LTR subfamily (data/tss_in_te.tsv, order=LTR) crossed
with motif presence (results/motif_presence.tsv):
  A. LTR family enrichment by birth class
  B. LTR family x TF motif association within LTR-derived promoters
  C. class prediction from family + TF motifs (5-fold CV)
  D. LTR family age (12b rules) and phyloP
Outputs -> results/ltr_family/*.tsv + step summary.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "ltr_family")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
AGE_RULES = [
    ("ancient", ("L2", "L3", "MIR", "CR1", "MER2", "tRNA")),
    ("primate_young", ("AluY", "L1PA1", "L1PA2", "L1HS", "LTR12", "HERVK",
                       "SVA", "LTR13", "ERVK", "LTR5", "MER51")),
    ("primate_old", ("AluJ", "AluS", "L1PA3", "L1PA4", "L1PA5", "L1PA6",
                     "L1PA7", "L1PB", "THE1", "MST", "Charlie", "MLT1", "MLT2",
                     "ERV1", "ERVL", "LTR7", "LTR2", "MER4", "MER61")),
]
FOCUS_TF = ["DUX4", "DUX4_HUMAN", "PHOX2B", "BATF3", "CTCF", "CTCF_HUMAN",
            "MAZ", "E2F6", "TFAP2A", "TFAP2C"]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def age_of(sub):
    s = str(sub)
    for age, prefixes in AGE_RULES:
        if any(s.startswith(p) for p in prefixes):
            return age
    return "unclassified"


def main():
    with step("36_ltr_family") as st:
        tt = pd.read_csv(os.path.join(WORKDIR, "data", "tss_in_te.tsv"),
                         sep="\t", header=None,
                         names=["c", "s", "e", "gene_id", "tc", "ts", "te",
                                "tid", "t_sub", "tstrand", "fam", "order"])
        tt = tt[tt.order == "LTR"].copy()
        tt["ov"] = tt.te.astype(int) - tt.ts.astype(int)
        tt = tt.sort_values("ov", ascending=False).drop_duplicates("gene_id")
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)[
            ["cls", "group"]]
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["tss_in_ltr"]]
        mot = pd.read_csv(os.path.join(RES, "motif_presence.tsv"), sep="\t",
                          index_col=0, low_memory=False)
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["phylop_mean"]]
        df = tt.set_index("gene_id").join(meta).join(te).join(ph)
        df = df[df.cls.isin(CLASSES)].copy()
        st.note("LTR-TSS genes: {}; families: {}".format(
            len(df), df.t_sub.nunique()))
        fam_n = df.t_sub.value_counts()
        keep = fam_n[fam_n >= 20].index
        df["family"] = np.where(df.t_sub.isin(keep), df.t_sub, "other")
        fams = [f for f in keep]
        st.note("families kept (n>=20): {}".format(len(fams)))

        # A. family x class enrichment
        rows = []
        for f in fams:
            in_f = df.family == f
            for c in CLASSES:
                a = int((in_f & (df.cls == c)).sum())
                b = int((in_f & (df.cls != c)).sum())
                cc = int((~in_f & (df.cls == c)).sum())
                dd = int((~in_f & (df.cls != c)).sum())
                if a + b < 20 or a + cc < 10:
                    continue
                orr, p = fisher_exact([[a, b], [cc, dd]])
                rows.append(dict(family=f, cls=c, n_family=a + b,
                                 n_class=a + cc,
                                 frac_in_family=a / max(a + b, 1),
                                 frac_elsewhere=cc / max(cc + dd, 1),
                                 log2OR=float(np.log2(orr)) if orr > 0 else np.nan,
                                 pval=p))
        A = pd.DataFrame(rows)
        A["fdr"] = bh_fdr(A.pval.values)
        A = A.sort_values(["cls", "log2OR"], ascending=[True, False])
        A.to_csv(safe_write_path(os.path.join(OUT, "family_by_class.tsv")),
                 sep="\t", index=False)

        # B. family x TF within LTR genes
        common = [c for c in mot.columns
                  if mot[c].notna().mean() > 0.5 and mot[c].mean() >= 0.05]
        msub = mot.loc[mot.index.intersection(df.index), common]
        d2 = df.join(msub)
        rows = []
        for f in fams:
            in_f = d2.family == f
            if in_f.sum() < 20:
                continue
            for tf in common:
                v = d2[tf].fillna(0).astype(int) == 1
                a = int((in_f & v).sum())
                b = int((in_f & ~v).sum())
                cc = int((~in_f & v).sum())
                dd = int((~in_f & ~v).sum())
                if a < 5:
                    continue
                orr, p = fisher_exact([[a, b], [cc, dd]])
                rows.append(dict(family=f, tf=tf, n_family=int(in_f.sum()),
                                 frac_family=a / max(a + b, 1),
                                 frac_other=cc / max(cc + dd, 1),
                                 log2OR=float(np.log2(orr)) if orr > 0 else np.nan,
                                 pval=p))
        B = pd.DataFrame(rows)
        B["fdr"] = bh_fdr(B.pval.values)
        B = B.sort_values(["family", "log2OR"], ascending=[True, False])
        B.to_csv(safe_write_path(os.path.join(OUT, "family_tf_pairs.tsv")),
                 sep="\t", index=False)

        # C. class prediction from family / TF / both
        d3 = d2.copy()
        fam_dum = pd.get_dummies(d3.family, prefix="fam")
        tf_mat = d3[common].fillna(0).astype(float)
        y = d3.cls.values
        skf = StratifiedKFold(5, shuffle=True, random_state=42)
        mrows = []
        for label, X in [("family", fam_dum), ("TF", tf_mat),
                         ("family+TF", pd.concat([fam_dum, tf_mat], axis=1))]:
            clf = make_pipeline(StandardScaler(with_mean=False),
                                LogisticRegression(max_iter=3000,
                                                   multi_class="multinomial"))
            try:
                sc = cross_val_score(clf, X.values, y, cv=skf,
                                     scoring="f1_macro")
                mrows.append(dict(feature_set=label, n_features=X.shape[1],
                                  macro_f1_mean=float(np.mean(sc)),
                                  macro_f1_sd=float(np.std(sc))))
            except Exception as e:
                mrows.append(dict(feature_set=label, n_features=X.shape[1],
                                  macro_f1_mean=np.nan, macro_f1_sd=np.nan))
        M = pd.DataFrame(mrows)
        M.to_csv(safe_write_path(os.path.join(OUT, "class_model.tsv")),
                 sep="\t", index=False)

        # D. family age
        d3["age"] = d3.t_sub.map(age_of)
        age_tab = pd.crosstab(d3.family, d3.age)
        age_tab.to_csv(safe_write_path(os.path.join(OUT, "family_age.tsv")),
                       sep="\t")
        age_med = d3.groupby("age").phylop_mean.median()
        r, pr = spearmanr(d3.family.map(
            {f: i for i, f in enumerate(fams)}), d3.phylop_mean,
            nan_policy="omit")
        st.note("top family-class: " + "; ".join(
            "{}:{} log2OR={:.2f} q={:.1e}".format(
                r.cls, r.family, r.log2OR, r.fdr)
            for r in A.sort_values("fdr").head(4).itertuples()))
        st.note("family-TF top: " + "; ".join(
            "{}+{} log2OR={:.2f} q={:.1e}".format(
                r.family, r.tf, r.log2OR, r.fdr)
            for r in B.sort_values("fdr").head(4).itertuples()))
        st.note("class model: " + "; ".join(
            "{} F1={:.3f}".format(r.feature_set, r.macro_f1_mean)
            for r in M.itertuples()))
        st.summary_md("results/ltr_family/step_summaries/36_ltr.md",
                      "36 LTR family x TF modules",
                      "## Family by class (top)\n\n"
                      + df_to_md(A.sort_values(["cls", "log2OR"],
                                               ascending=[True, False])
                                 .groupby("cls").head(4), ".3f")
                      + "\n\n## Family-TF pairs (top)\n\n"
                      + df_to_md(B.sort_values("fdr").head(25), ".3f")
                      + "\n\n## Class prediction\n\n" + df_to_md(M, ".3f")
                      + "\n\n## Family age\n\n"
                      + df_to_md(age_tab.reset_index(), "d"))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
