#!/usr/bin/env python
"""12b_te_subfamily_age.py — which TEs donate lncRNA TSSs, and how young are they?

Age bins curated from RepName consensus (documented in Methods):
  ancient      : L2, MIR, L3, CR1 (pre-eutherian)
  primate_old  : AluJ, AluS, L1PA3-7, L1PB, ERVL-MaLR (THE1, MLT), Charlie, ERV1
  primate_young: AluY, L1PA1/2, L1HS, LTR12/HERVK, SVA, MER51
Outputs: results/te_tss_donors_intergenic.tsv, results/te_age_by_class.tsv,
         results/te_ltr_donors.tsv, step summary
"""
import os
import sys
from collections import Counter

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data")

AGE_RULES = [
    ("ancient", ("L2", "L3", "MIR", "CR1", "MER2", "tRNA")),
    ("primate_young", ("AluY", "L1PA1", "L1PA2", "L1HS", "LTR12", "HERVK",
                       "SVA", "LTR13", "ERVK", "LTR5", "MER51")),
    ("primate_old", ("AluJ", "AluS", "L1PA3", "L1PA4", "L1PA5", "L1PA6",
                     "L1PA7", "L1PB", "THE1", "MST", "Charlie", "MLT1", "MLT2",
                     "ERV1", "ERVL", "LTR7", "LTR2", "MER4", "MER61")),
]


def age_of(sub):
    s = str(sub)
    for age, prefixes in AGE_RULES:
        if any(s.startswith(p) for p in prefixes):
            return age
    return "unclassified"


def main():
    with step("12b_te_subfamily_age") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                           index_col=0, low_memory=False)
        tt = pd.read_csv(os.path.join(DATA, "tss_in_te.tsv"), sep="\t",
                         names=["a_c", "a_s", "a_e", "gene_id", "t_c", "t_s",
                                "t_e", "t_id", "t_sub", "t_strand", "t_cls2",
                                "t_cls1"], header=None, usecols=[3, 8, 9, 10, 11],
                         dtype=str, engine="python")
        tt = tt.dropna().drop_duplicates("gene_id")
        tt["age"] = tt["t_sub"].map(age_of)
        cls_ser = meta["cls"].reindex(tt["gene_id"])
        cls_ser = cls_ser.fillna(meta["group"].reindex(tt["gene_id"])
                                 .map({"pcg": "PCG"}))
        tt["cls"] = cls_ser.values

        rows = []
        for c in ["divergent", "antisense", "intronic", "sense_overlapping",
                  "intergenic", "pseudogene", "transcribed_pseudogene", "PCG"]:
            sub = tt[tt.cls == c]
            if not len(sub):
                continue
            ac = Counter(sub["age"])
            rows.append(dict(cls=c, n_TSS_in_TE=len(sub),
                             young=ac.get("primate_young", 0),
                             old=ac.get("primate_old", 0),
                             ancient=ac.get("ancient", 0),
                             unclass=ac.get("unclassified", 0),
                             frac_young=ac.get("primate_young", 0) / len(sub),
                             frac_ancient=ac.get("ancient", 0) / len(sub)))
        age_tbl = pd.DataFrame(rows)
        age_tbl.to_csv(safe_write_path(os.path.join(RES, "te_age_by_class.tsv")),
                       sep="\t", index=False)
        st.note("age mix: " + "; ".join(
            f"{r.cls}: young={r.frac_young:.2f} ancient={r.frac_ancient:.2f}"
            for r in age_tbl.itertuples()))

        top = (tt[tt.cls == "intergenic"].groupby(["t_sub", "t_cls2"])
               .size().sort_values(ascending=False).head(20).reset_index(name="n"))
        top["age"] = top["t_sub"].map(age_of)
        top.to_csv(safe_write_path(
            os.path.join(RES, "te_tss_donors_intergenic.tsv")), sep="\t",
            index=False)
        st.note("top intergenic donors: " +
                ", ".join(f"{r.t_sub}({r.n})" for r in top.head(8).itertuples()))

        ltr = tt[tt.t_cls1 == "LTR"]
        ltr_rows = []
        for c in ["intergenic", "intronic", "antisense", "PCG"]:
            sub = ltr[ltr.cls == c]
            if not len(sub):
                continue
            fam = Counter(sub["t_cls2"])
            ltr_rows.append(dict(cls=c, n=len(sub),
                                 **{k: fam.get(k, 0) for k in ("ERVK", "ERV1", "ERVL", "LTR")}))
        ltr_tbl = pd.DataFrame(ltr_rows)
        ltr_tbl.to_csv(safe_write_path(os.path.join(RES, "te_ltr_donors.tsv")),
                       sep="\t", index=False)

        body = ("## TE age of TSS-donating elements\n\n" + df_to_md(age_tbl, ".3f")
                + "\n\n## Top TSS donors among intergenic lncRNAs\n\n"
                + df_to_md(top, ".0f", max_rows=20)
                + "\n\n## LTR superfamily donors\n\n" + df_to_md(ltr_tbl, ".0f"))
        st.summary_md("results/step_summaries/12b_te_subfamily_age.md",
                      "12b TE donors & age", body)


if __name__ == "__main__":
    main()
