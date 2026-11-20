#!/usr/bin/env python
"""88_gencode_biotype.py — benchmark the coordinate classes against GENCODE
transcript-level biotypes.

The five classes here are re-derived from coordinates. GENCODE v46 assigns
every transcript a transcript_type, and those labels are independent of this
pipeline. The gene-level gene_type is usually a single value, which is why the
manuscript describes annotation as collapsing subtypes, but the transcript
labels are informative. This module cross-tabulates the coordinate classes
against the dominant transcript type per gene and reports agreement, Cohen's
kappa and the sensitivity of each coordinate class to the matching label.

Outputs
  results/revision_checks/K_gencode_crosstab.tsv
  results/revision_checks/K_gencode_agreement.tsv
"""
import collections
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "revision_checks")
GTF = os.path.join(WORKDIR, "data", "gencode_v27",
                   "gencode.v27.annotation.gtf.gz")

# transcript_type to coordinate class
MAP = {
    "antisense": "antisense",
    "antisense_RNA": "antisense",
    "sense_intronic": "intronic",
    "sense_overlapping": "sense_overlapping",
    "lincRNA": "lincRNA",
    "lncRNA": "lincRNA",
    "macro_lncRNA": "lincRNA",
    "3prime_overlapping_ncRNA": "other_ncrna",
    "processed_transcript": "other_ncrna",
    "retained_intron": "other_ncrna",
    "misc_RNA": "other_ncrna",
    "TEC": "other_ncrna",
    "protein_coding": "protein_coding",
}


def attr(field, key):
    i = field.find(key + ' "')
    if i < 0:
        return None
    start = i + len(key) + 2
    j = field.find('"', start)
    return field[start:j]


def main():
    with step("88_gencode_biotype") as st:
        tx = collections.defaultdict(set)
        gt = {}
        with gzip.open(GTF, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                a = f[8]
                gi = attr(a, "gene_id")
                if gi is not None:
                    gi = gi.split(".")[0]
                if f[2] == "gene":
                    gt[gi] = attr(a, "gene_type")
                elif f[2] == "transcript":
                    tx[gi].add(attr(a, "transcript_type"))
        st.note("GENCODE genes={} with transcripts={}".format(len(gt), len(tx)))

        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m.index = m.index.astype(str).str.split(".").str[0]
        m = m[m.cls.isin(["divergent", "antisense", "intronic",
                          "sense_overlapping", "intergenic", "PCG",
                          "pseudogene", "transcribed_pseudogene"])]
        rows = []
        for gid, r in m.iterrows():
            t = tx.get(gid, set())
            if not t:
                rows.append((gid, r.cls, np.nan, np.nan))
                continue
            c = collections.Counter(MAP.get(x, "other_ncrna") for x in t)
            dom = c.most_common(1)[0][0]
            rows.append((gid, r.cls, dom, len(t)))
        D = pd.DataFrame(rows, columns=["gene_id", "cls", "gencode_dom",
                                        "n_types"])
        D = D.dropna(subset=["gencode_dom"])
        st.note("genes joined to GENCODE: {}".format(len(D)))

        ct = pd.crosstab(D.cls, D.gencode_dom)
        ct = ct.reindex(["PCG", "divergent", "antisense", "intronic",
                         "sense_overlapping", "intergenic", "pseudogene",
                         "transcribed_pseudogene"]).fillna(0).astype(int)
        frac = ct.div(ct.sum(1), axis=0)
        ct.to_csv(safe_write_path(os.path.join(
            OUT, "K_gencode_crosstab.tsv")), sep="\t")
        frac.to_csv(safe_write_path(os.path.join(
            OUT, "K_gencode_crosstab_frac.tsv")), sep="\t")

        # Cohen's kappa between coordinate class and GENCODE type, over
        # lncRNA genes with a GENCODE non-coding label. Divergent and
        # intergenic both map to lincRNA in GENCODE, so they are pooled.
        sub = D[D.gencode_dom != "protein_coding"].copy()
        sub["ours"] = sub.cls.replace({"divergent": "lincRNA",
                                       "intergenic": "lincRNA",
                                       "pseudogene": "other_ncrna",
                                       "transcribed_pseudogene": "other_ncrna"})
        labs = ["lincRNA", "antisense", "intronic", "sense_overlapping",
                "other_ncrna"]
        a = sub.ours.values
        b = sub.gencode_dom.values
        po = float((a == b).mean())
        pe = sum((a == L).mean() * (b == L).mean() for L in labs)
        kappa = (po - pe) / (1 - pe) if pe < 1 else np.nan

        # sensitivity and specificity of each coordinate class
        rows = []
        for ours, gtype in (("antisense", "antisense"),
                            ("intronic", "intronic"),
                            ("sense_overlapping", "sense_overlapping")):
            ours_n = sub[sub.ours == ours]
            spec = sub[sub.ours != ours]
            rows.append(dict(
                coordinate_class=ours, gencode_type=gtype,
                n=len(ours_n),
                sensitivity=float((ours_n.gencode_dom == gtype).mean()),
                specificity=float((spec.gencode_dom != gtype).mean()),
                prevalence=float((sub.gencode_dom == gtype).mean())))
        A = pd.DataFrame(rows + [dict(
            coordinate_class="all_lncRNA", gencode_type="any",
            n=len(sub), sensitivity=po, specificity=np.nan,
            prevalence=kappa)])
        A.to_csv(safe_write_path(os.path.join(
            OUT, "K_gencode_agreement.tsv")), sep="\t", index=False)
        for r in A.itertuples():
            if r.gencode_type == "any":
                st.note("overall agreement {:.1%}, kappa {:.3f}".format(
                    r.sensitivity, r.prevalence))
            else:
                st.note("{}: n={} sensitivity {:.1%} specificity {:.1%}"
                        .format(r.coordinate_class, r.n, r.sensitivity,
                                r.specificity))
        st.summary_md(os.path.join(
            RES, "step_summaries", "88_gencode_biotype.md"),
            "88 GENCODE transcript-type benchmark",
            df_to_md(ct) + "\n\n" + df_to_md(A, ".3f"))


if __name__ == "__main__":
    main()
