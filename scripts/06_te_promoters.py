#!/usr/bin/env python
"""06_te_promoters.py — TE content & exaptation in promoters per class.

Uses LOCAL Dfam-based TE annotation (hg38.nrph.TE.anno.sorted.bed, cols:
chrom start end id subfamily subclass class) + HERV_noncoding.bed.
Outputs: results/te_promoter_stats.tsv (per gene), te_by_class.tsv,
         te_subfamily.tsv, step summary.
"""
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR, df_to_md

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data")
TE_BED = ("/media/share/resource/reference/homo_sapiens/GRCh38/annotation/"
          "hg38.nrph.TE.anno.sorted.bed")
HERV_BED = ("/media/share/resource/reference/homo_sapiens/GRCh38/annotation/"
            "HERV_noncoding.bed")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"

TE_ORDER_COL = {"LINE": "LINE", "SINE": "SINE", "DNA": "DNA", "LTR": "LTR",
                "RC": "RC", "Unknown": "Unknown", "Other": "Other"}


def bt_to_file(outfile, *args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}", SIF_BT, "bedtools", *args]
    with open(safe_write_path(outfile), "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bedtools failed: {' '.join(args)}\n{r.stderr[:1000]}")


def main():
    with step("06_te_promoters") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                           index_col=0)
        prom = pd.read_csv(os.path.join(RES, "promoters_feat.bed"), sep="\t",
                           names=["chrom", "start", "end", "gene_id", "score", "strand"])

        # ---- per-promoter TE coverage & subfamily hits (-500/+200 window)
        # NOTE: separate HERV_noncoding.bed dropped — the Dfam (nrph) annotation
        # already covers ERV/LTR families (that file is BED12 anyway).
        ov_file = os.path.join(DATA, "te_prom_overlaps.tsv")
        bt_to_file(ov_file, "intersect", "-a", os.path.join(RES, "promoters_feat.bed"),
                   "-b", TE_BED, "-wa", "-wb")

        cols = ["chrom", "start", "end", "gene_id", "score", "strand",
                "t_chrom", "t_start", "t_end", "t_id", "t_sub", "t_strand",
                "t_cls2", "t_cls1"]
        ov = pd.read_csv(ov_file, sep="\t", names=cols, usecols=cols, header=None,
                         dtype={"t_sub": str, "t_cls2": str, "t_cls1": str},
                         on_bad_lines="skip")
        st.note(f"TE overlap rows={len(ov)}")

        # TE order mapping: col8 subclass (e.g. L1, Alu, ERVK...) col9 class (LINE/SINE/...)
        def order_of(row):
            c1 = str(row["t_cls1"])
            if c1 in TE_ORDER_COL:
                return c1
            return "Other"

        ov["order"] = ov.apply(order_of, axis=1)
        # overlap length per row
        ov["ov_len"] = ov.apply(lambda r: min(r["end"], r["t_end"]) -
                                max(r["start"], r["t_start"]), axis=1)
        ov = ov[ov["ov_len"] > 0]

        # per-promoter coverage fraction (700bp window) & orders present
        prom_len = prom.set_index("gene_id")["end"] - prom.set_index("gene_id")["start"]
        cov = ov.groupby("gene_id")["ov_len"].sum() / prom_len.clip(lower=1)
        orders = ov.groupby(["gene_id", "order"])["ov_len"].sum().unstack(fill_value=0)
        orders = orders.div(orders.sum(axis=1), axis=0).fillna(0)
        ltr_frac = orders.get("LTR", pd.Series(dtype=float))

        # TSS inside TE (1bp TSS point): intersect promoter-START point with TEs
        # The point must be the major TSS itself. The earlier version used the
        # promoter-window edge (start+1 or end), which sits up to 500 bp away
        # and gave a TSS-in-TE rate that was not a TSS property.
        tsspt = prom[["chrom", "gene_id", "strand"]].copy()
        tsspt["p"] = tsspt.gene_id.map(meta.tss)
        tsspt = tsspt.dropna(subset=["p"])
        tsspt["p"] = tsspt.p.astype(int)
        tsspt["s"] = tsspt["p"] - 1
        tsspt[["chrom", "s", "p", "gene_id"]].to_csv(
            safe_write_path(os.path.join(DATA, "tss_points.bed")),
            sep="\t", index=False, header=False)
        tss_te_file = os.path.join(DATA, "tss_in_te.tsv")
        bt_to_file(tss_te_file, "intersect", "-a", os.path.join(DATA, "tss_points.bed"),
                   "-b", TE_BED, "-wa", "-wb")
        tt = pd.read_csv(tss_te_file, sep="\t", header=None,
                         names=["a_c", "a_s", "a_e", "gene_id",
                                "t_c", "t_s", "t_e", "t_id", "t_sub",
                                "t_strand", "t_cls2", "t_cls1"],
                         usecols=[3, 8, 10, 11], dtype=str, on_bad_lines="skip")
        tt.columns = ["gene_id", "t_sub", "t_cls2", "t_cls1"]
        tss_in_te = tt["gene_id"].nunique() if len(tt) else 0
        st.note(f"TSS inside TE (unique genes)={tss_in_te}")
        # donor subfamily for TSS-in-TE
        if len(tt):
            donor = tt.groupby(["t_sub"]).size().sort_values(ascending=False)
            st.note("top TSS-donor TEs: " + ", ".join(
                f"{k}={v}" for k, v in donor.head(8).items()))

        # ---- assemble per-gene table
        te = pd.DataFrame(index=meta.index)
        te["te_cov"] = cov.reindex(te.index).fillna(0)
        for o in ("LINE", "SINE", "DNA", "LTR", "RC", "Other"):
            te[f"te_frac_{o}"] = orders[o].reindex(te.index).fillna(0) if o in orders else 0
        tt_set = set(tt["gene_id"]) if len(tt) else set()
        te["tss_in_te"] = [1 if g in tt_set else 0 for g in te.index]
        te["tss_in_ltr"] = [
            1 if g in set(tt[tt["t_cls1"] == "LTR"]["gene_id"]) else 0
            for g in te.index] if len(tt) else 0

        out = meta.join(te, how="left")
        out.to_csv(safe_write_path(os.path.join(RES, "te_promoter_stats.tsv")), sep="\t")

        # by class
        groups = ["PCG", "divergent", "antisense", "intronic", "sense_overlapping",
                  "intergenic", "pseudogene", "transcribed_pseudogene"]
        rows = []
        for g in groups:
            sub = out[(out["group"] == "pcg")] if g == "PCG" else out[out["cls"] == g]
            if not len(sub):
                continue
            rows.append(dict(cls=g, n=len(sub), te_cov=sub["te_cov"].mean(),
                             frac_te=sub["te_cov"].ge(0.05).mean(),
                             LTR=sub["te_frac_LTR"].mean(), LINE=sub["te_frac_LINE"].mean(),
                             SINE=sub["te_frac_SINE"].mean(), DNA=sub["te_frac_DNA"].mean(),
                             tss_in_te=sub["tss_in_te"].mean(),
                             tss_in_ltr=sub["tss_in_ltr"].mean()))
        sm = pd.DataFrame(rows)
        sm.to_csv(safe_write_path(os.path.join(RES, "te_by_class.tsv")),
                  sep="\t", index=False)
        st.note("TE by class: " + "; ".join(
            f"{r.cls}: cov={r.te_cov:.3f} tssTE={r.tss_in_te:.3f} LTR={r.LTR:.3f}"
            for r in sm.itertuples()))
        st.summary_md("results/step_summaries/06_te_promoters.md",
                      "06 TE exaptation", df_to_md(sm, ".3f"))


if __name__ == "__main__":
    main()
