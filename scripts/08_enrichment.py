#!/usr/bin/env python
"""08_enrichment.py — motif & TF-binding enrichment per lncRNA class vs PCG.

Inputs : results/gene_master_table.tsv, data/fimo/chunk_*/fimo.hits.tsv.gz,
         data/remap_hits.tsv.gz (optional)
Outputs: results/motif_presence.npy/.tsv   (gene x motif presence matrix)
         results/motif_enrichment.tsv      (class vs PCG Fisher + BH-FDR)
         results/motif_position_bias.tsv
         results/tf_binding_enrichment.tsv
         results/candidate_TFs.tsv         (motif AND binding supported)
"""
import glob
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR, df_to_md

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping", "intergenic"]


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    ranked = p[o] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(ranked, 0, 1)
    return out


def norm_tf_symbol(motif_id, alt):
    """Extract TF symbol from HOCOMOCO/JASPAR identifiers."""
    s = "" if alt is None or (isinstance(alt, float)) else str(alt).strip()
    if s and not s.startswith("NA"):
        s = s.split("_")[0]
        return s
    mid = "" if motif_id is None else str(motif_id)
    return mid.split(".")[0]


def load_fimo():
    frames = []
    for f in sorted(glob.glob(os.path.join(DATA, "fimo", "chunk_*", "fimo.hits.tsv.gz"))):
        with gzip.open(f, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
        df = pd.read_csv(f, sep="\t", compression="gzip", skiprows=1,
                         names=head[:9],
                         usecols=[0, 1, 2, 3],   # motif_id, alt, gene, start only
                         dtype={"motif_id": str, "motif_alt_id": str,
                                "sequence_name": str, "start": "float64"})
        df = df.dropna(subset=["sequence_name", "start"])
        df["start"] = df["start"].astype("int32")
        frames.append(df)
    if not frames:
        raise FileNotFoundError("no FIMO chunk outputs found")
    allhits = pd.concat(frames, ignore_index=True)
    allhits = allhits.rename(columns={"sequence_name": "gene_id"})
    allhits["tf"] = [norm_tf_symbol(m, a) for m, a in
                     zip(allhits["motif_id"], allhits["motif_alt_id"])]
    return allhits


def fisher_vs_pcg(presence, meta, index_genes, min_hits=5):
    """presence: bool Series gene->has-hit. Returns per class rows."""
    rows = []
    pcg = meta[meta["group"] == "pcg"]
    a_p = presence.reindex(pcg.index).fillna(False)
    for cls in CLASSES + ["pseudogene", "transcribed_pseudogene"]:
        sub = meta[meta["cls"] == cls]
        if not len(sub):
            continue
        a_c = presence.reindex(sub.index).fillna(False)
        n11 = int((a_c & True).sum()); n10 = int((~a_c).sum())
        n01 = int(a_p.sum());        n00 = int((~a_p).sum())
        if n11 < min_hits and n01 < min_hits:
            continue
        orr, p = fisher_exact([[n11, n10], [n01, n00]])
        rows.append(dict(cls=cls, n_class=len(sub), hits_class=n11,
                         frac_class=n11 / max(len(sub), 1), hits_pcg=n01,
                         frac_pcg=n01 / max(len(pcg), 1),
                         log2OR=np.log2(orr) if orr > 0 else np.nan, pval=p))
    return rows


def main():
    with step("08_enrichment") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                           index_col=0)

        # ---------------- motifs ------------------------------------------------
        hits = load_fimo()
        st.note(f"fimo hits={len(hits):,} motifs={hits['motif_id'].nunique()}")
        # presence per (gene, tf_symbol) — union over motif variants
        pres = hits[["gene_id", "tf"]].drop_duplicates()
        piv = pres.assign(v=1).pivot_table(index="gene_id", columns="tf",
                                           values="v", fill_value=0)
        st.note(f"presence matrix: {piv.shape[0]} genes x {piv.shape[1]} TFs")
        piv.to_csv(safe_write_path(os.path.join(RES, "motif_presence.tsv")),
                   sep="\t")   # sparse-ish; acceptable at this size

        rows = []
        for tf in piv.columns:
            presence = piv[tf].astype(bool)
            rows.extend([dict(tf=tf, **r) for r in fisher_vs_pcg(presence, meta, None)])
        me = pd.DataFrame(rows)
        if len(me):
            me["fdr"] = bh_fdr(me["pval"].values)
        me.to_csv(safe_write_path(os.path.join(RES, "motif_enrichment.tsv")),
                  sep="\t", index=False)
        st.note(f"motif_enrichment rows={len(me)}, "
                f"sig(FDR<0.05)={int((me['fdr'] < 0.05).sum())}")

        # position bias for top-5 enriched motifs per class (KS vs uniform)
        from scipy.stats import kstest
        pb_rows = []
        seq_len = 600  # fimo window -500/+100
        top = (me[me["fdr"] < 0.05].sort_values("log2OR", ascending=False)
               .groupby("cls").head(5))
        for _, r in top.iterrows():
            sub = hits[(hits["tf"] == r["tf"])]
            genes = set(meta[(meta["cls"] == r["cls"])].index) | \
                    set(meta[(meta["group"] == "pcg")].index) if r["cls"] != "PCG" else \
                    set(meta[(meta["group"] == "pcg")].index)
            offs = (sub[sub["gene_id"].isin(genes)]["start"] - 1 - 500).values
            offs = offs[(offs > -500) & (offs < 100)]
            if len(offs) < 50:
                continue
            # central enrichment: distances should cluster near TSS if functional
            centered = np.abs(offs) / 500.0
            stat, p = kstest(centered, "uniform")
            pb_rows.append(dict(cls=r["cls"], tf=r["tf"], n=len(offs),
                                mean_abs_offset=float(np.mean(np.abs(offs))),
                                ks_stat=stat, pval=p))
        pd.DataFrame(pb_rows).to_csv(
            safe_write_path(os.path.join(RES, "motif_position_bias.tsv")),
            sep="\t", index=False)

        # ---------------- real binding (ReMap) ----------------------------------
        remap_path = os.path.join(DATA, "remap_hits.tsv.gz")
        cand_rows = []
        if os.path.exists(remap_path):
            rh = pd.read_csv(remap_path, sep="\t",
                             names=["gene_id", "tf", "cell", "peak"], dtype=str)
            st.note(f"remap hits={len(rh):,} tfs={rh['tf'].nunique()}")
            rp = rh[["gene_id", "tf"]].drop_duplicates()
            rpiv = rp.assign(v=1).pivot_table(index="gene_id", columns="tf",
                                              values="v", fill_value=0)
            rows = []
            for tf in rpiv.columns:
                rows.extend([dict(tf=tf, **r) for r in
                             fisher_vs_pcg(rpiv[tf].astype(bool), meta, None)])
            be = pd.DataFrame(rows)
            be["fdr"] = bh_fdr(be["pval"].values)
            be.to_csv(safe_write_path(os.path.join(RES, "tf_binding_enrichment.tsv")),
                      sep="\t", index=False)
            st.note(f"binding rows={len(be)} sig={int((be['fdr'] < 0.05).sum())}")

            # candidates: motif AND binding concordant
            for cls in CLASSES:
                m1 = me[(me.cls == cls) & (me.fdr < 0.05) & (me.log2OR > 0.2)]
                b1 = be[(be.cls == cls) & (be.fdr < 0.05) & (be.log2OR > 0.2)]
                both = set(m1.tf) & set(b1.tf)
                monly = set(m1.tf) - both
                for tf in sorted(both):
                    cand_rows.append(dict(cls=cls, tf=tf, evidence="motif+binding",
                                          log2OR_motif=float(m1[m1.tf == tf].log2OR.iloc[0]),
                                          log2OR_binding=float(b1[b1.tf == tf].log2OR.iloc[0])))
                for tf in sorted(monly):
                    cand_rows.append(dict(cls=cls, tf=tf, evidence="motif_only",
                                          log2OR_motif=float(m1[m1.tf == tf].log2OR.iloc[0]),
                                          log2OR_binding=np.nan))
        cand = pd.DataFrame(cand_rows)
        cand.to_csv(safe_write_path(os.path.join(RES, "candidate_TFs.tsv")),
                    sep="\t", index=False)
        st.note(f"candidate TFs (rows)={len(cand)}; "
                f"dual-evidence={int((cand.evidence == 'motif+binding').sum())}")

        # summary
        top_me = (me[(me.fdr < 0.05) & (me.log2OR.abs() > 0.2)]
                  .sort_values(["cls", "log2OR"], ascending=[True, False])
                  .groupby("cls").head(15))
        body = ("## Motif enrichment (vs PCG, FDR<0.05, |log2OR|>0.2)\n\n"
                + (df_to_md(top_me, ".2f") if len(top_me) else "_(none significant)_")
                + "\n\n## Candidate TFs (motif AND ReMap binding)\n\n"
                + (df_to_md(cand[cand.evidence == "motif+binding"], ".2f")
                   if len(cand) else "_(none)_"))
        st.summary_md("results/step_summaries/08_enrichment.md",
                      "08 motif & binding enrichment", body)


if __name__ == "__main__":
    main()
