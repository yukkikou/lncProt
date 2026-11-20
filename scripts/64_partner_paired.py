#!/usr/bin/env python
"""64_partner_paired.py — divergent lncRNA versus its own partner, paired.

The Results compare divergent promoters with the genome-wide pool of
protein-coding promoters, and the gap is read as evidence that divergent
lncRNAs borrow strong initiation machinery. That comparison is confounded by
how the class is built: a lncRNA is called divergent only when a protein-coding
TSS is annotated within 1 kb, and well-annotated TSSs belong to well-expressed
genes. The comparison is therefore partly true by construction.

This module replaces the pool with each lncRNA's own partner. Every divergent
lncRNA is paired with the protein-coding gene that defines its class, and the
three promoter readouts are compared within the pair. A class that genuinely
borrows should resemble its partner. A class that merely sits next to strong
promoters will look weaker than its partner, because the lncRNA has no
independent initiation.

Readouts
  CAGE        is the promoter CAGE-active, McNemar on the discordant pairs
  chromatin   is the TSS in an active-promoter state, paired by cell type
  TF binding  fraction of MNase/ChIP-supported TFs bound, Wilcoxon on per-pair
              differences

Outputs
  results/partner_paired/cage.tsv
  results/partner_paired/chromatin.tsv
  results/partner_paired/binding.tsv
  results/step_summaries/64_partner_paired.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "partner_paired")
STATES = ["E003", "E017", "E067", "E071", "E088"]


def mcnemar(both, only_a, only_b, neither, label):
    """Exact McNemar on the two discordant cells."""
    n = only_a + only_b
    p = binomtest(only_a, n, 0.5).pvalue if n else np.nan
    return dict(readout=label, n_pairs=both + only_a + only_b + neither,
                both=both, lnc_only=only_a, partner_only=only_b,
                neither=neither,
                frac_lnc=float((both + only_a) / max(both + only_a + only_b
                                                     + neither, 1)),
                frac_partner=float((both + only_b) / max(both + only_a
                                                         + only_b + neither, 1)),
                pval=float(p) if n else np.nan)


def main():
    with step("64_partner_paired") as st:
        os.makedirs(OUT, exist_ok=True)
        c = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                        low_memory=False)
        c = c.set_index("gene_id")
        div = c[(c.cls == "divergent") & c.partner_pcg.notna()].copy()
        div = div[div.partner_pcg.isin(c.index)]
        st.note("divergent lncRNAs with an annotated partner in the table: "
                "{}".format(len(div)))

        # ---- CAGE
        lt = div.cage_refined.fillna(0).astype(int).values
        pt = c.loc[div.partner_pcg, "cage_refined"].fillna(0).astype(
            int).values
        row = mcnemar(int(((lt == 1) & (pt == 1)).sum()),
                      int(((lt == 1) & (pt == 0)).sum()),
                      int(((lt == 0) & (pt == 1)).sum()),
                      int(((lt == 0) & (pt == 0)).sum()),
                      "CAGE-active promoter")
        T1 = pd.DataFrame([row])
        T1.to_csv(safe_write_path(os.path.join(OUT, "cage.tsv")), sep="\t",
                  index=False)
        st.note("CAGE: lnc active {:.1%}, partner active {:.1%}, "
                "partner-only {} lnc-only {}, P={:.3g}".format(
                    row["frac_lnc"], row["frac_partner"], row["partner_only"],
                    row["lnc_only"], row["pval"]))

        # ---- chromatin state, paired within cell type
        ch = pd.read_csv(os.path.join(RES, "chromhmm_states.tsv"), sep="\t",
                         low_memory=False).set_index("gene_id")
        ch = ch.reindex(c.index)
        rows, cell_rows = [], []
        for e in STATES:
            if e not in ch.columns:
                continue
            s_lnc = ch.loc[div.index, e].astype(str)
            s_par = ch.loc[div.partner_pcg, e].astype(str).values
            s_lnc = s_lnc.values
            # Roadmap vocabulary. Active promoter means TssA or its flanking
            # state, the same definition the class-level table uses.
            act_lnc = np.isin(s_lnc, ["TssA", "TssAFlnk"])
            act_par = np.isin(s_par, ["TssA", "TssAFlnk"])
            m = (s_lnc != "nan") & (s_par != "nan")
            cell_rows.append(mcnemar(
                int((act_lnc & act_par & m).sum()),
                int((act_lnc & ~act_par & m).sum()),
                int((~act_lnc & act_par & m).sum()),
                int((~act_lnc & ~act_par & m).sum()), e))
        T2 = pd.DataFrame(cell_rows)
        if len(T2):
            T2.to_csv(safe_write_path(os.path.join(OUT, "chromatin.tsv")),
                      sep="\t", index=False)
            st.note("chromatin: lnc in a promoter state {:.1%} versus partner "
                    "{:.1%} across {} cell types".format(
                        T2.frac_lnc.mean(), T2.frac_partner.mean(), len(T2)))

        # ---- TF binding, paired per TF
        mp = pd.read_csv(os.path.join(RES, "motif_presence.tsv"), sep="\t",
                         index_col=0, low_memory=False)
        mp = mp.reindex(c.index)
        sub_lnc = mp.reindex(div.index)
        sub_par = mp.reindex(div.partner_pcg)
        n_tf = sub_lnc.notna().any(axis=0).sum()
        keep = [t for t in mp.columns
                if sub_lnc[t].notna().sum() > 100 and sub_par[t].notna().sum() > 100]
        rows = []
        for t in keep:
            a = sub_lnc[t].fillna(0).values
            b = sub_par[t].fillna(0).values
            rows.append(dict(tf=t, frac_lnc=float(a.mean()),
                             frac_partner=float(b.mean()),
                             diff=float(a.mean() - b.mean())))
        T3 = pd.DataFrame(rows).sort_values("diff")
        if len(T3):
            T3.to_csv(safe_write_path(os.path.join(OUT, "binding.tsv")),
                      sep="\t", index=False)
            w, pw = wilcoxon(T3.frac_lnc - T3.frac_partner)
            st.note("TF presence: {} TFs, median lnc {:.3f} versus partner "
                    "{:.3f}, paired Wilcoxon P={:.3g}".format(
                        len(T3), T3.frac_lnc.median(), T3.frac_partner.median(),
                        pw))

        # ---- ChIP-supported TF binding, from the ReMap hit list
        hits = {}
        want = set(div.index) | set(div.partner_pcg)
        path = os.path.join(WORKDIR, "data", "remap_hits.tsv.gz")
        if os.path.exists(path):
            for chunk in pd.read_csv(path, sep="\t", header=None,
                                     names=["gene", "tf", "exp", "acc"],
                                     usecols=[0, 1], chunksize=5_000_000):
                sub = chunk[chunk.gene.isin(want)]
                for g, t in zip(sub.gene.values, sub.tf.values):
                    hits.setdefault(g, set()).add(t)
            st.note("ReMap: {} of {} target genes carry a binding site".format(
                len(hits), len(want)))
            a = np.array([len(hits.get(g, ())) for g in div.index])
            b = np.array([len(hits.get(g, ())) for g in div.partner_pcg])
            ntf = max(len({t for s in hits.values() for t in s}), 1)
            T4 = pd.DataFrame([dict(
                readout="ChIP TF count per promoter",
                n_pairs=len(a), lnc_median=float(np.median(a)),
                partner_median=float(np.median(b)),
                lnc_mean=float(a.mean()), partner_mean=float(b.mean()),
                pval=float(wilcoxon(a, b).pvalue),
                n_tf_universe=int(ntf))])
            T4.to_csv(safe_write_path(os.path.join(OUT, "binding_chip.tsv")),
                      sep="\t", index=False)
            st.note("ChIP: lnc median {} TFs versus partner {}, mean {:.2f} "
                    "versus {:.2f}, P={:.3g}".format(
                        np.median(a), np.median(b), a.mean(), b.mean(),
                        T4.pval.iloc[0]))
        else:
            T4 = pd.DataFrame()

        st.summary_md(os.path.join(RES, "step_summaries",
                                   "64_partner_paired.md"),
                      "64 divergent lncRNA versus its own partner",
                      df_to_md(T1, ".3f") + "\n\n"
                      + (df_to_md(T2, ".3f") if len(T2) else "")
                      + "\n\nTF presence, 10 largest gaps\n\n"
                      + (df_to_md(pd.concat([T3.head(5), T3.tail(5)]), ".3f")
                         if len(T3) else ""))


if __name__ == "__main__":
    main()
