#!/usr/bin/env python
"""78_anage_lifespan.py — can AnAge support the cancer-conservation axis?

The Nature Communications paper on longer-lived mammals compares per-species
evolutionary rates at human cancer-associated positions. Reproducing that needs
per-species alignments or precomputed branch rates. Neither is on this cluster:
the only conservation resource is an aggregate 241-mammal phyloP track, one
score per base rather than one per species.

This module records what is possible with the seven species that have syntenic
lncRNA expression data. For each species it recomputes the human-versus-species
developmental breadth correlation, joins the AnAge maximum lifespan, and tests
whether breadth conservation rises with lifespan. This is an honest power check,
not a substitute for the cross-species rate analysis.

Outputs
  results/comparative/anage_species_panel.tsv
  results/step_summaries/78_anage_lifespan.md
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
SYN = os.path.join(RES, "syntenic")
ANAGE = os.path.join(WORKDIR, "data", "anage", "anage_data.txt")
OUT = os.path.join(RES, "comparative")
SPECIES = {"mouse": "Mus musculus", "macaque": "Macaca mulatta",
           "rat": "Rattus norvegicus", "rabbit": "Oryctolagus cuniculus",
           "chicken": "Gallus gallus", "opossum": "Monodelphis domestica"}


def main():
    with step("78_anage_lifespan") as st:
        aa = pd.read_csv(ANAGE, sep="\t", low_memory=False)
        lcol = [c for c in aa.columns if "Maximum longevity" in c][0]
        aa = aa[["Genus", "Species", lcol]].dropna()
        aa["binomial"] = aa.Genus + " " + aa.Species
        life = dict(zip(aa.binomial, aa[lcol]))
        st.note("AnAge rows with maximum longevity={}".format(len(aa)))

        rows = []
        for sp, name in SPECIES.items():
            f = os.path.join(SYN, "map_{}.tsv".format(sp))
            if not os.path.exists(f):
                continue
            d = pd.read_csv(f, sep="\t").dropna(
                subset=["human_breadth", "species_breadth"])
            if len(d) < 10:
                continue
            r, p = spearmanr(d.human_breadth, d.species_breadth)
            rows.append(dict(species=sp, binomial=name,
                             max_lifespan_years=life.get(name, np.nan),
                             n_genes=len(d), rho_breadth=r, pval=p))
        T = pd.DataFrame(rows)

        tests = []
        for lab, sub in (("all species", T),
                         ("n>=40", T[T.n_genes >= 40])):
            sub = sub.dropna(subset=["max_lifespan_years", "rho_breadth"])
            if len(sub) >= 4:
                r, p = spearmanr(sub.max_lifespan_years, sub.rho_breadth)
                tests.append(dict(subset=lab, n_species=len(sub),
                                  rho=float(r), pval=float(p)))
        U = pd.DataFrame(tests)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "anage_species_panel.tsv")), sep="\t", index=False)

        for r in T.itertuples():
            st.note("{}: lifespan={:.1f} y, n={}, breadth rho={:.3f}".format(
                r.species, r.max_lifespan_years, r.n_genes, r.rho_breadth))
        for r in U.itertuples():
            st.note("{}: {} species, Spearman rho={:.3f}, P={:.3f}".format(
                r.subset, r.n_species, r.rho, r.pval))
        st.note("Per-species evolutionary rates are not available locally, so "
                "the 58-species cancer-locus analysis of the reference paper "
                "cannot be reproduced.")
        st.summary_md(os.path.join(
            RES, "step_summaries", "78_anage_lifespan.md"),
            "78 AnAge lifespan panel",
            df_to_md(T, ".3f") + "\n\n" + df_to_md(U, ".3f"))


if __name__ == "__main__":
    main()
