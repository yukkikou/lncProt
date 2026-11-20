#!/usr/bin/env python
"""38_developmental.py — broad developmental-expression analysis across
Cardoso-Moreira 2019 (7 species, Expression Atlas).

Human series (E-MTAB-6814): expression breadth / organ breadth / stage
breadth / tau by birth class, and association with transcript age.
Cross-species: for one-to-one symbol orthologs, conservation of expression
breadth between human and mouse/rat/macaque/rabbit/chicken/opossum, by class.
Outputs -> results/developmental/*.tsv + step summary.
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
EA = os.path.join(WORKDIR, "data", "expression_atlas")
OUT = os.path.join(RES, "developmental")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
SPECIES = [("human", "E-MTAB-6814"), ("mouse", "E-MTAB-6798"),
           ("rat", "E-MTAB-6811"), ("macaque", "E-MTAB-6813"),
           ("rabbit", "E-MTAB-6782"), ("chicken", "E-MTAB-6769"),
           ("opossum", "E-MTAB-6833")]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def mwu(y, x):
    u, p = mannwhitneyu(y, x, alternative="two-sided")
    return float(2 * u / (len(y) * len(x)) - 1), float(p)


def tau_spec(v):
    v = np.asarray(v, float)
    if v.max() <= 0 or len(v) < 2:
        return np.nan
    x = v / v.max()
    return float((1 - x).sum() / (len(v) - 1))


def parse_config(exp):
    p = os.path.join(EA, exp, exp + "-configuration.xml")
    txt = open(p, encoding="utf-8", errors="ignore").read()
    return {m.group(1): m.group(2) for m in
            re.finditer(r'id="(g\d+)"\s+label="([^"]*)"', txt)}


def load_species(exp):
    tp = pd.read_csv(os.path.join(EA, exp, exp + "-tpms.tsv"), sep="\t",
                     low_memory=False)
    gid = tp.columns[0]
    gname = tp.columns[1]
    cfg = parse_config(exp)
    samples = [c for c in tp.columns[2:] if c in cfg]
    if not samples:
        samples = list(tp.columns[2:])
    def col_mean(s):
        parts = s.astype(str).str.split(",", expand=True).astype(float)
        return parts.mean(axis=1)
    X = pd.DataFrame({c: col_mean(tp[c]) for c in samples})
    meta = pd.DataFrame({"gid": tp[gid].astype(str).values,
                         "symbol": tp[gname].astype(str).values})
    labels = [cfg.get(s, "") for s in samples]
    stage = [l.split(";")[0].strip() if l else "" for l in labels]
    organ = [l.split(";")[1].strip().lower() if ";" in l else ""
             for l in labels]
    summary = meta.copy()
    summary["n_samples"] = X.shape[1]
    summary["max_tpm"] = X.max(axis=1).values
    summary["mean_tpm"] = X.mean(axis=1).values
    summary["breadth_frac"] = (X > 1).mean(axis=1).values
    summary["tau"] = X.apply(tau_spec, axis=1).values
    # organ / stage means
    org = pd.Series(organ)
    stg = pd.Series(stage)
    organ_means = X.T.groupby(org.values).mean().T
    stage_means = X.T.groupby(stg.values).mean().T
    summary["n_organs"] = int(org[org != ""].nunique())
    summary["organ_breadth"] = (organ_means > 1).sum(axis=1).values
    summary["stage_breadth"] = (stage_means > 1).sum(axis=1).values
    summary["n_stages"] = int(stg[stg != ""].nunique())
    return summary, organ_means, stage, organ


def main():
    with step("38_developmental") as st:
        summaries = {}
        organ_mats = {}
        for sp, exp in SPECIES:
            s, om, stage, organ = load_species(exp)
            summaries[sp] = s
            organ_mats[sp] = om
            st.note("{}: genes={} samples={} organs={} stages={}".format(
                sp, len(s), int(s.n_samples.iloc[0]),
                len(set(o for o in organ if o)),
                len(set(x for x in stage if x))))
        allsum = pd.concat([v.assign(species=k) for k, v in
                            summaries.items()], ignore_index=True)
        allsum = allsum.rename(columns={"gid": "gene"})
        allsum.to_csv(safe_write_path(os.path.join(
            OUT, "per_gene_summary.tsv")), sep="\t", index=False)

        # ---- human class analysis -------------------------------------
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        base = pd.Series(meta.index.str.split(".").str[0].values,
                         index=meta.index)
        m2 = meta.copy()
        m2["base"] = base.values
        h = summaries["human"].copy()
        h["base"] = h.gid.str.split(".").str[0]
        h = h.merge(m2[["base", "cls", "group"]].drop_duplicates("base"),
                    on="base", how="left")
        h.to_csv(safe_write_path(os.path.join(OUT, "human_by_gene.tsv")),
                 sep="\t", index=False)
        rows = []
        pcg = h[h.group == "pcg"]
        for feat in ["breadth_frac", "organ_breadth", "stage_breadth",
                     "tau", "max_tpm"]:
            x = pcg[feat].dropna().values
            for cls in CLASSES + ["pseudogene", "transcribed_pseudogene"]:
                y = h.loc[h.cls == cls, feat].dropna().values
                if len(y) < 20:
                    continue
                d, p = mwu(y, x)
                rows.append(dict(feature=feat, cls=cls, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=d, pval=p))
        T = pd.DataFrame(rows)
        T["fdr"] = bh_fdr(T.pval.values)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "human_class_tests.tsv")), sep="\t", index=False)

        # age association
        age = pd.read_csv(os.path.join(RES, "comparative",
                                       "E_transcript_age.tsv"), sep="\t",
                          index_col=0)[["age_level"]]
        age.index = age.index.str.split(".").str[0]
        age = age[~age.index.duplicated()]
        hh = h.set_index("base").join(age, how="left")
        arows = []
        for feat in ["breadth_frac", "organ_breadth", "stage_breadth", "tau"]:
            ok = hh[["age_level", feat]].dropna()
            if len(ok) > 100:
                r, p = spearmanr(ok.age_level, ok[feat])
                arows.append(dict(feature=feat, n=len(ok), rho=float(r),
                                  pval=float(p)))
        A = pd.DataFrame(arows)
        A.to_csv(safe_write_path(os.path.join(OUT, "age_vs_breadth.tsv")),
                 sep="\t", index=False)

        # ---- cross-species breadth conservation ----------------------
        hs = summaries["human"].copy()
        hs["sym"] = hs.symbol.str.upper()
        hs = hs.dropna(subset=["sym"]).drop_duplicates("sym")
        hmeta = m2[["base", "cls", "group"]].drop_duplicates("base")
        hs["base"] = hs.gid.str.split(".").str[0]
        hs = hs.merge(hmeta, on="base", how="left")
        rows = []
        for sp, exp in SPECIES:
            if sp == "human":
                continue
            ss = summaries[sp].copy()
            ss["sym"] = ss.symbol.str.upper()
            ss = ss.dropna(subset=["sym"]).drop_duplicates("sym")
            mg = hs.merge(ss[["sym", "breadth_frac", "max_tpm", "tau"]],
                          on="sym", suffixes=("_human", "_sp"))
            for cls in ["PCG"] + CLASSES:
                sub = mg[(mg.cls == cls) if cls != "PCG" else
                         (mg.group == "pcg")]
                if len(sub) < 20:
                    continue
                r1, p1 = spearmanr(sub.breadth_frac_human,
                                   sub.breadth_frac_sp)
                r2, p2 = spearmanr(sub.max_tpm_human, sub.max_tpm_sp)
                rows.append(dict(species=sp, cls=cls, n=len(sub),
                                 rho_breadth=float(r1), p_breadth=float(p1),
                                 rho_max=float(r2), p_max=float(p2)))
        X = pd.DataFrame(rows)
        X.to_csv(safe_write_path(os.path.join(
            OUT, "cross_species.tsv")), sep="\t", index=False)

        st.note("human class tests: " + "; ".join(
            "{} {} d={:.3f} p={:.1e}".format(
                r.feature, r.cls, r.cliffs_delta, r.pval)
            for r in T.itertuples() if r.cls in ("divergent", "intergenic")
            and r.feature == "breadth_frac"))
        st.note("cross-species rho(breadth): " + "; ".join(
            "{} {} {:.2f}".format(r.species, r.cls, r.rho_breadth)
            for r in X.itertuples() if r.cls in ("PCG", "divergent",
                                                 "intergenic")
            and r.species in ("mouse", "macaque")))
        st.summary_md("results/developmental/step_summaries/38_dev.md",
                      "38 developmental expression (7 species)",
                      "## Human class vs PCG\n\n" + df_to_md(T, ".3f")
                      + "\n\n## Age vs developmental breadth\n\n"
                      + df_to_md(A, ".3f")
                      + "\n\n## Cross-species breadth conservation\n\n"
                      + df_to_md(X, ".3f", max_rows=60))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
