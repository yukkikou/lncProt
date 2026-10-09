#!/usr/bin/env python
"""31_tf_contexts.py — TF requirements per birth class + tissue/cancer context.

Modes:
  A (local)  - per-class top-TF pairwise motif co-occurrence (Fisher OR vs
               PCG background) + cofactor binding occupancy per class
               (EP300/SMARCA4/etc from remap_hits pivot).
  T (matrix) - GTEx tissues: Spearman(TF expression across 49 tissues,
               class activity) for a candidate-TF panel; TCGA pan-cancer:
               per-sample TF expression vs LTR-lnc score and class activity,
               aggregated pan-cancer. One npz load.
Outputs -> results/tf_contexts/*.tsv + step summaries.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
TD = os.path.join(WORKDIR, "data", "tissue_disease")
OUT = os.path.join(RES, "tf_contexts")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
COFACTORS = ["EP300", "CREBBP", "SMARCA4", "ARID1A", "MED1", "BRD4",
             "TAF1", "YY1", "CTCF"]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def mode_a():
    with step("31_tf_A") as st:
        mp = pd.read_csv(os.path.join(RES, "motif_presence.tsv"), sep="\t",
                         index_col=0, low_memory=False)
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        mp = mp.join(meta[["cls", "group"]])
        me = pd.read_csv(os.path.join(RES, "motif_enrichment.tsv"), sep="\t")
        top = (me.sort_values("log2OR", ascending=False)
                 .groupby("cls").head(25))
        pcg = mp[mp.group == "pcg"]
        rows = []
        for cls in CLASSES:
            sub = mp[mp.cls == cls]
            tfs = [t for t in top[top.cls == cls].tf.unique()
                   if t in mp.columns][:20]
            for i in range(len(tfs)):
                for j in range(i + 1, len(tfs)):
                    a, b = tfs[i], tfs[j]
                    x = (sub[a] == 1) & (sub[b] == 1)
                    pa = (sub[a] == 1)
                    pb = (sub[b] == 1)
                    # co-occurrence OR within class vs its own marginals
                    t11 = int((pa & pb).sum()); t10 = int((pa & ~pb).sum())
                    t01 = int((~pa & pb).sum())
                    t00 = int((~pa & ~pb).sum())
                    orr, p = fisher_exact([[t11, t10], [t01, t00]])
                    rows.append(dict(cls=cls, tf_a=a, tf_b=b, n_class=len(sub),
                                     both=t11, orr=float(np.log2(orr))
                                     if orr > 0 else np.nan, pval=p))
        P = pd.DataFrame(rows)
        P["fdr"] = bh_fdr(P.pval.values)
        P = P.sort_values(["cls", "orr"], ascending=[True, False])
        P.to_csv(safe_write_path(os.path.join(OUT, "A_tf_pairs.tsv")),
                 sep="\t", index=False)

        # cofactor occupancy per class from remap_hits (gene x tf presence)
        hits = pd.read_csv(os.path.join(WORKDIR, "data",
                                        "remap_hits.tsv.gz"), sep="\t",
                           header=None, names=["gene_id", "tf", "cell",
                                               "study"],
                           usecols=["gene_id", "tf"])
        hits = hits[hits.tf.isin(COFACTORS)]
        binm = hits.drop_duplicates(["gene_id", "tf"]).assign(v=1) \
            .pivot(index="gene_id", columns="tf", values="v").fillna(0)
        m = meta[["cls", "group"]].join(binm, how="left").fillna(0)
        crows = []
        for cf in [c for c in COFACTORS if c in m.columns]:
            for cls in CLASSES + ["PCG"]:
                v = m.loc[(m.cls == cls) if cls != "PCG"
                          else (m.group == "pcg"), cf]
                crows.append(dict(cls=cls, cofactor=cf, frac_bound=float(v.mean()),
                                  n=len(v)))
        C = pd.DataFrame(crows)
        C.to_csv(safe_write_path(
            os.path.join(OUT, "A_cofactor_occupancy.tsv")), sep="\t",
            index=False)
        st.note("A: pairs={} (top: {}); cofactor divergent EP300={:.2f} "
                "intergenic EP300={:.2f}".format(
                    len(P),
                    "; ".join("{}:{}+{}".format(r.cls, r.tf_a, r.tf_b)
                              for r in P.head(3).itertuples()),
                    C[(C.cls == "divergent") &
                      (C.cofactor == "EP300")].frac_bound.iloc[0],
                    C[(C.cls == "intergenic") &
                      (C.cofactor == "EP300")].frac_bound.iloc[0]))
        st.summary_md("results/tf_contexts/step_summaries/31_a.md",
                      "31A TF pairs + cofactor occupancy",
                      "## Top motif co-occurrence pairs per class\n\n"
                      + df_to_md(P.groupby("cls").head(10), ".3f")
                      + "\n\n## Cofactor occupancy per class\n\n"
                      + df_to_md(C.pivot(index="cls", columns="cofactor",
                                         values="frac_bound"), ".2f"))


def mode_t():
    with step("31_tf_T") as st:
        d = np.nan_to_num(np.load(os.path.join(TD, "toil_class_expr.npz"))
                          ["M"])
        genes = pd.read_csv(os.path.join(TD, "toil_genes.tsv"), sep="\t",
                            index_col=0)
        samp = pd.read_csv(os.path.join(TD, "toil_samples.tsv"), sep="\t")
        gsym = genes["symbol"]
        cand = pd.read_csv(os.path.join(RES, "candidate_TFs.tsv"), sep="\t")
        panel = sorted(set(cand.tf.unique())
                       | {"DUX4", "PHOX2B", "BATF3", "CTCF", "MAZ"})
        pos_of = {g: k for k, g in enumerate(genes.index)}
        tf_idx = {t: pos_of[g] for g, t in zip(genes.index, gsym)
                  if t in panel}
        st.note("TF panel mapped: {}/{}".format(len(tf_idx), len(panel)))

        # ---- GTEx tissues: TF expr x class activity across 49 tissues -----
        gtex = samp[samp.cohort == "GTEx"]
        tissues = sorted(t for t, n in gtex.tissue.value_counts().items()
                         if n >= 20)
        gcols = gtex.index.to_numpy()
        W = np.zeros((len(tissues), len(gtex)), dtype=np.float32)
        for k, t in enumerate(tissues):
            m = (gtex.tissue == t).values
            W[k, m] = 1.0
        W /= W.sum(axis=1, keepdims=True)
        cls_mask = {c: (genes.cls == c).values for c in CLASSES}
        need_rows = np.array(sorted(set(tf_idx.values())
                                    | set(np.where(genes.cls.isin(
                                        CLASSES).values)[0].tolist())))
        sub = d[need_rows][:, gcols] @ W.T          # genes x tissues
        pos = {g: k for k, g in enumerate(need_rows)}
        tf_expr = pd.DataFrame({t: sub[pos[i]] for t, i in tf_idx.items()},
                               index=tissues)
        cls_act = pd.DataFrame(
            {c: (sub[[pos[i] for i in np.where(cls_mask[c])[0]]]
                 >= 1).mean(axis=0) for c in CLASSES}, index=tissues)
        rows = []
        for t in tf_expr.columns:
            for c in CLASSES:
                r, p = spearmanr(tf_expr[t], cls_act[c])
                rows.append(dict(scope="gtex_tissue", tf=t, cls=c, rho=r,
                                 pval=p))
        T1 = pd.DataFrame(rows)
        T1["fdr"] = bh_fdr(T1.pval.values)

        # ---- TCGA pan-cancer: per-sample TF expr vs LTR score -------------
        tcga = samp[(samp.cohort == "TCGA") & (samp.condition == "tumor")]
        ltr_rows = genes.cls.isin(CLASSES).values
        ltr = ltr_rows & (genes.tss_in_ltr == 1).fillna(False).values
        oth = ltr_rows & ~ltr
        ltr_score = d[ltr].mean(axis=0) - d[oth].mean(axis=0)
        score_s = pd.Series(ltr_score, index=samp.index)
        rows = []
        key_tfs = sorted(set(list(tf_idx)[:60]))
        tf_vec = {t: pd.Series(d[i], index=samp.index)
                  for t, i in list(tf_idx.items())[:60]}
        for cancer, subg in tcga.groupby("cancer"):
            if len(subg) < 30:
                continue
            sc = score_s.reindex(subg.index)
            for t in key_tfs:
                r, p = spearmanr(tf_vec[t].reindex(subg.index), sc)
                rows.append(dict(scope="tcga_" + cancer, tf=t, cls="LTR_score",
                                 rho=r, pval=p, n=len(subg)))
        T2 = pd.DataFrame(rows)
        T2["fdr"] = bh_fdr(T2.pval.values)
        T = pd.concat([T1, T2], ignore_index=True)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "T_tf_class_coupling.tsv")), sep="\t", index=False)
        top = T1[(T1.fdr < 0.05)].sort_values("rho", ascending=False)
        st.note("T: gtex couplings sig={}; top: ".format(len(top))
                + "; ".join("{} {} rho={:.2f}".format(r.tf, r.cls, r.rho)
                            for r in top.head(5).itertuples()))
        st.summary_md("results/tf_contexts/step_summaries/31_t.md",
                      "31T TF-class coupling (tissues + pan-cancer)",
                      "## GTEx tissue couplings (top by rho)\n\n"
                      + df_to_md(top.head(30), ".3f")
                      + "\n\n## Pan-cancer LTR-score couplings (top |rho|)\n\n"
                      + df_to_md(T2.reindex(T2.rho.abs().sort_values(
                          ascending=False).index).head(30), ".3f"))


if __name__ == "__main__":
    mode = (sys.argv[1] if len(sys.argv) > 1 else "").upper()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    if mode == "A":
        mode_a()
    elif mode == "T":
        mode_t()
    else:
        raise SystemExit("usage: 31_tf_contexts.py A|T")
