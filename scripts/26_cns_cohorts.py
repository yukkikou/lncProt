#!/usr/bin/env python
"""26_cns_cohorts.py — extra CNS tumour cohorts (MB, SFT; PitNET precomputed).

MB   1610305236_CJ_MB: 206 salmon libs, heterogeneous stems (S-lib, YY-NNNNNN-lib,
     plain numeric); libs aggregated to samples by stem core regex.
SFT  1610305236_CJ_SFT/RNA-seq: 107 libs; grade crosswalk via 6-digit stem
     prefix -> SFT_WGS_Sample_Info.tsv Sample_ID (37/107 resolved).
PitNET: analysed from the project's merged TPM matrix -> results/pitnet/.

Per cohort: class-detection hierarchy, LTR-relative activation; SFT grade
trend (Spearman + MWU III vs I). Outputs -> results/cns_cohorts/.
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR

QD = {"mb": os.path.join(WORKDIR, "data", "mb_own", "quants"),
      "sft": os.path.join(WORKDIR, "data", "sft_own", "quants")}
OUT = os.path.join(WORKDIR, "results", "cns_cohorts")
RES = os.path.join(WORKDIR, "results")
GTF = ("/media/share/resource/reference/homo_sapiens/GRCh38/annotation/"
       "gencode.v46.annotation.gtf")
SFT_INFO = ("/media/share/data8/volume1/1610305236_CJ_SFT/"
            "SFT_WGS_Sample_Info.tsv")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def tx2gene():
    import gzip
    tx = {}
    pat = re.compile(r'transcript_id "([^"]+)";.*?gene_id "([^"]+)"')
    pat2 = re.compile(r'gene_id "([^"]+)";.*?transcript_id "([^"]+)"')
    op = gzip.open if GTF.endswith(".gz") else open
    with op(GTF, "rt") as fh:
        for line in fh:
            if not line.startswith("#") and "\ttranscript\t" in line:
                m = pat.search(line) or pat2.search(line)
                if m:
                    a, b = m.group(1), m.group(2)
                    if a.startswith("ENST"):
                        tx[a] = b
                    else:
                        tx[b] = a
    return tx


def sample_core(stem):
    m = re.match(r"S?(\d{2}-\d{6})", stem)
    return m.group(1) if m else stem


def per_sample_table(meta, te, quants, t2g, aggregate=True):
    sym2cls = dict(zip(meta.gene_name.astype(str), meta.cls))
    g2sym = {g.split(".")[0]: str(s) for g, s in
             zip(meta.gene_id, meta.gene_name)}
    sym2ltr = {str(s): int(v) for s, v in zip(
        te.gene_name, pd.to_numeric(te.tss_in_ltr, errors="coerce")
        .fillna(0)) if pd.notna(s)}
    cols = {}
    for q in sorted(quants):
        s = os.path.basename(os.path.dirname(q))
        d = pd.read_csv(q, sep="\t", usecols=["Name", "TPM"])
        d["gene"] = d.Name.map(lambda t: t2g.get(t)
                               or t2g.get(t.split(".")[0]))
        g = d.dropna(subset=["gene"]).groupby("gene")["TPM"].sum()
        cols[s] = g
    M = pd.DataFrame(cols).fillna(0)
    M.index = [g.split(".")[0] for g in M.index]
    M = M.groupby(M.index).sum()
    if aggregate:
        core = {s: sample_core(s) for s in M.columns}
        M = M.groupby(core, axis=1).mean()
    cls = np.array([sym2cls.get(g2sym.get(g, ""), None) for g in M.index])
    ltr = np.array([sym2ltr.get(g2sym.get(g, ""), 0) for g in M.index])
    LM = np.log2(M + 1)
    lnc = np.isin(cls, CLASSES)
    rows = []
    for s in M.columns:
        v = LM[s]
        r = dict(sample=s)
        for c in CLASSES + ["PCG"]:
            idx = np.isin(cls, c if c != "PCG" else ["PCG"])
            r["{}_frac".format(c)] = float((v[idx] >= 1).sum()
                                           / max(int(idx.sum()), 1))
        r["ltr_rel"] = float(v[lnc & (ltr == 1)].mean()
                             - v[lnc & (ltr == 0)].mean())
        rows.append(r)
    T = pd.DataFrame(rows)
    return T, M


def main():
    with step("26_cns_cohorts") as st:
        os.makedirs(OUT, exist_ok=True)
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", low_memory=False)
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                         sep="\t", low_memory=False)
        te["gene_name"] = te.gene_name.astype(str)
        t2g = tx2gene()

        # ---- MB ----------------------------------------------------------
        qm = glob.glob(os.path.join(QD["mb"], "*", "quant.sf"))
        st.note("MB quants: {}".format(len(qm)))
        Tm, Mm = per_sample_table(meta, te, qm, t2g, aggregate=True)
        Tm.to_csv(safe_write_path(os.path.join(OUT, "mb_per_sample.tsv")),
                  sep="\t", index=False)
        hier = ", ".join("{} {:.0f}%".format(c, Tm["{}_frac".format(c)]
                                             .mean() * 100)
                         for c in ["PCG"] + CLASSES)
        st.note("MB samples={} hierarchy: {}".format(len(Tm), hier))
        st.note("MB LTR-delta median {:+.3f} (IQR {:+.3f}..{:+.3f})".format(
            Tm.ltr_rel.median(), Tm.ltr_rel.quantile(.25),
            Tm.ltr_rel.quantile(.75)))

        # ---- SFT ----------------------------------------------------------
        qs = glob.glob(os.path.join(QD["sft"], "*", "quant.sf"))
        st.note("SFT quants: {}".format(len(qs)))
        n_grade = 0
        Ts, Ms = per_sample_table(meta, te, qs, t2g, aggregate=False)
        info = pd.read_csv(SFT_INFO, sep="\t")
        info["sid"] = info.Sample_ID.astype(str)
        Ts["sid"] = Ts["sample"].str.split("-").str[0]
        Ts = Ts.merge(info[["sid", "Grade", "Sample_Name"]], on="sid",
                      how="left")
        Ts.to_csv(safe_write_path(os.path.join(OUT, "sft_per_sample.tsv")),
                  sep="\t", index=False)
        hier = ", ".join("{} {:.0f}%".format(c, Ts["{}_frac".format(c)]
                                             .mean() * 100)
                         for c in ["PCG"] + CLASSES)
        st.note("SFT samples={} hierarchy: {}".format(len(Ts), hier))
        st.note("SFT LTR-delta median {:+.3f}".format(Ts.ltr_rel.median()))
        g = pd.to_numeric(Ts.Grade, errors="coerce")
        ok = Ts.ltr_rel.notna() & g.notna() & g.isin([1, 2, 3])
        if int(ok.sum()) >= 10:
            for gv in (1, 2, 3):
                sub = Ts.ltr_rel[ok & (g == gv)]
                st.note("SFT grade {}: LTR-delta {:+.3f} (n={})".format(
                    gv, sub.median(), len(sub)))
            lo, hi = Ts.ltr_rel[ok & (g == 1)], Ts.ltr_rel[ok & (g == 3)]
            rho, pr = spearmanr(g[ok], Ts.ltr_rel[ok])
            _, p = mannwhitneyu(hi, lo)
            n_grade = int(ok.sum())
            st.note("SFT III vs I MWU p={:.2e}; trend rho={:.2f} "
                    "p={:.1e} (n={})".format(p, rho, pr, n_grade))

        st.summary_md(
            "results/cns_cohorts/step_summaries/26_cns_cohorts.md",
            "26 Extra CNS cohorts (MB / SFT; PitNET separately)",
            "- MB (medulloblastoma): {} libs -> {} samples; hierarchy {}; "
            "LTR-delta median {:+.3f}\n"
            "- SFT: {} samples; grade-resolved {}/{}; see notes above\n"
            "- PitNET: results/pitnet/pitnet_per_sample.tsv (207 tumours vs "
            "10 normals, MWU p=3.9e-02)".format(
                len(qm), len(Tm),
                ", ".join("{} {:.0f}%".format(c, Tm["{}_frac".format(c)]
                                              .mean() * 100)
                          for c in ["PCG", "divergent", "antisense",
                                    "intergenic"]),
                Tm.ltr_rel.median(), len(Ts), n_grade, len(Ts)))


if __name__ == "__main__":
    main()
