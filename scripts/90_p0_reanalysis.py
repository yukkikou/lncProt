#!/usr/bin/env python
"""90_p0_reanalysis.py — reviewer P0 reanalyses (round 3).

Modes
  p2  orientation-resolved classes and headline statistics
  p1  separation of TSS placement from activity measurement
  p3  window-spillover controls for RAMPAGE and GRO-cap
  p4  age identifiability: mouse-primary estimate and permutation null

Outputs under results/review_p0/.
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "review_p0")
CLASSES = ["divergent_h2h", "antisense", "intronic", "sense_overlapping",
           "intergenic", "pseudogene", "transcribed_pseudogene", "PCG"]


def load_master():
    m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                    index_col=0, low_memory=False)
    m["tss"] = pd.to_numeric(m["tss"], errors="coerce")
    return m


def genes_meta():
    return pd.read_csv(os.path.join(WORKDIR, "data", "genes", "genes.tsv"),
                       sep="\t").set_index("gene_id")


def is_h2h(m, lnc_gid, pcg_gid):
    ps = m.strand.get(pcg_gid)
    lt = m.tss.get(lnc_gid)
    pt = m.tss.get(pcg_gid)
    if pd.isna(lt) or pd.isna(pt):
        return None
    return (lt < pt) if ps == "+" else (lt > pt)


def reassign(m, genes):
    """Return a Series mapping each lncRNA gene to its revised class."""
    exon = pd.read_csv(os.path.join(RES, "exon_overlap_pairs.tsv"), sep="\t")
    ex = set(zip(exon.lnc_gid, exon.pcg_gid))
    cls = m.cls.copy()
    dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
    for r in dp.itertuples():
        h = is_h2h(m, r.lnc_gid, r.pcg_gid)
        if h is None:
            continue
        if h:
            cls[r.lnc_gid] = "divergent_h2h"
            continue
        # convergent: reassign by the rules that would have applied
        if (r.lnc_gid, r.pcg_gid) in ex:
            cls[r.lnc_gid] = "antisense"
            continue
        try:
            ls, le = genes.start[r.lnc_gid], genes.end[r.lnc_gid]
            ps, pe = genes.start[r.pcg_gid], genes.end[r.pcg_gid]
        except KeyError:
            cls[r.lnc_gid] = "intergenic"
            continue
        if ls >= ps and le <= pe:
            cls[r.lnc_gid] = "intronic"
        else:
            cls[r.lnc_gid] = "intergenic"
    return cls


def mode_p2():
    with step("90_p0_p2") as st:
        os.makedirs(OUT, exist_ok=True)
        m = load_master()
        genes = genes_meta()
        newcls = reassign(m, genes)
        m["cls_new"] = newcls
        m["cls_new"] = m["cls_new"].fillna(m["cls"])

        counts = m.groupby(["group", "cls", "cls_new"]).size().reset_index(
            name="n")
        counts.to_csv(safe_write_path(os.path.join(
            OUT, "P2_class_reassignment.tsv")), sep="\t", index=False)

        # per-gene features
        cage = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                           low_memory=False)
        cage["base"] = cage.gene_id.astype(str).str.split(".").str[0]
        cage = cage.drop_duplicates("base").set_index("base")
        ch = pd.read_csv(os.path.join(RES, "chromhmm_states.tsv"), sep="\t",
                         low_memory=False).set_index("gene_id")
        act = ch.isin(["TssA", "TssAFlnk"]).mean(axis=1)
        sf = pd.read_csv(os.path.join(RES, "seq_features.tsv"), sep="\t",
                         low_memory=False)
        sf["base"] = sf.gene_id.astype(str).str.split(".").str[0]
        sf = sf.drop_duplicates("base").set_index("base")
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)
        cv = pd.read_csv(os.path.join(RES, "conservation.tsv"), sep="\t",
                         low_memory=False)
        cv["base"] = cv.gene_id.astype(str).str.split(".").str[0]
        cv = cv.drop_duplicates("base").set_index("base")
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)
        ph["base"] = ph.index.astype(str).str.split(".").str[0]
        ph = ph.drop_duplicates("base").set_index("base")
        base = pd.Series(m.index.astype(str).str.split(".").str[0],
                         index=m.index)

        rows = []
        for c, g in m.groupby("cls_new"):
            if len(g) < 20:
                continue
            idx = g.index
            rows.append(dict(
                cls_new=c, n=len(g),
                cage=float(base.loc[idx].map(cage.active).mean()),
                chrom_active=float(act.reindex(idx).mean()),
                gc=float(base.loc[idx].map(sf.gc).mean()),
                cpgoe=float(base.loc[idx].map(sf.cpg_oe).mean()),
                te_cov=float(te.reindex(idx).te_cov.mean()),
                phast=float(base.loc[idx].map(cv.phast_mean).mean()),
                phylop=float(base.loc[idx].map(ph.phylop_mean).mean())))
        S = pd.DataFrame(rows).set_index("cls_new")
        S.to_csv(safe_write_path(os.path.join(
            OUT, "P2_class_stats_by_orientation.tsv")), sep="\t")

        # pair-level table with revised geometry and per-side features
        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        dp["lt"] = dp.lnc_gid.map(m.tss)
        dp["pt"] = dp.pcg_gid.map(m.tss)
        dp["ps"] = dp.pcg_gid.map(m.strand)
        dp["ls"] = dp.lnc_gid.map(m.strand)
        dp = dp.dropna(subset=["lt", "pt", "ps", "ls"])
        dp["h2h"] = np.where(dp["ps"] == "+", dp["lt"] < dp["pt"], dp["lt"] > dp["pt"])
        dp["prom_overlap"] = np.clip(
            np.minimum(dp["lt"] + 200, dp["pt"] + 200)
            - np.maximum(dp["lt"] - 500, dp["pt"] - 500), 0, None)
        dp["lact"] = dp["lnc_gid"].astype(str).str.split(".").str[0].map(
            cage.active)
        dp["pact"] = dp["pcg_gid"].astype(str).str.split(".").str[0].map(
            cage.active)
        dp["l_chrom"] = dp["lnc_gid"].map(act)
        dp["p_chrom"] = dp["pcg_gid"].map(act)
        # CGC partner
        cgc = pd.read_csv(os.path.join(
            WORKDIR, "data", "cosmic", "extracted",
            "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz"), sep="\t",
            usecols=["GENE_SYMBOL"])
        cgc = set(cgc.GENE_SYMBOL.astype(str))
        dp["pcg_name"] = dp["pcg_gid"].map(m.gene_name)
        dp["cgc"] = dp.pcg_name.isin(cgc).astype(int)
        # age
        age = pd.read_csv(os.path.join(
            RES, "comparative", "E_transcript_age.tsv"), sep="\t",
            index_col=0, low_memory=False)
        dp["mouse"] = dp["lnc_gid"].map(age.mouse)
        dp["age_level"] = dp["lnc_gid"].map(age.age_level)
        dp.to_csv(safe_write_path(os.path.join(
            OUT, "P2_pair_table.tsv.gz")), sep="\t", index=False)

        rows = []
        for lab, sub in (("head_to_head", dp[dp.h2h]),
                         ("convergent", dp[~dp.h2h])):
            s2 = sub.dropna(subset=["lact", "pact"])
            a = int(((s2.lact == 1) & (s2.pact == 0)).sum())
            b2 = int(((s2.lact == 0) & (s2.pact == 1)).sum())
            rows.append(dict(
                subset=lab, n=len(sub),
                median_spacing=float(sub.distance.median()),
                median_prom_overlap=float(sub.prom_overlap.median()),
                lnc_cage=float(s2.lact.mean()),
                partner_cage=float(s2.pact.mean()),
                lnc_only=a, partner_only=b2,
                lnc_chrom=float(sub.l_chrom.mean()),
                partner_chrom=float(sub.p_chrom.mean()),
                cgc_frac=float(sub.cgc.mean()),
                mouse_present=float(sub.mouse.mean()),
                older_than_mouse=float((sub.age_level >= 2).mean())))
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "P2_pair_stats_by_orientation.tsv")), sep="\t", index=False)

        for r in S.reset_index().itertuples():
            st.note("{}: n={} CAGE={:.3f} chrom={:.3f} phast={:.3f}".format(
                r.cls_new, r.n, r.cage, r.chrom_active, r.phast))
        for r in T.itertuples():
            st.note("{}: n={} spacing={:.0f} overlap={:.0f} CAGE {:.3f} v "
                    "{:.3f} CGC {:.3f} mouse {:.3f}".format(
                        r.subset, r.n, r.median_spacing,
                        r.median_prom_overlap, r.lnc_cage, r.partner_cage,
                        r.cgc_frac, r.mouse_present))
        st.summary_md(os.path.join(
            RES, "step_summaries", "90_p0_p2.md"),
            "P0-2 orientation-resolved classes",
            df_to_md(S.reset_index(), ".3f") + "\n\n" + df_to_md(T, ".3f"))




def load_rampage():
    import glob, gzip
    out = []
    for path in sorted(glob.glob(os.path.join(
            WORKDIR, "data", "encode_rampage", "bed", "*.bed.gz"))):
        iv = {}
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                c = f[0].replace("chr", "")
                iv.setdefault(c, []).append((int(f[1]), int(f[2])))
        for c in iv:
            iv[c] = np.array(sorted(iv[c]))
        out.append(iv)
    return out


def rag_hit(iv, chrom, pos, half=250):
    t = iv.get(str(chrom).replace("chr", ""))
    if t is None or pd.isna(pos):
        return 0
    s, e = t[:, 0], t[:, 1]
    j = np.searchsorted(e, pos + half, "right") - 1
    if j < 0:
        return 0
    return int(s[j] <= pos + half and e[j] >= pos - half)


def load_fantom_peaks():
    import gzip
    d = {}
    with gzip.open(os.path.join(WORKDIR, "data", "fantom",
                                "hg38_cage_peaks.bed.gz"), "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                score = float(f[4])
            except ValueError:
                continue
            d.setdefault(f[0].replace("chr", ""), []).append(
                (int(f[1]), int(f[2]), score))
    out = {}
    for c, v in d.items():
        v.sort()
        out[c] = (np.array([x[0] for x in v]), np.array([x[1] for x in v]),
                  np.array([x[2] for x in v]))
    return out


def fant_max(peaks, chrom, pos, half=500):
    t = peaks.get(str(chrom).replace("chr", ""))
    if t is None or pd.isna(pos):
        return 0.0
    s, e, sc = t
    lo = np.searchsorted(e, pos - half, "left")
    hi = np.searchsorted(s, pos + half, "right")
    if hi <= lo:
        return 0.0
    return float(sc[lo:hi].max())


def unref_tss(m):
    tx = pd.read_csv(os.path.join(WORKDIR, "data", "genes",
                                  "transcripts.tsv"), sep="\t",
                     low_memory=False)
    tx["tss_u"] = np.where(tx.strand == "+", tx.tx_start, tx.tx_end)
    tx["score"] = (tx.mane.fillna(0).astype(float) * 4
                   + tx.basic.fillna(0).astype(float) * 2
                   + tx.length.fillna(0).astype(float) / 1e6)
    tx = tx.sort_values("score", ascending=False).drop_duplicates("gene_id")
    return tx.set_index("gene_id").tss_u


def mode_p1():
    with step("90_p0_p1") as st:
        os.makedirs(OUT, exist_ok=True)
        m = load_master()
        tu = unref_tss(m)
        m["tss_unref"] = m.index.map(tu)
        st.note("unrefined TSS recovered for {}/{}".format(
            m.tss_unref.notna().sum(), len(m)))
        rag = load_rampage()
        peaks = load_fantom_peaks()

        cage = pd.read_csv(os.path.join(RES, "cage_activity.tsv"), sep="\t",
                           low_memory=False)
        cage["base"] = cage.gene_id.astype(str).str.split(".").str[0]
        cage = cage.drop_duplicates("base").set_index("base")
        base = pd.Series(m.index.astype(str).str.split(".").str[0],
                         index=m.index)
        m["cage_active"] = base.map(cage.active)

        ref_hits = np.zeros(len(m), dtype=int)
        unref_hits = np.zeros(len(m), dtype=int)
        for k, iv in enumerate(rag):
            for i, r in enumerate(m.itertuples()):
                ref_hits[i] += rag_hit(iv, r.chrom, r.tss)
                unref_hits[i] += rag_hit(iv, r.chrom, r.tss_unref)
        m["ramp_ref"] = (ref_hits > 0).astype(int)
        m["ramp_unref"] = (unref_hits > 0).astype(int)

        ref_f = np.array([fant_max(peaks, r.chrom, r.tss)
                          for r in m.itertuples()])
        unref_f = np.array([fant_max(peaks, r.chrom, r.tss_unref)
                            for r in m.itertuples()])
        m["fant_ref"] = (ref_f > 1).astype(int)
        m["fant_unref"] = (unref_f > 1).astype(int)
        m["fant_any2k"] = [fant_max(peaks, r.chrom, r.tss_unref, 2000) > 1
                           for r in m.itertuples()]
        m.to_csv(safe_write_path(os.path.join(
            OUT, "P1_per_gene_scenarios.tsv.gz")), sep="\t")

        rows = []
        for cls in ["divergent_h2h", "antisense", "intronic",
                    "sense_overlapping", "intergenic", "pseudogene",
                    "transcribed_pseudogene", "PCG"]:
            g = m[m.cls == cls] if cls != "divergent_h2h" else \
                m[m.cls == "divergent"]
            if len(g) < 20:
                continue
            rows.append(dict(cls=cls, n=len(g),
                             A_refined_cage=float(g.cage_active.mean()),
                             B_refined_rampage=float(g.ramp_ref.mean()),
                             C_unref_rampage=float(g.ramp_unref.mean()),
                             D_unref_fantom=float(g.fant_unref.mean()),
                             E_unref_fantom_any2k=float(
                                 g[g.fant_any2k].fant_unref.mean())))
        S = pd.DataFrame(rows).set_index("cls")
        S.to_csv(safe_write_path(os.path.join(
            OUT, "P1_class_scenarios.tsv")), sep="\t")

        from scipy.stats import spearmanr
        rho = {}
        for col in ["B_refined_rampage", "C_unref_rampage",
                    "D_unref_fantom", "E_unref_fantom_any2k"]:
            rho[col] = float(spearmanr(S.A_refined_cage, S[col]).statistic)
        # within CAGE-refined subset
        sub = m[m.cage_active == 1]
        rows2 = []
        for cls in ["divergent", "antisense", "intronic",
                    "sense_overlapping", "intergenic"]:
            g = sub[sub.cls == cls]
            if len(g) < 20:
                continue
            rows2.append(dict(cls=cls, n=len(g),
                              refined_cage=float(g.cage_active.mean()),
                              unref_rampage=float(g.ramp_unref.mean()),
                              unref_fantom=float(g.fant_unref.mean())))
        S2 = pd.DataFrame(rows2).set_index("cls")
        S2.to_csv(safe_write_path(os.path.join(
            OUT, "P1_refined_subset.tsv")), sep="\t")
        for r in S.reset_index().itertuples():
            st.note("{}: A={:.3f} B={:.3f} C={:.3f} D={:.3f} E={:.3f}".format(
                r.cls, r.A_refined_cage, r.B_refined_rampage,
                r.C_unref_rampage, r.D_unref_fantom, r.E_unref_fantom_any2k))
        st.note("class-ordering rho vs A: " + ", ".join(
            "{}={:.2f}".format(k, v) for k, v in rho.items()))
        st.summary_md(os.path.join(
            RES, "step_summaries", "90_p0_p1.md"),
            "P0-1 placement versus activity",
            df_to_md(S.reset_index(), ".3f") + "\n\n" + df_to_md(
                S2.reset_index(), ".3f") + "\n\n" + str(rho))




def load_rampage_full():
    import glob, gzip
    out = []
    for path in sorted(glob.glob(os.path.join(
            WORKDIR, "data", "encode_rampage", "bed", "*.bed.gz"))):
        d = {}
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 6:
                    continue
                try:
                    d.setdefault(f[0].replace("chr", ""), []).append(
                        (int(f[1]), int(f[2]), f[5]))
                except ValueError:
                    continue
        iv = {}
        for c, v in d.items():
            v.sort()
            iv[c] = (np.array([x[0] for x in v]),
                     np.array([x[1] for x in v]),
                     np.array([x[2] for x in v]))
        out.append(iv)
    return out


def strict_hit(iv, chrom, pos, strand, partner_pos, half=100):
    t = iv.get(str(chrom).replace("chr", ""))
    if t is None or pd.isna(pos):
        return 0
    st, en, sd = t
    lo = np.searchsorted(st, pos - half - 1000, "left")
    hi = np.searchsorted(st, pos + half, "right")
    for k in range(lo, hi):
        if en[k] < pos - half:
            continue
        if sd[k] != strand:
            continue
        mid = (st[k] + en[k]) / 2.0
        if abs(mid - pos) <= abs(mid - partner_pos):
            return 1
    return 0


def mode_p3():
    with step("90_p0_p3") as st:
        os.makedirs(OUT, exist_ok=True)
        m = load_master()
        rag = load_rampage_full()
        m["strict_hits"] = 0
        for iv in rag:
            for i, r in enumerate(m.itertuples()):
                if pd.isna(r.tss):
                    continue
                pp = m.tss.get(r.partner_pcg) if pd.notna(r.partner_pcg) else r.tss
                if pd.isna(pp):
                    pp = r.tss
                if strict_hit(iv, r.chrom, r.tss, r.strand, pp):
                    m.iloc[i, m.columns.get_loc("strict_hits")] += 1
        m["ramp_strict"] = (m.strict_hits > 0).astype(int)
        m.to_csv(safe_write_path(os.path.join(
            OUT, "P3_rampage_strict.tsv.gz")), sep="\t")

        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        dp["lt"] = dp.lnc_gid.map(m.tss)
        dp["pt"] = dp.pcg_gid.map(m.tss)
        dp["ls"] = dp.lnc_gid.map(m.strand)
        dp["ps"] = dp.pcg_gid.map(m.strand)
        dp = dp.dropna(subset=["lt", "pt", "ls", "ps"])
        dp["h2h"] = np.where(dp["ps"] == "+", dp["lt"] < dp["pt"],
                             dp["lt"] > dp["pt"])
        lh = np.zeros(len(dp), dtype=int)
        ph = np.zeros(len(dp), dtype=int)
        for iv in rag:
            for k, r in enumerate(dp.itertuples()):
                lc = m.chrom.get(r.lnc_gid)
                pc = m.chrom.get(r.pcg_gid)
                lh[k] += strict_hit(iv, lc, r.lt, r.ls, r.pt)
                ph[k] += strict_hit(iv, pc, r.pt, r.ps, r.lt)
        dp["lnc_strict"] = (lh > 0).astype(int)
        dp["pcg_strict"] = (ph > 0).astype(int)

        # GRO-cap: w50 ratios, and spacing-adjusted window
        sig = pd.read_csv(os.path.join(
            RES, "proseq", "divergent_pairs_signal.tsv"), sep="\t")
        sig = sig[sig.label == "K562_GROcap_wTAP"].drop_duplicates(
            ["lnc_gid", "pcg_gid"])
        sig = sig.merge(dp[["lnc_gid", "pcg_gid", "h2h"]],
                        on=["lnc_gid", "pcg_gid"], how="inner")
        rows = []
        for lab, sub in (("head_to_head", sig[sig.h2h]),
                         ("convergent", sig[~sig.h2h])):
            for w in ("w50", "w100"):
                lo = sub[f"lnc_own_{w}"] > 0
                po = sub[f"pcg_own_{w}"] > 0
                both = lo & po
                ratio = (sub.loc[both, f"lnc_own_{w}"]
                         / sub.loc[both, f"pcg_own_{w}"])
                own = ((sub[f"lnc_own_{w}"] > sub[f"lnc_opp_{w}"])
                       & (sub[f"lnc_own_{w}"] + sub[f"lnc_opp_{w}"] > 0))
                rows.append(dict(subset=lab, window=w, n=len(sub),
                                 both=float(both.mean()),
                                 median_ratio=float(ratio.median()),
                                 own_strand=float(own.mean())))
        # spacing-adjusted: use w50 when spacing > 200, else w100
        adj = []
        for r in sig.itertuples():
            w = "w50" if r.distance > 200 else "w100"
            adj.append((getattr(r, f"lnc_own_{w}"),
                        getattr(r, f"pcg_own_{w}"),
                        getattr(r, f"lnc_opp_{w}")))
        adj = pd.DataFrame(adj, columns=["l", "p", "o"])
        adj["h2h"] = sig.h2h.values
        for lab, sub in (("head_to_head", adj[adj.h2h]),
                         ("convergent", adj[~adj.h2h])):
            both = (sub.l > 0) & (sub.p > 0)
            rows.append(dict(subset=lab, window="spacing_adjusted",
                             n=len(sub), both=float(both.mean()),
                             median_ratio=float((sub.l[both]
                                                 / sub.p[both]).median()),
                             own_strand=float(((sub.l > sub.o)
                                               & (sub.l + sub.o > 0)).mean())))
        G = pd.DataFrame(rows)
        G.to_csv(safe_write_path(os.path.join(
            OUT, "P3_grocap_windows.tsv")), sep="\t", index=False)
        R = []
        for lab, sub in (("head_to_head", dp[dp.h2h]),
                         ("convergent", dp[~dp.h2h])):
            R.append(dict(subset=lab, n=len(sub),
                          lnc_strict=float(sub.lnc_strict.mean()),
                          pcg_strict=float(sub.pcg_strict.mean()),
                          lnc_only=int(((sub.lnc_strict == 1)
                                        & (sub.pcg_strict == 0)).sum()),
                          pcg_only=int(((sub.lnc_strict == 0)
                                        & (sub.pcg_strict == 1)).sum())))
        R = pd.DataFrame(R)
        R.to_csv(safe_write_path(os.path.join(
            OUT, "P3_rampage_pairs.tsv")), sep="\t", index=False)
        for r in R.itertuples():
            st.note("RAMPAGE strict {}: lnc {:.1%} partner {:.1%} only {} v "
                    "{}".format(r.subset, r.lnc_strict, r.pcg_strict,
                                r.lnc_only, r.pcg_only))
        for r in G.itertuples():
            st.note("GRO-cap {} {}: both {:.1%} ratio {:.2f} own {:.1%}"
                    .format(r.subset, r.window, r.both, r.median_ratio,
                            r.own_strand))
        st.summary_md(os.path.join(
            RES, "step_summaries", "90_p0_p3.md"),
            "P0-3 window controls",
            df_to_md(R, ".3f") + "\n\n" + df_to_md(G, ".3f"))




def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def mode_p4():
    with step("90_p0_p4") as st:
        os.makedirs(OUT, exist_ok=True)
        cmp34 = _load_module(os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "34_comparative.py"), "cmp34")
        CMP = os.path.join(WORKDIR, "data", "comparative")
        db = __import__("pylib.chainlift", fromlist=["ChainDB"]).ChainDB(
            cmp34.open_chain(os.path.join(
                CMP, "chains", "hg38ToMm39.over.chain.gz")))
        mouse = cmp34.load_tss(os.path.join(
            CMP, "gtf", "Mus_musculus.GRCm39.116.gtf.gz"))
        for c in mouse:
            for s in ("+", "-"):
                mouse[c][s] = np.array(sorted(mouse[c][s]))

        m = load_master()
        m = m[m.group == "lnc"].dropna(subset=["tss"]).copy()
        age = pd.read_csv(os.path.join(
            RES, "comparative", "E_transcript_age.tsv"), sep="\t",
            index_col=0, low_memory=False)
        m["mouse_obs"] = m.index.map(age.mouse)

        # gene-density decile from PCG TSS within +/-500 kb
        pcg = load_master()
        pcg = pcg[(pcg.group == "pcg") & pcg.tss.notna()]
        dens = pd.Series(0, index=m.index, dtype=float)
        for c, g in m.groupby("chrom"):
            pv = np.sort(pcg.tss[pcg.chrom == c].values)
            if not len(pv):
                continue
            pos = g.tss.values
            lo = np.searchsorted(pv, pos - 500000, "left")
            hi = np.searchsorted(pv, pos + 500000, "right")
            dens.loc[g.index] = hi - lo
        m["dens_bin"] = pd.qcut(dens, 10, labels=False, duplicates="drop").astype(int)

        def match(chrom, pos, strand, win=1000):
            lf = db.lift(str(chrom), int(pos) - 1, int(pos))
            if lf is None:
                return 0
            qn = lf["chrom"]
            qp = lf["start"]
            qs = lf["strand"]
            arr = mouse.get(qn.replace("chr", ""), {}).get(qs)
            if arr is None or not len(arr):
                return 0
            j = np.searchsorted(arr, qp)
            for k in (j - 1, j, j + 1):
                if 0 <= k < len(arr) and abs(int(arr[k]) - qp) <= win:
                    return 1
            return 0

        obs = np.array([match(r.chrom, r.tss, r.strand)
                        for r in m.itertuples()])
        m["match_obs"] = obs
        st.note("observed mouse match {:.1%}; table mouse presence {:.1%}"
                .format(obs.mean(), m.mouse_obs.mean()))

        rng = np.random.default_rng(42)
        nperm = 50
        null_rates = {c: [] for c in m.cls.unique()}
        keys = list(m.groupby(["chrom", "dens_bin"]).groups.items())
        base_strand = m.strand.values.copy()
        base_tss = m.tss.values.copy()
        cls = m.cls.values
        chroms = m.chrom.values
        bins = m.dens_bin.values
        for _ in range(nperm):
            perm_tss = base_tss.copy()
            for (c, b), idx in keys:
                pos_idx = m.index.get_indexer(pd.Index(idx))
                perm_tss[pos_idx] = rng.permutation(base_tss[pos_idx])
            rate = {}
            for i in range(len(m)):
                rate.setdefault(cls[i], []).append(
                    match(chroms[i], perm_tss[i], base_strand[i]))
            for c in rate:
                null_rates[c].append(float(np.mean(rate[c])))

        rows = []
        for c in sorted(null_rates):
            nl = np.array(null_rates[c])
            g = m[m.cls == c]
            rows.append(dict(cls=c, n=len(g),
                             observed=float(g.match_obs.mean()),
                             null_mean=float(nl.mean()),
                             null_lo=float(np.percentile(nl, 2.5)),
                             null_hi=float(np.percentile(nl, 97.5)),
                             excess=float(g.match_obs.mean() - nl.mean()),
                             z=float((g.match_obs.mean() - nl.mean())
                                     / (nl.std() + 1e-9))))
        T = pd.DataFrame(rows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "P4_age_null.tsv")), sep="\t", index=False)
        m[["cls", "chrom", "tss", "strand", "dens_bin", "mouse_obs",
           "match_obs"]].to_csv(safe_write_path(os.path.join(
               OUT, "P4_age_per_gene.tsv.gz")), sep="\t")
        for r in T.itertuples():
            st.note("{}: n={} obs={:.1%} null={:.1%} [{:.1%},{:.1%}] "
                    "excess={:+.1%} z={:.1f}".format(
                        r.cls, r.n, r.observed, r.null_mean, r.null_lo,
                        r.null_hi, r.excess, r.z))
        st.summary_md(os.path.join(
            RES, "step_summaries", "90_p0_p4.md"),
            "P0-4 age identifiability",
            df_to_md(T, ".4f"))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "p2"
    {"p2": mode_p2, "p1": mode_p1, "p3": mode_p3, "p4": mode_p4}[mode]()
