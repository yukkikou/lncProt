#!/usr/bin/env python
"""37_mouse_taxonomy.py — mouse five-class lncRNA taxonomy (mm10 / GRCm38).

Reconstructs the coordinate-defined lncRNA birth classes used for human
(scripts/04_classify_loci.py) in mouse, using FANTOM5 mm10 CAGE peaks and
Ensembl GRCm38 release 102 (the last Ensembl annotation on mm10, chosen so the
annotation and CAGE peaks share one assembly; verified with Actb/Gapdh).

Classes (same priority order as human 04):
  divergent         lnc major TSS opposite-strand <=1 kb from a PCG TSS
  antisense         lnc exons overlap PCG exons on the opposite strand
  intronic          lnc gene body fully inside a PCG body, no exon overlap
  sense_overlapping lnc exons overlap PCG exons on the same strand
  intergenic        none of the above

CLI modes
  prepare    parse GTF -> gene table, transcripts, exons, PCG exon/TSS sets
  classify   CAGE-refined major TSS + five classes
  activity   map CAGE peaks/TPM to genes; activity/breadth/tau by class
  compare    class composition + activity vs the human table
  all        prepare + classify (activity is run separately when TPM is ready)

Examples
  PY=/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/bin/python
  $PY scripts/37_mouse_taxonomy.py prepare
  $PY scripts/37_mouse_taxonomy.py classify
  $PY scripts/37_mouse_taxonomy.py activity      # needs the complete TPM matrix
  $PY scripts/37_mouse_taxonomy.py compare
"""
import argparse
import gzip
import os
import re
import subprocess
import sys
import zlib

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR, df_to_md

DATA = os.path.join(WORKDIR, "data")
RES = os.path.join(WORKDIR, "results", "mouse_taxonomy")
TMP = os.path.join(DATA, "mouse_taxonomy_tmp")
GTF = os.path.join(DATA, "comparative", "gtf",
                   "Mus_musculus.GRCm38.102.gtf.gz")
CAGE_BED = os.path.join(DATA, "fantom_mouse", "mm10_CAGE_peaks.bed.gz")
CAGE_TPM = os.path.join(DATA, "fantom_mouse", "mm10_CAGE_tpm.osc.txt.gz")
HUMAN_MASTER = os.path.join(WORKDIR, "results", "gene_master_table.tsv")
HUMAN_ACT = os.path.join(WORKDIR, "results", "cage_activity_by_class.tsv")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"

DIV_WIN = 1000
MATCH_WIN = 250
INTERG_MIN = 10000
TPM_MIN = 1.0
EXPECTED_TPM_BYTES = 450_000_000   # official file ~459.6 MB (tpm.osc)
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
LNC_TYPES = {
    "lincRNA", "antisense", "processed_transcript", "sense_intronic",
    "sense_overlapping", "bidirectional_promoter_lncRNA", "macro_lncRNA",
    "3prime_overlapping_ncRNA", "non_coding",
}
QC_GENES = ["Actb", "Gapdh", "Neat1", "Malat1", "Xist", "Tsix", "H19",
            "Gas5", "Snhg1", "Meg3", "Kcnq1ot1", "Airn", "Pvt1"]

_ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
_TAG_RE = re.compile(r'tag "([^"]*)"')


# ------------------------------------------------------------- utilities ---
def group_of(biotype):
    if biotype == "protein_coding":
        return "pcg"
    if biotype in LNC_TYPES:
        return "lnc"
    if "pseudogene" in biotype:
        return "pseudo"
    return "other"


def norm_chrom(c):
    if c == "MT":
        return "chrM"
    return c if c.startswith("chr") else "chr" + c


def parse_attrs(s):
    return dict(_ATTR_RE.findall(s))


def bt(*args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", WORKDIR + ":" + WORKDIR, SIF_BT, "bedtools", *args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("bedtools failed: " + " ".join(args) +
                           "\n" + r.stderr[:2000])
    return r.stdout


def write_bed(df, path, cols):
    path = safe_write_path(path)
    d = df.sort_values([cols[0], cols[1]], kind="mergesort")
    d[cols].to_csv(path, sep="\t", index=False, header=False)


def tmp(name):
    os.makedirs(TMP, exist_ok=True)
    return os.path.join(TMP, name)


# --------------------------------------------------------------- prepare ---
def cmd_prepare():
    with step("37_mouse_prepare") as st:
        if not os.path.exists(GTF):
            raise FileNotFoundError(GTF)
        genes, txs, exons = [], [], []
        with gzip.open(GTF, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                feat = f[2]
                if feat not in ("gene", "transcript", "exon"):
                    continue
                a = parse_attrs(f[8])
                gid = a.get("gene_id")
                if not gid:
                    continue
                chrom = norm_chrom(f[0])
                start, end, strand = int(f[3]), int(f[4]), f[6]
                if feat == "gene":
                    genes.append(dict(
                        gene_id=gid, gene_name=a.get("gene_name", gid),
                        gene_type=a.get("gene_biotype", ""), chrom=chrom,
                        start=start, end=end, strand=strand))
                elif feat == "transcript":
                    tags = _TAG_RE.findall(f[8])
                    tss = start if strand == "+" else end
                    txs.append(dict(
                        gene_id=gid, transcript_id=a.get("transcript_id", ""),
                        chrom=chrom, start=start, end=end, strand=strand,
                        tss=tss, length=end - start + 1,
                        tsl=a.get("transcript_support_level", "NA"),
                        basic=int(bool({"basic", "gencode_basic"} & set(tags))),
                        canonical=int(bool({"Ensembl_canonical", "canonical"}
                                           & set(tags)))))
                else:
                    exons.append((chrom, start - 1, end,
                                  gid + "|" + a.get("transcript_id", ""),
                                  0, strand))
        st.note("GTF parsed: genes=" + str(len(genes)) +
                " transcripts=" + str(len(txs)) +
                " exons=" + str(len(exons)))

        g = pd.DataFrame(genes).set_index("gene_id")
        g.index.name = "gene_id"
        g["tss"] = np.where(g["strand"] == "+", g["start"], g["end"])
        g["length"] = g["end"] - g["start"] + 1
        g["group"] = g["gene_type"].map(group_of)
        ex = pd.DataFrame(exons, columns=["chrom", "start", "end", "gid_tx",
                                          "score", "strand"])
        ex["gene_id"] = ex["gid_tx"].str.split("|").str[0]
        g["n_exons"] = g.index.map(
            ex.groupby("gene_id").size()).fillna(0).astype(int)
        tx = pd.DataFrame(txs)

        os.makedirs(RES, exist_ok=True)
        g.to_csv(safe_write_path(os.path.join(RES, "gene_table.tsv")),
                 sep="\t")
        tx.to_csv(safe_write_path(os.path.join(RES, "transcripts.tsv")),
                  sep="\t", index=False)
        ex[["chrom", "start", "end", "gid_tx", "score", "strand"]].to_csv(
            safe_write_path(os.path.join(RES, "exons.bed.gz")),
            sep="\t", index=False, header=False, compression="gzip")

        counts = g["group"].value_counts()
        st.note("groups: " + ", ".join(f"{k}={v}" for k, v in counts.items()))
        st.note("lnc biotypes: " + ", ".join(
            f"{k}={v}" for k, v in
            g[g["group"] == "lnc"]["gene_type"].value_counts().items()))

        # PCG exon / TSS sets requested by the analysis design
        pcg = g[g["group"] == "pcg"]
        pex = ex[ex["gene_id"].isin(pcg.index)].copy()
        write_bed(pex, os.path.join(RES, "pcg_exons.bed"),
                  ["chrom", "start", "end", "gene_id", "score", "strand"])
        pt = pcg.reset_index()
        pt["s0"] = pt["tss"] - 1
        pt["e0"] = pt["tss"]
        pt["score"] = 0
        write_bed(pt, os.path.join(RES, "pcg_tss.bed"),
                  ["chrom", "s0", "e0", "gene_id", "score", "strand"])
        st.note("PCG exon/ TSS sets written: " + str(len(pcg)) + " genes, " +
                str(len(pex)) + " exons")

        # ---- assembly sanity check against known mm10 coordinates
        for name, chrom, tss in (("Actb", "chr5", 142903115),
                                 ("Gapdh", "chr6", 125161715)):
            r = g[g["gene_name"] == name]
            if not len(r):
                st.note("QC " + name + ": NOT FOUND")
                continue
            r = r.iloc[0]
            ok = (r["chrom"] == chrom and abs(int(r["start"]) - tss) < 2000)
            st.note(("QC " + name + ": " + str(r["chrom"]) + ":" +
                     str(int(r["start"])) + "-" + str(int(r["end"])) +
                     " strand=" + r["strand"] + " -> ") +
                    ("OK (mm10)" if ok else "MISMATCH"))
            if not ok:
                raise AssertionError("assembly mismatch for " + name)

        body = ("| group | n |\n|---|---|\n" +
                "\n".join(f"| {k} | {v} |" for k, v in counts.items()) +
                "\n\nlncRNA biotypes:\n\n" +
                df_to_md(g[g["group"] == "lnc"]["gene_type"]
                         .value_counts().rename_axis("biotype")
                         .reset_index(name="n"), max_rows=20) + "\n")
        st.summary_md("results/step_summaries/37_mouse_prepare.md",
                      "37 mouse gene table (GRCm38.102)", body)


# --------------------------------------------------------------- classify ---
def load_gene_table():
    return pd.read_csv(os.path.join(RES, "gene_table.tsv"), sep="\t",
                       index_col=0, low_memory=False)


def load_cage():
    rows = []
    with gzip.open(CAGE_BED, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[0].startswith(("#", "track", "browser")):
                continue
            chrom, s, e = f[0], int(f[1]), int(f[2])
            pid, score, strand = f[3], float(f[4]), f[5]
            if len(f) >= 8 and f[6].isdigit() and f[7].isdigit():
                summit = (int(f[6]) + int(f[7])) // 2
            else:
                summit = (s + e) // 2
            rows.append((chrom, summit, pid, score, strand))
    return pd.DataFrame(rows, columns=["chrom", "summit", "peak_id",
                                       "score", "strand"])


def major_tss(genes, tx, cage, st):
    """Per gene: CAGE-refined major TSS if a same-strand peak summit lies
    within +/-MATCH_WIN of any transcript TSS; else best-transcript TSS."""
    tx = tx.copy()
    tx["tsl_n"] = pd.to_numeric(tx["tsl"], errors="coerce").fillna(9)
    tx["rank"] = ((tx["canonical"] * 4 + tx["basic"] * 2) * 1e9
                  - tx["tsl_n"] * 1e6 + tx["length"])
    best = tx.sort_values("rank", ascending=False).groupby("gene_id").first()
    gid2best_tss = best["tss"].to_dict()

    cage_by_gene = {}
    if cage is not None and len(cage):
        cand = tx[["chrom", "tss", "gene_id", "strand"]].copy()
        cand["sc"] = 0
        cand["w_start"] = (cand["tss"] - 1 - MATCH_WIN).clip(lower=0)
        cand["w_end"] = cand["tss"] - 1 + MATCH_WIN
        write_bed(cand, tmp("cand_win.bed"),
                  ["chrom", "w_start", "w_end", "gene_id", "sc", "strand"])
        c = cage.copy()
        c["s"] = c["summit"] - 1
        c["e"] = c["summit"]
        write_bed(c, tmp("cage_tss.bed"),
                  ["chrom", "s", "e", "peak_id", "score", "strand"])
        pairs = bt("intersect", "-a", tmp("cand_win.bed"),
                   "-b", tmp("cage_tss.bed"), "-wa", "-wb", "-s")
        hits = []
        for ln in pairs.splitlines():
            f = ln.split("\t")
            if len(f) < 12:
                continue
            hits.append((f[3], f[9], float(f[10]), int(f[7]) + 1))
        hd = pd.DataFrame(hits, columns=["gid", "peak_id", "score", "summit"])
        if len(hd):
            idx = hd.groupby("gid")["score"].idxmax()
            cage_by_gene = hd.loc[idx].set_index("gid")[
                ["peak_id", "score", "summit"]].to_dict("index")

    out = {}
    for gid, row in genes.iterrows():
        tss0 = gid2best_tss.get(
            gid, row["start"] if row["strand"] == "+" else row["end"])
        c = cage_by_gene.get(gid)
        if c:
            out[gid] = (c["summit"], 1, c["peak_id"], c["score"])
        else:
            out[gid] = (tss0, 0, "", np.nan)
    return out


def cmd_classify():
    with step("37_mouse_classify") as st:
        genes = load_gene_table()
        tx = pd.read_csv(os.path.join(RES, "transcripts.tsv"), sep="\t")
        ex = pd.read_csv(os.path.join(RES, "exons.bed.gz"), sep="\t",
                         header=None, compression="gzip",
                         names=["chrom", "start", "end", "gid_tx", "score",
                                "strand"])
        ex["gid"] = ex["gid_tx"].str.split("|").str[0]
        st.note("loaded genes=" + str(len(genes)) +
                " transcripts=" + str(len(tx)) + " exons=" + str(len(ex)))

        cage = load_cage()
        st.note("CAGE peaks=" + str(len(cage)))
        tssmap = major_tss(genes, tx, cage, st)
        n_ref = sum(1 for v in tssmap.values() if v[1] == 1)
        st.note("CAGE-refined major TSS: " + str(n_ref) + "/" +
                str(len(tssmap)) + " (" + f"{n_ref / len(tssmap):.1%}" + ")")

        lnc = genes[genes["group"] == "lnc"]
        pcg = genes[genes["group"] == "pcg"]
        pseudo = genes[genes["group"] == "pseudo"]
        lnc_ex = ex[ex["gid"].isin(lnc.index)].copy()
        pcg_ex = ex[ex["gid"].isin(pcg.index)].copy()

        # ---- exon-level lnc x PCG overlap pairs
        write_bed(lnc_ex, tmp("lnc_ex.bed"),
                  ["chrom", "start", "end", "gid", "score", "strand"])
        write_bed(pcg_ex, tmp("pcg_ex.bed"),
                  ["chrom", "start", "end", "gid", "score", "strand"])
        ov = bt("intersect", "-a", tmp("lnc_ex.bed"), "-b", tmp("pcg_ex.bed"),
                "-wa", "-wb")
        pairs, strands = set(), {}
        for ln in ov.splitlines():
            f = ln.split("\t")
            if len(f) < 12:
                continue
            a, b = f[3], f[9]
            pairs.add((a, b))
            strands[(a, b)] = (f[5], f[11])
        pd.DataFrame(sorted(pairs), columns=["lnc_gid", "pcg_gid"]).to_csv(
            safe_write_path(os.path.join(RES, "exon_overlap_pairs.tsv")),
            sep="\t", index=False)
        st.note("exon_overlap_pairs=" + str(len(pairs)))
        pairs_by_lnc = {}
        for a, b in pairs:
            sa, sb = strands[(a, b)]
            pairs_by_lnc.setdefault(a, []).append((b, sa == sb))

        # ---- gene-body containment
        lb = lnc.reset_index().rename(columns={"chrom": "c", "start": "s",
                                               "end": "e"})
        lb["s0"] = lb["s"] - 1
        lb["sc"] = 0
        write_bed(lb, tmp("lnc_g.bed"),
                  ["c", "s0", "e", "gene_id", "sc", "strand"])
        pb = pcg.reset_index().rename(columns={"chrom": "c", "start": "s",
                                               "end": "e"})
        pb["s0"] = pb["s"] - 1
        pb["sc"] = 0
        write_bed(pb, tmp("pcg_g.bed"),
                  ["c", "s0", "e", "gene_id", "sc", "strand"])
        cont = bt("intersect", "-a", tmp("lnc_g.bed"), "-b", tmp("pcg_g.bed"),
                  "-wa", "-wb")
        contain = {}
        for ln in cont.splitlines():
            f = ln.split("\t")
            if len(f) < 12:
                continue
            if int(f[1]) >= int(f[7]) and int(f[2]) <= int(f[8]):
                contain.setdefault(f[3], f[9])
        st.note("contained_in_pcg=" + str(len(contain)))

        # ---- TSS BEDs
        def tss_bed(df, path):
            d = df.reset_index()
            d["tss"] = [tssmap[g][0] for g in d["gene_id"]]
            d["s0"] = d["tss"] - 1
            d["e0"] = d["tss"]
            d["sc"] = 0
            write_bed(d, path,
                      ["chrom", "s0", "e0", "gene_id", "sc", "strand"])

        tss_bed(lnc, tmp("lnc_tss.bed"))
        tss_bed(pcg, tmp("pcg_tss.bed"))
        tss_bed(pseudo, tmp("pseudo_tss.bed"))

        cl = bt("closest", "-a", tmp("lnc_tss.bed"), "-b", tmp("pcg_tss.bed"),
                "-S", "-d", "-t", "first")
        div, div_str = {}, {}
        for ln in cl.splitlines():
            f = ln.split("\t")
            if len(f) < 13 or f[6] == ".":
                continue
            div[f[3]] = int(f[12])
            div_str[f[3]] = f[9]
        cs = bt("closest", "-a", tmp("lnc_tss.bed"), "-b", tmp("pcg_tss.bed"),
                "-s", "-d", "-t", "first")
        same_d = {}
        for ln in cs.splitlines():
            f = ln.split("\t")
            if len(f) >= 13 and f[6] != ".":
                same_d[f[3]] = int(f[12])
        gb = bt("closest", "-a", tmp("lnc_g.bed"), "-b", tmp("pcg_g.bed"),
                "-d", "-t", "first")
        body_d = {}
        for ln in gb.splitlines():
            f = ln.split("\t")
            if len(f) >= 13 and f[6] != ".":
                body_d[f[3]] = int(f[12])

        # ---- pseudogene association
        pdb = pseudo.reset_index().rename(columns={"chrom": "c", "start": "s",
                                                  "end": "e"})
        pdb["s0"] = pdb["s"] - 1
        pdb["sc"] = 0
        write_bed(pdb, tmp("pseudo_g.bed"),
                  ["c", "s0", "e", "gene_id", "sc", "strand"])
        ov2 = bt("intersect", "-a", tmp("lnc_g.bed"), "-b", tmp("pseudo_g.bed"),
                 "-wa", "-wb")
        pseudo_ov = {}
        for ln in ov2.splitlines():
            f = ln.split("\t")
            if len(f) >= 12:
                pseudo_ov.setdefault(f[3], f[9])
        pc = bt("closest", "-a", tmp("lnc_tss.bed"),
                "-b", tmp("pseudo_tss.bed"), "-d", "-t", "first")
        pseudo_tss_d = {}
        for ln in pc.splitlines():
            f = ln.split("\t")
            if len(f) >= 13 and f[6] != ".":
                pseudo_tss_d[f[3]] = int(f[12])

        # ---- assemble classes
        recs, div_rows = [], []
        for gid, row in lnc.iterrows():
            tss, refined, pid, pscore = tssmap[gid]
            if gid in div and div[gid] <= DIV_WIN:
                cls, partner = "divergent", div_str[gid]
                nm = (genes.at[partner, "gene_name"]
                      if partner in genes.index else "?")
                div_rows.append((gid, partner, div[gid], row["strand"], nm))
            else:
                plist = pairs_by_lnc.get(gid, [])
                anti = [b for b, same in plist if not same]
                sense = [b for b, same in plist if same]
                if anti:
                    cls, partner = "antisense", anti[0]
                elif gid in contain:
                    cls, partner = "intronic", contain[gid]
                elif sense:
                    cls, partner = "sense_overlapping", sense[0]
                else:
                    cls, partner = "intergenic", ""
            recs.append(dict(
                gene_id=gid, gene_name=row["gene_name"],
                gene_type=row["gene_type"], group="lnc", cls=cls,
                chrom=row["chrom"], tss=tss, strand=row["strand"],
                cage_refined=refined, cage_peak_id=pid,
                cage_peak_score=pscore, partner_pcg=partner,
                div_dist_opp=div.get(gid), same_tss_dist=same_d.get(gid),
                pcg_body_dist=body_d.get(gid),
                pseudo_assoc=1 if (gid in pseudo_ov or
                                   pseudo_tss_d.get(gid, 1e9) <= 2000) else 0,
                pseudo_overlap_gid=pseudo_ov.get(gid, ""),
                pseudo_tss_dist=pseudo_tss_d.get(gid), transcribed_pseudo=0))

        for gid, row in pcg.iterrows():
            tss, refined, pid, pscore = tssmap[gid]
            recs.append(dict(
                gene_id=gid, gene_name=row["gene_name"],
                gene_type=row["gene_type"], group="pcg", cls="PCG",
                chrom=row["chrom"], tss=tss, strand=row["strand"],
                cage_refined=refined, cage_peak_id=pid,
                cage_peak_score=pscore, partner_pcg="", transcribed_pseudo=0))
        for gid, row in pseudo.iterrows():
            tss, refined, pid, pscore = tssmap[gid]
            tp = int(str(row["gene_type"]).startswith("transcribed_"))
            recs.append(dict(
                gene_id=gid, gene_name=row["gene_name"],
                gene_type=row["gene_type"], group="pseudo",
                cls="transcribed_pseudogene" if tp else "pseudogene",
                chrom=row["chrom"], tss=tss, strand=row["strand"],
                cage_refined=refined, cage_peak_id=pid,
                cage_peak_score=pscore, partner_pcg="",
                transcribed_pseudo=tp))

        m = pd.DataFrame(recs).set_index("gene_id")
        m["div_500"] = [
            (1 if (d is not None and not pd.isna(d) and d <= 500) else 0)
            if g == "lnc" else np.nan
            for d, g in zip(m["div_dist_opp"], m["group"])]
        m["div_2000"] = [
            (1 if (d is not None and not pd.isna(d) and d <= 2000) else 0)
            if g == "lnc" else np.nan
            for d, g in zip(m["div_dist_opp"], m["group"])]
        m["intergenic_near"] = [
            1 if (g == "intergenic" and (b is None or pd.isna(b) or
                                         b <= INTERG_MIN))
            else (0 if g == "intergenic" else np.nan)
            for g, b in zip(m["cls"], m["pcg_body_dist"])]
        pd.DataFrame(div_rows,
                     columns=["lnc_gid", "pcg_gid", "distance", "lnc_strand",
                              "pcg_name"]).to_csv(
            safe_write_path(os.path.join(RES, "divergent_pairs.tsv")),
            sep="\t", index=False)

        dist = m[m["group"] == "lnc"]["cls"].value_counts()
        st.note("lnc classes: " + ", ".join(f"{k}={v}" for k, v in
                                            dist.items()))
        if int(dist.sum()) != len(lnc):
            raise AssertionError("class sum " + str(dist.sum()) +
                                 " != lnc " + str(len(lnc)))
        for c in CLASSES:
            if c not in dist.index:
                st.note("WARNING: class " + c + " has zero members")
        pcg_ref = m[(m["group"] == "pcg") & (m["cage_refined"] == 1)].shape[0]
        lnc_ref = m[(m["group"] == "lnc") & (m["cage_refined"] == 1)].shape[0]
        st.note("CAGE-refined: PCG " + f"{pcg_ref / len(pcg):.1%}" +
                ", lnc " + f"{lnc_ref / len(lnc):.1%}")

        m.to_csv(safe_write_path(os.path.join(RES, "gene_master_table.tsv")),
                 sep="\t")
        os.makedirs(os.path.join(RES, "classes"), exist_ok=True)
        for cls, sub in m[m["group"] == "lnc"].groupby("cls"):
            d = sub.reset_index()
            d["s0"], d["e0"] = d["tss"] - 1, d["tss"]
            write_bed(d, os.path.join(RES, "classes", cls + ".bed"),
                      ["chrom", "s0", "e0", "gene_id", "strand"])

        # ---- named-locus QC
        for nm in QC_GENES:
            r = m[m["gene_name"] == nm]
            if not len(r):
                st.note("QC " + nm + ": not in pcg/lnc/pseudo table")
                continue
            r = r.iloc[0]
            st.note("QC " + nm + ": group=" + str(r["group"]) + " class=" +
                    str(r["cls"]) + " " + str(r["chrom"]) + ":" +
                    str(int(r["tss"])) + r["strand"] + " refined=" +
                    str(int(r["cage_refined"])) +
                    (" partner=" + str(r["partner_pcg"])
                     if r["partner_pcg"] else ""))

        body = ("| class | n | fraction |\n|---|---|---|\n" +
                "\n".join(f"| {k} | {v} | {v / len(lnc):.3f} |"
                          for k, v in dist.items()) +
                f"\n\nTotal lncRNA: {len(lnc)}; PCG: {len(pcg)}; "
                f"pseudogene: {len(pseudo)}.\n\n"
                f"CAGE-refined major TSS: PCG {pcg_ref / len(pcg):.1%}, "
                f"lncRNA {lnc_ref / len(lnc):.1%}.\n")
        st.summary_md("results/step_summaries/37_mouse_classify.md",
                      "37 mouse five-class taxonomy", body)


# --------------------------------------------------------------- activity ---
def osc_peak_stats(path, tpm_min=TPM_MIN, chunk=4096):
    """Streaming per-peak stats from a FANTOM .osc TPM matrix.

    Returns DataFrame indexed by peak id with max_tpm, mean_tpm, breadth
    (# samples > tpm_min) and tau. Memory stays bounded (chunked)."""
    ids, blocks = [], {"max_tpm": [], "mean_tpm": [], "breadth": [],
                       "tau": []}
    header, sel = None, None
    buf_ids, buf_vals = [], []
    n_samples = 0

    def flush():
        if not buf_ids:
            return
        arr = np.asarray(buf_vals, dtype=np.float32)
        mx = arr.max(axis=1)
        ids.extend(buf_ids)
        blocks["max_tpm"].append(mx)
        blocks["mean_tpm"].append(arr.mean(axis=1))
        blocks["breadth"].append((arr > tpm_min).sum(axis=1).astype(np.int32))
        denom = np.maximum(mx, 1e-9)
        x = arr / denom[:, None]
        tau = (1.0 - x).sum(axis=1) / max(arr.shape[1] - 1, 1)
        tau[mx <= 0] = np.nan
        blocks["tau"].append(tau.astype(np.float32))
        buf_ids.clear()
        buf_vals.clear()

    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                sel = [i for i, nm in enumerate(header)
                       if i > 0 and nm.startswith("tpm.")]
                if not sel:
                    sel = list(range(1, len(header)))
                n_samples = len(sel)
                continue
            pid = f[0]
            if "::" not in pid:          # 01STAT/02STAT metadata rows
                continue
            vals = []
            for i in sel:
                v = f[i] if i < len(f) else "0"
                vals.append(v if v not in ("", "NA") else "0")
            buf_ids.append(pid)
            buf_vals.append(vals)
            if len(buf_ids) >= chunk:
                flush()
    flush()
    out = pd.DataFrame({"peak_id": ids})
    for k in ("max_tpm", "mean_tpm", "breadth", "tau"):
        out[k] = np.concatenate(blocks[k]) if blocks[k] else np.array([])
    out = out.set_index("peak_id")
    out.attrs["n_samples"] = n_samples
    return out


def cmd_activity():
    with step("37_mouse_activity") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        mode = "peaks-only"
        stats = None
        if os.path.exists(CAGE_TPM):
            sz = os.path.getsize(CAGE_TPM)
            st.note("TPM file size=" + f"{sz / 1e6:.1f} MB")
            if sz >= EXPECTED_TPM_BYTES:
                try:
                    stats = osc_peak_stats(CAGE_TPM)
                    st.note("TPM matrix parsed: " + str(len(stats)) +
                            " peaks x " + str(stats.attrs.get("n_samples")) +
                            " samples")
                    mode = "TPM"
                except (EOFError, OSError, ValueError, zlib.error) as e:
                    st.note("TPM parse failed (" + type(e).__name__ +
                            ") — falling back to peaks-only")
            else:
                st.note("TPM still incomplete (" +
                        f"{sz / EXPECTED_TPM_BYTES:.0%}" +
                        " of expected) — falling back to peaks-only")

        act = pd.DataFrame(index=meta.index)
        if mode == "TPM":
            st_ = stats.reindex(meta["cage_peak_id"].values)
            act["active"] = (st_["max_tpm"] > TPM_MIN).astype(int).values
            act["max_tpm"] = st_["max_tpm"].values
            act["mean_tpm"] = st_["mean_tpm"].values
            act["breadth"] = st_["breadth"].values
            act["tau"] = st_["tau"].values
            linked = meta["cage_peak_id"].isin(stats.index)
            st.note("promoters linked to TPM rows: " + str(int(linked.sum())) +
                    " (" + f"{linked.mean():.0%}" + ")")
            stats.to_csv(safe_write_path(
                os.path.join(RES, "cage_peak_stats.tsv.gz")), sep="\t",
                compression="gzip")
        else:
            act["active"] = (meta["cage_refined"] == 1).astype(int).values
            act["max_tpm"] = meta["cage_peak_score"].values
            act["mean_tpm"] = np.nan
            act["breadth"] = np.nan
            act["tau"] = np.nan

        act["activity_mode"] = mode
        out = meta.join(act)
        out.to_csv(safe_write_path(os.path.join(RES, "cage_activity.tsv")),
                   sep="\t")

        rows = []
        for g in ["PCG"] + CLASSES + ["pseudogene", "transcribed_pseudogene"]:
            sub = out[out["group"] == "pcg"] if g == "PCG" else \
                out[out["cls"] == g]
            if not len(sub):
                continue
            rows.append(dict(
                cls=g, n=len(sub), active_frac=sub["active"].mean(),
                log10_maxtpm=np.log10(sub["max_tpm"].fillna(0) + 1).mean(),
                tau=sub["tau"].mean() if mode == "TPM" else np.nan,
                breadth=sub["breadth"].mean() if mode == "TPM" else np.nan))
        sm = pd.DataFrame(rows)
        sm["mode"] = mode
        sm.to_csv(safe_write_path(
            os.path.join(RES, "cage_activity_by_class.tsv")), sep="\t",
            index=False)
        st.note("mode=" + mode + "; activity by class: " + "; ".join(
            f"{r.cls}: active={r.active_frac:.2f}" for r in sm.itertuples()))
        st.summary_md("results/step_summaries/37_mouse_activity.md",
                      "37 mouse CAGE activity (" + mode + ")",
                      df_to_md(sm, ".2f"))


# ---------------------------------------------------------------- compare ---
def cmd_compare():
    with step("37_mouse_compare") as st:
        mm = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                         index_col=0, low_memory=False)
        hu = pd.read_csv(HUMAN_MASTER, sep="\t", index_col=0,
                         low_memory=False)
        rows = []
        for sp, df in (("human", hu), ("mouse", mm)):
            lnc = df[df["group"] == "lnc"]
            for c in CLASSES:
                n = int((lnc["cls"] == c).sum())
                rows.append(dict(species=sp, cls=c, n=n,
                                 fraction=n / len(lnc)))
            rows.append(dict(species=sp, cls="TOTAL_LNC", n=len(lnc),
                             fraction=1.0))
            for g in ("pcg", "pseudo"):
                rows.append(dict(species=sp, cls="TOTAL_" + g.upper(),
                                 n=int((df["group"] == g).sum()),
                                 fraction=np.nan))
        comp = pd.DataFrame(rows)
        counts = {}
        for sp, df in (("human", hu), ("mouse", mm)):
            lnc = df[df["group"] == "lnc"]
            counts[sp] = {c: int((lnc["cls"] == c).sum()) for c in CLASSES}
        n_h = sum(counts["human"].values())
        n_m = sum(counts["mouse"].values())
        wide = pd.DataFrame({
            "human_n": [counts["human"][c] for c in CLASSES],
            "human_fraction": [counts["human"][c] / n_h for c in CLASSES],
            "mouse_n": [counts["mouse"][c] for c in CLASSES],
            "mouse_fraction": [counts["mouse"][c] / n_m for c in CLASSES],
        }, index=CLASSES)
        wide.index.name = "cls"
        wide.to_csv(safe_write_path(
            os.path.join(RES, "cross_species_class_comparison.tsv")), sep="\t")
        comp.to_csv(safe_write_path(
            os.path.join(RES, "cross_species_class_counts_long.tsv")),
            sep="\t", index=False)

        # chi-square on the 5 class counts
        try:
            from scipy.stats import chi2_contingency
            tab = np.array([[int(((hu[hu.group == "lnc"]["cls"] == c).sum()))
                             for c in CLASSES],
                            [int(((mm[mm.group == "lnc"]["cls"] == c).sum()))
                             for c in CLASSES]])
            chi2, p, dof, _ = chi2_contingency(tab)
            st.note("human-vs-mouse class composition chi2=" +
                    f"{chi2:.1f}" + " dof=" + str(dof) + " p=" + f"{p:.2e}")
        except Exception as e:                     # noqa: BLE001
            st.note("chi-square skipped: " + type(e).__name__)

        # activity comparison (only if both activity tables exist)
        act_rows = []
        paths = {"human": HUMAN_ACT,
                 "mouse": os.path.join(RES, "cage_activity_by_class.tsv")}
        for sp, p in paths.items():
            if not os.path.exists(p):
                continue
            a = pd.read_csv(p, sep="\t")
            if sp == "mouse" and "mode" in a.columns and \
                    not (a["mode"] == "TPM").any():
                st.note("mouse activity is peaks-only fallback; activity "
                        "comparison deferred until the TPM matrix is complete")
                continue
            a["species"] = sp
            act_rows.append(a)
        if act_rows:
            ac = pd.concat(act_rows, ignore_index=True)
            ac.to_csv(safe_write_path(
                os.path.join(RES, "cross_species_activity_comparison.tsv")),
                sep="\t", index=False)
            st.note("activity comparison rows=" + str(len(ac)))
        else:
            ac = pd.DataFrame()
            st.note("activity comparison skipped (run 'activity' for mouse)")

        body = ("## lncRNA class composition\n\n" +
                df_to_md(wide.reset_index(), ".3f") +
                "\n\n## CAGE activity by class\n\n" +
                (df_to_md(ac, ".2f") if len(ac) else
                 "_mouse activity not yet run_") + "\n")
        st.summary_md("results/step_summaries/37_mouse_compare.md",
                      "37 human vs mouse taxonomy", body)


# ------------------------------------------------------- data availability ---
DATA_AVAIL = """Data availability — mouse five-class lncRNA taxonomy
================================================================

FANTOM5 mouse CAGE (reprocessed to mm10):
  CAGE peaks: mm10_fair+new_CAGE_peaks_phase1and2.bed.gz
    https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks/mm10_fair+new_CAGE_peaks_phase1and2.bed.gz
  Per-sample TPM: mm10_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz
    https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_expression/mm10_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz
  FANTOM5 (Forrest et al. 2014, Nature 507:462-470; Lizio et al. 2015,
  Genome Biol 16:22), mm10 reprocessed release.

Mouse gene annotation:
  Ensembl Mus musculus GRCm38, release 102 (GTF)
    https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
  (Ensembl 2020; last Ensembl release on the mm10/GRCm38 assembly, chosen to
  match the mm10 CAGE peaks. Actb/Gapdh TSS coordinates verified against mm10.)

Assembly: mm10 / GRCm38.p6 (NCBI GCA_000001635.8). No liftOver chain was used;
  peaks and annotation are natively mm10.

Human comparison (for context):
  Ensembl/GENCODE v46 GRCh38 annotation and FANTOM5 hg38 CAGE, as described in
  the main methods.
"""


def cmd_data_availability():
    with step("37_mouse_data_availability") as st:
        p = safe_write_path(os.path.join(RES, "data_availability.txt"))
        with open(p, "w") as fh:
            fh.write(DATA_AVAIL)
        st.note("data-availability block -> " + p)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["prepare", "classify", "activity",
                                     "compare", "data-availability", "all"])
    args = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    if args.mode == "prepare":
        cmd_prepare()
    elif args.mode == "classify":
        cmd_classify()
    elif args.mode == "activity":
        cmd_activity()
    elif args.mode == "compare":
        cmd_compare()
    elif args.mode == "data-availability":
        cmd_data_availability()
    else:
        cmd_prepare()
        cmd_classify()
        cmd_data_availability()


if __name__ == "__main__":
    main()
