#!/usr/bin/env python
"""59_cross_species_taxonomy.py — five-class taxonomy for non-human species.

The human classes (scripts/04) use a CAGE-refined major TSS. No comparable CAGE
peak set exists for the non-human scBaseCount samples, so this module uses the
coordinate rules alone, with the major TSS taken from the best-supported
transcript of each gene (canonical and basic first, then transcript support,
then length). The same simplification was validated on mouse, where the
CAGE-refined and purely coordinate TSS assignments give the same class for
94.6% of lncRNAs and the divergent fraction moves by 1.4 points.

Classes use the human priority order so that class definitions are comparable
  divergent         lnc TSS opposite-strand within 1 kb of a PCG TSS
  antisense         lnc exons overlap PCG exons on the opposite strand
  intronic          lnc body inside a PCG body with no exon overlap
  sense_overlapping lnc exons overlap PCG exons on the same strand
  intergenic        none of the above

CLI
  prepare <species>   parse the Ensembl GTF into one table per species
  classify <species>  assign the five classes
  all                 every species in SPECIES

Outputs
  results/cross_species/<species>/gene_table.tsv
  results/cross_species/<species>/gene_master_table.tsv
  results/cross_species/taxonomy_summary.tsv
  results/step_summaries/59_cross_species.md
"""
import argparse
import gzip
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

DATA = os.path.join(WORKDIR, "data")
OUT = os.path.join(WORKDIR, "results", "cross_species")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"
DIV_WIN = 1000
INTERG_MIN = 10000

# Ensembl release 116 except where the atlas reference is fixed older
SPECIES = {
    "macaque": ("Macaca_mulatta.Mmul_10.116", "Mmul_10"),
    "rat": ("Rattus_norvegicus.mRatBN7.2.110", "mRatBN7.2"),
    "rabbit": ("Oryctolagus_cuniculus.OryCun2.0.116", "OryCun2.0"),
    "opossum": ("Monodelphis_domestica.ASM229v1.116", "ASM229v1"),
    "mouse": ("Mus_musculus.GRCm38.102", "GRCm38"),
}
LNC_TYPES = {
    "lncRNA", "lincRNA", "antisense", "processed_transcript", "sense_intronic",
    "sense_overlapping", "bidirectional_promoter_lncRNA", "macro_lncRNA",
    "3prime_overlapping_ncRNA", "non_coding",
}
_ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
_TAG_RE = re.compile(r'tag "([^"]*)"')


def group_of(biotype):
    if biotype == "protein_coding":
        return "pcg"
    if biotype in LNC_TYPES:
        return "lnc"
    if "pseudogene" in biotype:
        return "pseudo"
    return "other"


def norm_chrom(c):
    if c in ("MT", "M"):
        return "chrM"
    return c if c.startswith("chr") else "chr" + c


def bt(*args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", WORKDIR + ":" + WORKDIR, SIF_BT, "bedtools", *args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("bedtools failed: " + " ".join(args) + "\n"
                           + r.stderr[:2000])
    return r.stdout


def tmp(species, name):
    d = os.path.join(DATA, "cross_species_tmp", species)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


def write_bed(df, path, cols):
    d = df.sort_values([cols[0], cols[1]], kind="mergesort")
    d[cols].to_csv(safe_write_path(path), sep="\t", index=False, header=False)


# --------------------------------------------------------------- prepare ---
def cmd_prepare(sp):
    gtf_rel, _ = SPECIES[sp]
    gtf = os.path.join(DATA, "comparative", "gtf", gtf_rel + ".gtf.gz")
    out = os.path.join(OUT, sp)
    os.makedirs(out, exist_ok=True)
    with step("59_{}_prepare".format(sp)) as st:
        genes, txs, exons = [], [], []
        with gzip.open(gtf, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                feat = f[2]
                if feat not in ("gene", "transcript", "exon"):
                    continue
                a = _ATTR_RE.findall(f[8])
                a = dict(a) if a else {}
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
                    tags = set(_TAG_RE.findall(f[8]))
                    tss = start if strand == "+" else end
                    txs.append(dict(
                        gene_id=gid,
                        transcript_id=a.get("transcript_id", ""),
                        chrom=chrom, start=start, end=end, strand=strand,
                        tss=tss, length=end - start + 1,
                        tsl=a.get("transcript_support_level", "NA"),
                        basic=int(bool({"basic", "gencode_basic"} & tags)),
                        canonical=int(bool({"Ensembl_canonical", "canonical",
                                            "MANE_Select"} & tags))))
                else:
                    exons.append((chrom, start - 1, end,
                                  gid + "|" + a.get("transcript_id", ""),
                                  0, strand))
        st.note("GTF parsed: genes={} transcripts={} exons={}".format(
            len(genes), len(txs), len(exons)))

        g = pd.DataFrame(genes).set_index("gene_id")
        g["tss"] = np.where(g["strand"] == "+", g["start"], g["end"])
        g["length"] = g["end"] - g["start"] + 1
        g["group"] = g["gene_type"].map(group_of)
        tx = pd.DataFrame(txs)
        ex = pd.DataFrame(exons, columns=["chrom", "start", "end", "gid_tx",
                                          "score", "strand"])
        ex["gene_id"] = ex["gid_tx"].str.split("|").str[0]
        g["n_exons"] = g.index.map(
            ex.groupby("gene_id").size()).fillna(0).astype(int)
        g.to_csv(safe_write_path(os.path.join(out, "gene_table.tsv")),
                 sep="\t")
        tx.to_csv(safe_write_path(os.path.join(out, "transcripts.tsv")),
                  sep="\t", index=False)
        ex[["chrom", "start", "end", "gid_tx", "score", "strand"]].to_csv(
            safe_write_path(os.path.join(out, "exons.bed.gz")), sep="\t",
            index=False, header=False, compression="gzip")
        counts = g["group"].value_counts()
        st.note("groups: " + ", ".join(
            "{}={}".format(k, v) for k, v in counts.items()))


# -------------------------------------------------------------- classify ---
def cmd_classify(sp):
    out = os.path.join(OUT, sp)
    with step("59_{}_classify".format(sp)) as st:
        g = pd.read_csv(os.path.join(out, "gene_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        tx = pd.read_csv(os.path.join(out, "transcripts.tsv"), sep="\t")
        ex = pd.read_csv(os.path.join(out, "exons.bed.gz"), sep="\t",
                         header=None, compression="gzip",
                         names=["chrom", "start", "end", "gid_tx", "score",
                                "strand"])
        ex["gene_id"] = ex["gid_tx"].str.split("|").str[0]

        # Best transcript per gene. Without CAGE this is the closest available
        # proxy for the major TSS.
        tx["tsl_n"] = pd.to_numeric(tx["tsl"], errors="coerce").fillna(9)
        tx["rank"] = ((tx["canonical"] * 4 + tx["basic"] * 2) * 1e9
                      - tx["tsl_n"] * 1e6 + tx["length"])
        best = tx.sort_values("rank", ascending=False).groupby("gene_id").first()
        tss = best["tss"].reindex(g.index).fillna(g["tss"])
        g["tss_major"] = tss.values
        st.note("major TSS from a canonical or basic transcript: {}/{}".format(
            int(best["canonical"].sum() + best["basic"].sum()), len(g)))

        lnc = g[g.group == "lnc"]
        pcg = g[g.group == "pcg"]
        pseudo = g[g.group == "pseudo"]
        st.note("lnc={} pcg={} pseudo={}".format(len(lnc), len(pcg),
                                                 len(pseudo)))

        ex["gid"] = ex["gene_id"]
        lx = ex[ex.gid.isin(lnc.index)]
        px = ex[ex.gid.isin(pcg.index)]
        write_bed(lx, tmp(sp, "lnc_ex.bed"),
                  ["chrom", "start", "end", "gid", "score", "strand"])
        write_bed(px, tmp(sp, "pcg_ex.bed"),
                  ["chrom", "start", "end", "gid", "score", "strand"])
        ov = bt("intersect", "-a", tmp(sp, "lnc_ex.bed"),
                "-b", tmp(sp, "pcg_ex.bed"), "-wa", "-wb")
        pairs, strands = set(), {}
        for ln in ov.splitlines():
            f = ln.split("\t")
            if len(f) < 12:
                continue
            key = (f[3], f[9])
            pairs.add(key)
            strands[key] = (f[5], f[11])
        pairs_by_lnc = {}
        for a, b in pairs:
            sa, sb = strands[(a, b)]
            pairs_by_lnc.setdefault(a, []).append((b, sa == sb))
        st.note("exon_overlap_pairs={}".format(len(pairs)))

        # ---- gene body containment
        for name, df in (("lnc", lnc), ("pcg", pcg)):
            d = df.reset_index().rename(columns={"chrom": "c", "start": "s",
                                                 "end": "e"})
            d["s0"] = d["s"] - 1
            d["sc"] = 0
            write_bed(d, tmp(sp, name + "_g.bed"),
                      ["c", "s0", "e", "gene_id", "sc", "strand"])
        cont = bt("intersect", "-a", tmp(sp, "lnc_g.bed"),
                  "-b", tmp(sp, "pcg_g.bed"), "-wa", "-wb")
        contain = {}
        for ln in cont.splitlines():
            f = ln.split("\t")
            if len(f) < 12:
                continue
            if int(f[1]) >= int(f[7]) and int(f[2]) <= int(f[8]):
                contain.setdefault(f[3], f[9])
        st.note("contained_in_pcg={}".format(len(contain)))

        # ---- TSS sets
        def tss_bed(df, path):
            d = df.reset_index()
            d["s0"] = d["tss_major"] - 1
            d["e0"] = d["tss_major"]
            d["sc"] = 0
            write_bed(d, path,
                      ["chrom", "s0", "e0", "gene_id", "sc", "strand"])

        tss_bed(lnc, tmp(sp, "lnc_tss.bed"))
        tss_bed(pcg, tmp(sp, "pcg_tss.bed"))
        cl = bt("closest", "-a", tmp(sp, "lnc_tss.bed"),
                "-b", tmp(sp, "pcg_tss.bed"), "-S", "-d", "-t", "first")
        div, div_str = {}, {}
        for ln in cl.splitlines():
            f = ln.split("\t")
            if len(f) < 13 or f[6] == ".":
                continue
            div[f[3]] = int(f[12])
            div_str[f[3]] = f[9]

        # ---- assemble
        recs, drows = [], []
        for gid, row in lnc.iterrows():
            if gid in div and div[gid] <= DIV_WIN:
                cls, partner = "divergent", div_str[gid]
                nm = (g.at[partner, "gene_name"] if partner in g.index else "?")
                drows.append(dict(lnc_gid=gid, pcg_gid=partner,
                                  distance=div[gid], lnc_strand=row["strand"],
                                  pcg_name=nm))
            else:
                pl = pairs_by_lnc.get(gid, [])
                anti = [b for b, same in pl if not same]
                sense = [b for b, same in pl if same]
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
                chrom=row["chrom"], tss=int(row["tss_major"]),
                strand=row["strand"], partner_pcg=partner))
        L = pd.DataFrame(recs)
        pd.DataFrame(drows).to_csv(
            safe_write_path(os.path.join(out, "divergent_pairs.tsv")),
            sep="\t", index=False)
        P = g[g.group.isin(["pcg", "pseudo", "other"])].reset_index()
        P["cls"] = P["group"].map({"pcg": "PCG", "pseudo": "pseudogene"}
                                  ).fillna("other")
        # the gene table already carries a tss column, so the CAGE-free major
        # TSS has to replace it rather than be renamed onto it
        P = P.drop(columns=["tss"]).rename(columns={"tss_major": "tss"})[
            ["gene_id", "gene_name", "gene_type", "group", "cls", "chrom",
             "tss", "strand", "n_exons", "length"]]
        L["n_exons"] = L["gene_id"].map(g["n_exons"])
        L["length"] = L["gene_id"].map(g["length"])
        M = pd.concat([L, P[["gene_id", "gene_name", "gene_type", "group",
                             "cls", "chrom", "tss", "strand", "n_exons",
                             "length"]]], ignore_index=True)
        M.to_csv(safe_write_path(os.path.join(
            out, "gene_master_table.tsv")), sep="\t", index=False)
        cnt = L.cls.value_counts()
        st.note("classes: " + ", ".join(
            "{}={}".format(k, v) for k, v in cnt.items()))
        st.note("divergent fraction of lncRNA = {:.1%}".format(
            cnt.get("divergent", 0) / max(len(L), 1)))
        st.note("divergent pairs={}".format(len(drows)))


def summarise():
    rows = []
    for sp in SPECIES:
        p = os.path.join(OUT, sp, "gene_master_table.tsv")
        if not os.path.exists(p):
            continue
        m = pd.read_csv(p, sep="\t", low_memory=False)
        lnc = m[m.group == "lnc"]
        row = dict(species=sp, lnc=len(lnc), pcg=int((m.group == "pcg").sum()))
        for c in ("divergent", "antisense", "intronic", "sense_overlapping",
                  "intergenic"):
            row[c] = int((lnc.cls == c).sum())
            row[c + "_pct"] = round(100 * (lnc.cls == c).mean(), 1) if len(lnc) else np.nan
        rows.append(row)
    S = pd.DataFrame(rows)
    S.to_csv(safe_write_path(os.path.join(OUT, "taxonomy_summary.tsv")),
             sep="\t", index=False)
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["prepare", "classify", "all"])
    ap.add_argument("species", nargs="*", default=[])
    args = ap.parse_args()
    targets = args.species or list(SPECIES)
    if args.mode in ("prepare", "all"):
        for sp in targets:
            cmd_prepare(sp)
    if args.mode in ("classify", "all"):
        for sp in targets:
            cmd_classify(sp)
        S = summarise()
        with step("59_cross_species_summary") as st:
            st.summary_md("results/step_summaries/59_cross_species.md",
                          "59 cross-species taxonomy", df_to_md(S, ".3g"))


if __name__ == "__main__":
    main()
