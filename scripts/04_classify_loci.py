#!/usr/bin/env python
"""04_classify_loci.py — CORE: locus classification + CAGE-refined major TSS
+ promoter sets. Handles the special cases explicitly:

  divergent         lnc major TSS opposite-strand <=1kb from PCG major TSS
  antisense         lnc exons overlap PCG exons on opposite strand
  intronic          lnc gene body fully inside a PCG gene body (no exon overlap)
  sense_overlapping lnc exons overlap PCG exons on same strand
  intergenic        none of the above (near/far by 10kb from any PCG)

PCG and pseudogenes get their own comparator promoter sets. Every lncRNA also
carries flags: pseudogene association, lnc-lnc divergent partner, FANTOM
enhancer overlap at the promoter.

Inputs : data/genes/* (02), data/fantom/hg38_cage_peaks.bed.gz, F5.hg38.enhancers.bed.gz
Outputs: results/gene_master_table.tsv, results/classes/*.bed,
         results/promoters_{feat,fimo}.{bed,fa}, results/promoter_meta.tsv,
         results/divergent_pairs.tsv, results/exon_overlap_pairs.tsv
"""
import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR, df_to_md

DATA = os.path.join(WORKDIR, "data")
RES = os.path.join(WORKDIR, "results")
GENOME = "/media/share/resource/reference/homo_sapiens/GRCh38/genome/GRCh38.CHR.fa"
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"
CAGE_PEAKS = os.path.join(DATA, "fantom", "hg38_cage_peaks.bed.gz")
ENH = os.path.join(DATA, "fantom", "F5.hg38.enhancers.bed.gz")
MATCH_WIN = 250
DIV_WIN = 1000
INTERG_MIN = 10000
FEAT_UP, FEAT_DN = 500, 200
FIMO_UP, FIMO_DN = 500, 100


# ---------------------------------------------------------------- bedtools ---
def bt(*args):
    """bedtools via singularity, reference tree mounted read-only."""
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}",
           SIF_BT, "bedtools", *args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bedtools failed: {' '.join(args)}\n{r.stderr[:2000]}")
    return r.stdout


def bt_to_file(outfile, *args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}",
           SIF_BT, "bedtools", *args]
    with open(safe_write_path(outfile), "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bedtools failed: {' '.join(args)}\n{r.stderr[:2000]}")


def sort_bed(df, chrom_col="chrom", start_col="start"):
    c = chrom_col if chrom_col in df.columns else df.columns[0]
    s = start_col if start_col in df.columns else df.columns[1]
    return df.sort_values([c, s], kind="mergesort")


def write_bed(df, path, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    sort_bed(df, cols[0], cols[1])[cols].to_csv(safe_write_path(path), sep="\t",
                                                index=False, header=False)


# ------------------------------------------------------------- 1. load 02 ---
def load_genes():
    g = pd.read_csv(os.path.join(DATA, "genes", "genes.tsv"), sep="\t")
    t = pd.read_csv(os.path.join(DATA, "genes", "transcripts.tsv"), sep="\t")
    ex = pd.read_csv(os.path.join(DATA, "genes", "exons.bed"), sep="\t",
                     names=["chrom", "start", "end", "gid_tx", "score", "strand"])
    ex["gid"] = ex["gid_tx"].str.split("|").str[0]
    ex["tx"] = ex["gid_tx"].str.split("|").str[1]
    return g, t, ex


# --------------------------------------------- 2. CAGE peaks & refinement ---
def load_cage():
    """FANTOM fair+new CAGE peaks (defensive column parse)."""
    if not os.path.exists(CAGE_PEAKS):
        return None
    rows = []
    op = gzip.open(CAGE_PEAKS, "rt")
    for line in op:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or f[0].startswith(("#", "track", "browser")):
            continue
        chrom, s, e, pid, score, strand = f[0], int(f[1]), int(f[2]), f[3], float(f[4]), f[5]
        # summit: thickStart/thickEnd midpoint if present else interval midpoint
        if len(f) >= 8 and f[6].isdigit() and f[7].isdigit():
            summit = (int(f[6]) + int(f[7])) // 2
        else:
            summit = (s + e) // 2
        rows.append((chrom, summit, pid, score, strand))
    df = pd.DataFrame(rows, columns=["chrom", "summit", "peak_id", "score", "strand"])
    return df


def major_tss(genes, tx, cage):
    """Per gene: CAGE-refined major TSS if a same-strand peak (summit) lies
    within +/-MATCH_WIN of any transcript TSS; else best-transcript TSS.
    Best transcript priority: MANE > basic > tsl(1<2<3<NA) > longest."""
    tx = tx.copy()
    tx["tsl_n"] = pd.to_numeric(tx["tsl"], errors="coerce").fillna(9)
    tx["rank"] = (tx["mane"] * 4 + tx["basic"] * 2) * 1e9 - tx["tsl_n"] * 1e6 + tx["length"]
    best = tx.sort_values(["rank"], ascending=False).groupby("gene_id").first()
    gid2best_tss = best["tss"].to_dict()

    out = {}
    cage_by_gene = {}
    if cage is not None and len(cage):
        # candidate TSS windows (BED6: strand must be col6 for bedtools -s)
        cand = tx[["chr", "tss", "gene_id", "strand"]].copy()
        cand["sc"] = 0
        cand["w_start"] = (cand["tss"] - 1 - MATCH_WIN).clip(lower=0)
        cand["w_end"] = cand["tss"] - 1 + MATCH_WIN
        write_bed(cand, os.path.join(DATA, "tmp_cand_win.bed"),
                  ["chr", "w_start", "w_end", "gene_id", "sc", "strand"])
        # CAGE summit points (BED6)
        c = cage.copy()
        c["s"] = c["summit"] - 1
        c["e"] = c["summit"]
        cbed = c[["chrom", "s", "e", "peak_id", "score", "strand"]].copy()
        cbed["chrom"] = cbed["chrom"]
        write_bed(cbed, os.path.join(DATA, "tmp_cage_tss.bed"),
                  ["chrom", "s", "e", "peak_id", "score", "strand"])
        # intersect windows x peaks, require same strand (col6 both sides)
        pairs = bt("intersect", "-a", os.path.join(DATA, "tmp_cand_win.bed"),
                   "-b", os.path.join(DATA, "tmp_cage_tss.bed"), "-wa", "-wb", "-s")
        hits = []
        for ln in pairs.splitlines():
            f = ln.split("\t")
            if len(f) < 12:
                continue
            # a: chrom ws we gid sc strand | b: chrom s e peak_id score strand
            gid, pid, score, summit = f[3], f[9], float(f[10]), int(f[7]) + 1
            hits.append((gid, pid, score, summit))
        hd = pd.DataFrame(hits, columns=["gid", "peak_id", "score", "summit"])
        if len(hd):
            idx = hd.groupby("gid")["score"].idxmax()
            cage_by_gene = hd.loc[idx].set_index("gid")[["peak_id", "score", "summit"]].to_dict("index")

    for gid, row in genes.iterrows():
        tss0 = gid2best_tss.get(gid, row["start"] if row["strand"] == "+" else row["end"])
        c = cage_by_gene.get(gid)
        if c:
            out[gid] = (c["summit"], 1, c["peak_id"], c["score"])
        else:
            out[gid] = (tss0, 0, "", np.nan)
    return out


# ------------------------------------------------------- 3. classification ---
def classify(genes, tx, ex, tssmap, st):
    lnc = genes[genes["group"] == "lnc"]
    pcg = genes[genes["group"] == "pcg"]
    pseudo = genes[genes["group"] == "pseudo"]

    # ---- exon-level lnc x pcg overlap pairs
    lnc_ex = ex[ex["gid"].isin(lnc.index)]
    pcg_ex = ex[ex["gid"].isin(pcg.index)]
    lnc_ex = lnc_ex.copy(); lnc_ex["sc"] = 0
    pcg_ex = pcg_ex.copy(); pcg_ex["sc"] = 0
    write_bed(lnc_ex, os.path.join(DATA, "tmp_lnc_ex.bed"),
              ["chrom", "start", "end", "gid", "sc", "strand"])
    write_bed(pcg_ex, os.path.join(DATA, "tmp_pcg_ex.bed"),
              ["chrom", "start", "end", "gid", "sc", "strand"])
    ov = bt("intersect", "-a", os.path.join(DATA, "tmp_lnc_ex.bed"),
            "-b", os.path.join(DATA, "tmp_pcg_ex.bed"), "-wa", "-wb")
    pairs = set()
    strands = {}
    for ln in ov.splitlines():
        f = ln.split("\t")
        if len(f) < 12:
            continue
        a, b, sa, sb = f[3], f[9], f[5], f[11]
        pairs.add((a, b))
        strands[(a, b)] = (sa, sb)
    pd.DataFrame(sorted(pairs), columns=["lnc_gid", "pcg_gid"]).to_csv(
        safe_write_path(os.path.join(RES, "exon_overlap_pairs.tsv")), sep="\t", index=False)
    st.note(f"exon_overlap_pairs={len(pairs)}")

    # index pairs by lnc gid for O(1) lookup: gid -> [(pcg_gid, same_strand?)]
    pairs_by_lnc = {}
    for a, b in pairs:
        sa, sb = strands[(a, b)]
        pairs_by_lnc.setdefault(a, []).append((b, sa == sb))

    # ---- gene-body containment (lnc inside pcg)
    lb = lnc.reset_index().rename(columns={"chr": "c", "start": "s", "end": "e"})
    lb["s0"] = lb["s"] - 1
    lb["sc"] = 0
    write_bed(lb, os.path.join(DATA, "tmp_lnc_g.bed"), ["c", "s0", "e", "gene_id", "sc", "strand"])
    pb = pcg.reset_index().rename(columns={"chr": "c", "start": "s", "end": "e"})
    pb["s0"] = pb["s"] - 1
    pb["sc"] = 0
    write_bed(pb, os.path.join(DATA, "tmp_pcg_g.bed"), ["c", "s0", "e", "gene_id", "sc", "strand"])
    cont = bt("intersect", "-a", os.path.join(DATA, "tmp_lnc_g.bed"),
              "-b", os.path.join(DATA, "tmp_pcg_g.bed"), "-wa", "-wb")
    contain = {}
    for ln in cont.splitlines():
        f = ln.split("\t")
        if len(f) < 12:
            continue
        a, b = f[3], f[9]
        if int(f[1]) >= int(f[7]) and int(f[2]) <= int(f[8]):  # a within b
            contain.setdefault(a, b)
    st.note(f"contained_in_pcg={len(contain)}")

    # ---- TSS-based distances (divergent)
    def tss_bed(df, path):
        d = df.reset_index()
        d["tss"] = [tssmap[g][0] for g in d["gene_id"]]
        d["s0"] = d["tss"] - 1
        d["e0"] = d["tss"]
        d["sc"] = 0
        write_bed(d, path, ["chr", "s0", "e0", "gene_id", "sc", "strand"])

    tss_bed(lnc, os.path.join(DATA, "tmp_lnc_tss.bed"))
    tss_bed(pcg, os.path.join(DATA, "tmp_pcg_tss.bed"))
    tss_bed(pseudo, os.path.join(DATA, "tmp_pseudo_tss.bed"))

    cl = bt("closest", "-a", os.path.join(DATA, "tmp_lnc_tss.bed"),
            "-b", os.path.join(DATA, "tmp_pcg_tss.bed"), "-S", "-d", "-t", "first")
    div, div_str = {}, {}
    for ln in cl.splitlines():
        f = ln.split("\t")
        if len(f) < 13 or f[6] == ".":
            continue
        gid, partner, dist = f[3], f[9], int(f[12])
        div[gid] = dist
        div_str[gid] = partner
    # same-strand nearest PCG TSS distance (context)
    cs = bt("closest", "-a", os.path.join(DATA, "tmp_lnc_tss.bed"),
            "-b", os.path.join(DATA, "tmp_pcg_tss.bed"), "-s", "-d", "-t", "first")
    same_d = {}
    for ln in cs.splitlines():
        f = ln.split("\t")
        if len(f) < 13 or f[6] == ".":
            continue
        same_d[f[3]] = int(f[12])
    # nearest PCG gene body distance (intergenic near/far)
    gb = bt("closest", "-a", os.path.join(DATA, "tmp_lnc_g.bed"),
            "-b", os.path.join(DATA, "tmp_pcg_g.bed"), "-d", "-t", "first")
    body_d = {}
    for ln in gb.splitlines():
        f = ln.split("\t")
        if len(f) < 13 or f[6] == ".":
            continue
        body_d[f[3]] = int(f[12])

    # ---- lnc-lnc divergent (for lncs NOT divergent vs PCG)
    ll = bt("closest", "-a", os.path.join(DATA, "tmp_lnc_tss.bed"),
            "-b", os.path.join(DATA, "tmp_lnc_tss.bed"), "-S", "-d", "-t", "first")
    ll_div = {}
    for ln in ll.splitlines():
        f = ln.split("\t")
        if len(f) < 13 or f[6] == ".":
            continue
        a, b, dist = f[3], f[9], int(f[12])
        if a != b and dist <= DIV_WIN:
            ll_div[a] = (b, dist)

    # ---- pseudogene association
    pe = bt("intersect", "-a", os.path.join(DATA, "tmp_lnc_g.bed"),
            "-b", os.path.join(DATA, "tmp_pseudo_tss.bed"), "-wa", "-wb", "-u")
    # simpler: gene body overlap with any pseudogene gene body:
    pdb = pseudo.reset_index().rename(columns={"chr": "c", "start": "s", "end": "e"})
    pdb["s0"] = pdb["s"] - 1
    pdb["sc"] = 0
    write_bed(pdb, os.path.join(DATA, "tmp_pseudo_g.bed"), ["c", "s0", "e", "gene_id", "sc", "strand"])
    ov2 = bt("intersect", "-a", os.path.join(DATA, "tmp_lnc_g.bed"),
             "-b", os.path.join(DATA, "tmp_pseudo_g.bed"), "-wa", "-wb")
    pseudo_ov = {}
    for ln in ov2.splitlines():
        f = ln.split("\t")
        if len(f) >= 12:
            pseudo_ov.setdefault(f[3], f[9])
    pc = bt("closest", "-a", os.path.join(DATA, "tmp_lnc_tss.bed"),
            "-b", os.path.join(DATA, "tmp_pseudo_tss.bed"), "-d", "-t", "first")
    pseudo_tss_d = {}
    for ln in pc.splitlines():
        f = ln.split("\t")
        if len(f) >= 13 and f[6] != ".":
            pseudo_tss_d[f[3]] = int(f[12])

    # ---- assemble classes (priority order)
    recs = []
    div_rows = []
    for gid, row in lnc.iterrows():
        tss, refined, pid, pscore = tssmap[gid]
        # 1 divergent
        if gid in div and div[gid] <= DIV_WIN:
            cls, partner = "divergent", div_str[gid]
            div_rows.append((gid, partner, div[gid], row["strand"],
                             genes.loc[partner, "gene_name"] if partner in genes.index else "?"))
        else:
            plist = pairs_by_lnc.get(gid, [])
            anti = [b for b, same in plist if not same]
            sense = [b for b, same in plist if same]
            if anti:                                   # 2 antisense
                cls, partner = "antisense", anti[0]
            elif gid in contain:                       # 3 intronic
                cls, partner = "intronic", contain[gid]
            elif sense:                                # 4 sense_overlapping
                cls, partner = "sense_overlapping", sense[0]
            else:
                cls, partner = "intergenic", ""
        recs.append(dict(
            gene_id=gid, gene_name=row["gene_name"], group="lnc", cls=cls,
            chrom=row["chr"], tss=tss, strand=row["strand"],
            cage_refined=refined, cage_peak_id=pid, cage_peak_score=pscore,
            partner_pcg=partner,
            div_dist_opp=div.get(gid), same_tss_dist=same_d.get(gid),
            pcg_body_dist=body_d.get(gid),
            lnc_lnc_div_partner=ll_div.get(gid, ("", 0))[0],
            lnc_lnc_div_dist=ll_div.get(gid, ("", 0))[1],
            pseudo_assoc=1 if (gid in pseudo_ov or pseudo_tss_d.get(gid, 1e9) <= 2000) else 0,
            pseudo_overlap_gid=pseudo_ov.get(gid, ""),
            pseudo_tss_dist=pseudo_tss_d.get(gid),
            transcribed_pseudo=0))

    for gid, row in pcg.iterrows():
        tss, refined, pid, pscore = tssmap[gid]
        recs.append(dict(gene_id=gid, gene_name=row["gene_name"], group="pcg", cls="PCG",
                         chrom=row["chr"], tss=tss, strand=row["strand"],
                         cage_refined=refined, cage_peak_id=pid, cage_peak_score=pscore,
                         partner_pcg="", transcribed_pseudo=0))
    for gid, row in pseudo.iterrows():
        tss, refined, pid, pscore = tssmap[gid]
        recs.append(dict(gene_id=gid, gene_name=row["gene_name"], group="pseudo",
                         cls="transcribed_pseudogene" if row["transcribed_pseudo"] else "pseudogene",
                         chrom=row["chr"], tss=tss, strand=row["strand"],
                         cage_refined=refined, cage_peak_id=pid, cage_peak_score=pscore,
                         partner_pcg="", transcribed_pseudo=int(row["transcribed_pseudo"])))

    m = pd.DataFrame(recs).set_index("gene_id")
    # divergent sensitivity flags
    m["div_500"] = [(1 if (not np.isnan(d) and d <= 500) else 0) if g == "lnc" else np.nan
                    for d, g in zip(m["div_dist_opp"].astype(float), m["group"])]
    m["div_2000"] = [(1 if (not np.isnan(d) and d <= 2000) else 0) if g == "lnc" else np.nan
                     for d, g in zip(m["div_dist_opp"].astype(float), m["group"])]
    m["intergenic_near"] = [1 if (g == "intergenic" and (b is None or b <= INTERG_MIN)) else
                            (0 if g == "intergenic" else np.nan)
                            for g, b in zip(m["cls"], m["pcg_body_dist"])]
    pd.DataFrame(div_rows, columns=["lnc_gid", "pcg_gid", "distance", "lnc_strand", "pcg_name"]).to_csv(
        safe_write_path(os.path.join(RES, "divergent_pairs.tsv")), sep="\t", index=False)
    return m


# ------------------------------------------------------- 4. enhancer flag ----
def enhancer_flag(meta, st):
    if not os.path.exists(ENH):
        meta["enhancer_overlap"] = np.nan
        st.note("enhancer file missing — flag = NA")
        return meta
    d = meta.reset_index()
    d["s0"] = [t - (FEAT_UP if s == "+" else FEAT_DN) - 1 for t, s in zip(d["tss"], d["strand"])]
    d["e0"] = [t + (FEAT_DN if s == "+" else FEAT_UP) - 1 for t, s in zip(d["tss"], d["strand"])]
    d["sc"] = 0
    write_bed(d, os.path.join(DATA, "tmp_prom_win.bed"), ["chrom", "s0", "e0", "gene_id", "sc", "strand"])
    ov = bt("intersect", "-a", os.path.join(DATA, "tmp_prom_win.bed"), "-b", ENH, "-u")
    hit = {ln.split("\t")[3] for ln in ov.splitlines() if ln}
    meta["enhancer_overlap"] = [1 if g in hit else 0 for g in meta.index]
    st.note(f"enhancer_overlap={int(meta['enhancer_overlap'].sum())} promoters")
    return meta


# ------------------------------------------------------- 5. promoter sets ----
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def make_promoters(meta):
    import pysam
    fa = pysam.FastaFile(GENOME)
    chrom_len = dict(zip(fa.references, fa.lengths))
    d = meta.reset_index()
    n_missing = 0

    for tag, (up, dn) in (("feat", (FEAT_UP, FEAT_DN)), ("fimo", (FIMO_UP, FIMO_DN))):
        bed_rows, seqs = [], []
        for _, r in d.iterrows():
            tss, strand, chrom = int(r["tss"]), r["strand"], r["chrom"]
            if strand == "+":
                s0, e0 = tss - up - 1, tss + dn - 1
            else:
                s0, e0 = tss - dn, tss + up
            if chrom not in chrom_len:                 # scaffold absent from CHR.fa
                n_missing += 1
                bed_rows.append((chrom, max(0, s0), e0, r["gene_id"], 0, strand))
                seqs.append((r["gene_id"], ""))
                continue
            s0c, e0c = max(0, s0), min(e0, chrom_len[chrom])
            seq = fa.fetch(chrom, s0c, e0c).upper() if s0c < e0c else ""
            if strand == "-":
                seq = revcomp(seq)
            bed_rows.append((chrom, s0c, e0c, r["gene_id"], 0, strand))
            seqs.append((r["gene_id"], seq))
        bedp = safe_write_path(os.path.join(RES, f"promoters_{tag}.bed"))
        with open(bedp, "w") as fh:
            for c, s, e, gid, sc, st_ in sorted(bed_rows, key=lambda x: (x[0], x[1])):
                fh.write(f"{c}\t{s}\t{e}\t{gid}\t0\t{st_}\n")
        fap = safe_write_path(os.path.join(RES, f"promoters_{tag}.fa"))
        with open(fap, "w") as fh:
            for gid, seq in seqs:
                fh.write(f">{gid}\n{seq}\n")
    fa.close()
    return n_missing


# ------------------------------------------------------------------- main ---
def main():
    with step("04_classify_loci") as st:
        genes, tx, ex = load_genes()
        genes = genes.set_index("gene_id")
        st.note(f"genes={len(genes)}")
        cage = load_cage()
        st.note(f"cage_peaks={0 if cage is None else len(cage)}")
        tssmap = major_tss(genes, tx, cage)
        n_ref = sum(1 for v in tssmap.values() if v[1] == 1)
        st.note(f"cage_refined={n_ref} ({n_ref/len(tssmap):.1%})")
        meta = classify(genes, tx, ex, tssmap, st)
        meta = enhancer_flag(meta, st)

        # class distribution + asserts
        dist = meta[meta["group"] == "lnc"]["cls"].value_counts()
        st.note("lnc classes: " + ", ".join(f"{k}={v}" for k, v in dist.items()))
        assert int(dist.sum()) == 19258, f"class sum {dist.sum()} != 19258"
        pcg_ref = meta[(meta.group == "pcg") & (meta.cage_refined == 1)].shape[0]
        st.note(f"PCG cage_refined={pcg_ref/ (meta.group=='pcg').sum():.1%}")

        meta.to_csv(safe_write_path(os.path.join(RES, "gene_master_table.tsv")), sep="\t")
        # per-class BEDs
        os.makedirs(os.path.join(RES, "classes"), exist_ok=True)
        for cls, sub in meta[meta["group"] == "lnc"].groupby("cls"):
            d = sub.reset_index()
            d["s0"], d["e0"] = d["tss"] - 1, d["tss"]
            write_bed(d, os.path.join(RES, "classes", f"{cls}.bed"),
                      ["chrom", "s0", "e0", "gene_id", "strand"])
        n_miss = make_promoters(meta)
        st.note(f"promoter sets written (feat + fimo); scaffold-missing-seq genes={n_miss}")

        # named-gene QC print
        named = meta[meta["gene_name"].isin(
            ["XIST", "TSIX", "GAS5", "NEAT1", "MALAT1", "H19", "SNHG1", "LINC-PINT", "PEG3"])]
        for _, r in named.iterrows():
            st.note(f"QC {r['gene_name']}: group={r['group']} class={r['cls']} refined={r['cage_refined']}")

        body = ("| class | n |\n|---|---|\n" +
                "\n".join(f"| {k} | {v} |" for k, v in dist.items()) +
                f"\n\nCAGE-refined TSS: PCG {pcg_ref/(meta.group=='pcg').sum():.1%}, "
                f"lnc {meta[(meta.group=='lnc') & (meta.cage_refined==1)].shape[0]/19258:.1%}\n")
        st.summary_md("results/step_summaries/04_classify_loci.md", "04 locus classification", body)


if __name__ == "__main__":
    main()
