#!/usr/bin/env python
"""41_syntenic_lnc.py — syntenic cross-species mapping of human lncRNAs.

BioMart/Compara orthologs are nearly absent for lncRNAs (0-5 per class).
Here we map each human lncRNA TSS to a species lncRNA TSS by chain liftover
(same strand, <=1 kb), then compare developmental expression breadth
(Cardoso-Moreira) between human and the syntenic species gene.

Modes: species key(s): mouse | macaque | rat | rabbit | chicken | opossum
Outputs -> results/syntenic/*.tsv + step summary.
"""
import gzip
import os
import sys
from bisect import bisect_right
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR
from pylib.chainlift import ChainDB

RES = os.path.join(WORKDIR, "results")
CMP = os.path.join(WORKDIR, "data", "comparative")
OUT = os.path.join(RES, "syntenic")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
SPECIES = {
    "mouse": ("gtf/Mus_musculus.GRCm39.116.gtf.gz",
              "chains/hg38ToMm39.over.chain.gz"),
    "rat": ("gtf/Rattus_norvegicus.mRatBN7.2.110.gtf.gz",
            "chains/hg38ToRn7.over.chain.gz"),
    "macaque": ("gtf/Macaca_mulatta.Mmul_10.116.gtf.gz",
                "chains/hg38ToRheMac10.over.chain.gz"),
    "rabbit": ("gtf/Oryctolagus_cuniculus.OryCun2.0.116.gtf.gz",
               "chains/hg38ToOryCun2.over.chain.gz"),
    "chicken": ("gtf/Gallus_gallus.GRCg6a.104.gtf.gz",
                "chains/hg38ToGalGal6.over.chain.gz"),
    "opossum": ("gtf/Monodelphis_domestica.ASM229v1.116.gtf.gz",
                "chains/hg38ToMonDom5.over.chain.gz"),
}
LNC_BIOTYPES = {"lncRNA", "lincRNA", "antisense", "processed_transcript",
                "sense_intronic", "sense_overlapping", "macro_lncRNA",
                "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA",
                "non_coding", "sense_overlapping"}


def get_cands(db, chrom):
    for key in (chrom, db._canon(chrom), "chr" + db._canon(chrom)):
        if key in db.by_tname:
            return db.by_tname[key]
    return None


def lift_point(db, chrom, pos):
    cands = get_cands(db, chrom)
    if not cands:
        return None
    for c in cands:
        if not (c.t_start <= pos < c.t_end):
            continue
        i = bisect_right(c.block_t, pos) - 1
        if i < 0:
            continue
        bt, bq, size = c.blocks[i]
        off = pos - bt
        if off < 0 or off >= size:
            continue
        qpos = (bq + off) if c.q_strand == "+" else \
            c.q_size - (bq + off) - 1
        return c.q_name, qpos, c.q_strand
    return None


def gtf_key(d, name):
    for key in (name, name.replace("chr", ""), "chr" + name.replace("chr", "")):
        if key in d:
            return key
    return None


def parse_gtf_lnc(path):
    d = defaultdict(lambda: {"+": [], "-": []})
    gid = defaultdict(dict)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            attr = p[8]
            bt = None
            for key in ("gene_biotype", "gene_type"):
                if key + ' "' in attr:
                    bt = attr.split(key + ' "')[1].split('"')[0]
                    break
            if bt is None or bt not in LNC_BIOTYPES:
                continue
            chrom = p[0]
            start, end, strand = int(p[3]), int(p[4]), p[6]
            g = attr.split('gene_id "')[1].split('"')[0]
            tss = start if strand == "+" else end
            d[chrom][strand].append(tss)
            gid[chrom][tss] = g
    for c in d:
        for s in ("+", "-"):
            d[c][s] = np.array(sorted(d[c][s]))
    return d, gid


def main():
    mode = (sys.argv[1] if len(sys.argv) > 1 else "mouse").lower()
    if mode not in SPECIES:
        raise SystemExit("usage: 41_syntenic_lnc.py " + "|".join(SPECIES))
    gtf_rel, chain_rel = SPECIES[mode]
    gtf = os.path.join(CMP, gtf_rel)
    chain = os.path.join(CMP, chain_rel)
    if not (os.path.exists(gtf) and os.path.exists(chain)):
        raise RuntimeError("missing {} or {}".format(gtf, chain))
    with step("41_syntenic_" + mode) as st:
        # chain may be gz; ChainDB needs plain text
        cp = chain
        if cp.endswith(".gz"):
            plain = cp[:-3]
            if not os.path.exists(plain):
                with gzip.open(cp, "rt") as fi, open(plain, "w") as fo:
                    for line in fi:
                        fo.write(line)
            cp = plain
        db = ChainDB(cp)
        tss_db, gid_db = parse_gtf_lnc(gtf)
        st.note("{}: lnc genes parsed={}".format(
            mode, sum(len(v["+"]) + len(v["-"]) for v in tss_db.values())))

        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        lnc = meta[meta.group == "lnc"][["chrom", "tss", "strand", "cls"]]
        # human breadth
        hs = pd.read_csv(os.path.join(RES, "developmental",
                                      "per_gene_summary.tsv"), sep="\t")
        hs = hs[hs.species == "human"][["gene", "breadth_frac"]]
        hs["base"] = hs.gene.str.split(".").str[0]
        hs = hs.drop_duplicates("base").set_index("base")
        # species breadth
        ss = pd.read_csv(os.path.join(RES, "developmental",
                                      "per_gene_summary.tsv"), sep="\t")
        ss = ss[ss.species == mode][["gene", "breadth_frac"]]
        ss["gene"] = ss.gene.str.split(".").str[0]
        ss = ss.drop_duplicates("gene").set_index("gene")

        rows = []
        for gid, r in lnc.iterrows():
            if pd.isna(r.tss) or pd.isna(r.chrom):
                continue
            pos = int(r.tss) - 1
            lf = lift_point(db, str(r.chrom), pos)
            if lf is None:
                continue
            qname, qpos, qstrand = lf
            strand = "+" if qstrand == "+" else "-"
            k = gtf_key(tss_db, qname)
            arr = tss_db.get(k, {}).get(strand) if k else None
            if arr is None or len(arr) == 0:
                continue
            i = np.searchsorted(arr, qpos)
            hit = None
            for j in (i - 1, i):
                if 0 <= j < len(arr) and abs(int(arr[j]) - qpos) <= 1000:
                    hit = int(arr[j])
                    break
            if hit is None:
                continue
            gk = gtf_key(gid_db, qname)
            sgid = gid_db.get(gk, {}).get(hit) if gk else None
            if sgid is None:
                continue
            sgid = sgid.split(".")[0]
            hb = hs.breadth_frac.get(gid.split(".")[0], np.nan)
            sb = ss.breadth_frac.get(sgid, np.nan)
            rows.append(dict(human_gene=gid, cls=r.cls, species_gene=sgid,
                             human_breadth=hb, species_breadth=sb))
        M = pd.DataFrame(rows, columns=["human_gene", "cls", "species_gene",
                                        "human_breadth", "species_breadth"])
        M.to_csv(safe_write_path(os.path.join(
            OUT, "map_{}.tsv".format(mode))), sep="\t", index=False)
        st.note("{}: mapped lnc={} (of {})".format(mode, len(M), len(lnc)))

        trows = []
        for cls in CLASSES:
            sub = M[(M.cls == cls) & M.human_breadth.notna() &
                    M.species_breadth.notna()]
            if len(sub) >= 20:
                r, p = spearmanr(sub.human_breadth, sub.species_breadth)
                trows.append(dict(species=mode, cls=cls, n=len(sub),
                                  rho=float(r), pval=float(p),
                                  median_human=float(sub.human_breadth.median()),
                                  median_species=float(
                                      sub.species_breadth.median())))
            else:
                trows.append(dict(species=mode, cls=cls, n=len(sub),
                                  rho=np.nan, pval=np.nan,
                                  median_human=np.nan, median_species=np.nan))
        T = pd.DataFrame(trows)
        T.to_csv(safe_write_path(os.path.join(
            OUT, "rho_{}.tsv".format(mode))), sep="\t", index=False)
        st.note("rho: " + "; ".join(
            "{} n={} rho={:.2f}".format(r.cls, r.n, r.rho)
            for r in T.itertuples()))
        st.summary_md("results/syntenic/step_summaries/41_{}.md".format(mode),
                      "41 syntenic lnc mapping ({})".format(mode),
                      "## Mapped lncRNAs by class\n\n"
                      + df_to_md(T, ".3f"))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    main()
