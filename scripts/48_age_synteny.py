#!/usr/bin/env python
"""48_age_synteny.py — strand-agnostic, synteny-based lncRNA age assignment.

Instead of matching TSS points (strand-biased and prone to neighbour-TSS
false positives), lift each human lncRNA gene body and require reciprocal
gene-body overlap (>=50% of both intervals) with a target-species gene.
Species: macaque, mouse, rat; deepest presence defines the age level.
Outputs -> results/comparative/E_transcript_age_synteny.tsv.
"""
import gzip
import os
import sys
from bisect import bisect_right
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR
from pylib.chainlift import ChainDB

RES = os.path.join(WORKDIR, "results")
CMP = os.path.join(WORKDIR, "data", "comparative")
OUT = os.path.join(RES, "comparative")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
SPECIES = [("macaque", 1, "chains/hg38ToRheMac10.over.chain.gz",
            "gtf/Macaca_mulatta.Mmul_10.116.gtf.gz"),
           ("mouse", 2, "chains/hg38ToMm39.over.chain.gz",
            "gtf/Mus_musculus.GRCm39.116.gtf.gz"),
           ("rat", 3, "chains/hg38ToRn7.over.chain.gz",
            "gtf/Rattus_norvegicus.mRatBN7.2.110.gtf.gz")]


def get_cands(db, chrom):
    for key in (chrom, db._canon(chrom), "chr" + db._canon(chrom)):
        if key in db.by_tname:
            return db.by_tname[key]
    return None


def open_chain(path):
    if path.endswith(".gz"):
        plain = path[:-3]
        if not os.path.exists(plain):
            with gzip.open(path, "rt") as fi, open(plain, "w") as fo:
                for line in fi:
                    fo.write(line)
        return plain
    return path


def load_genes(gtf):
    """chrom -> (starts, ends, ids, biotypes), sorted by start."""
    d = defaultdict(lambda: ([], [], [], []))
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attr = f[8]
            gid = attr.split('gene_id "')[1].split('"')[0]
            bt = ""
            for key in ("gene_biotype", "gene_type"):
                if key + ' "' in attr:
                    bt = attr.split(key + ' "')[1].split('"')[0]
                    break
            d[f[0]][0].append(int(f[3]) - 1)
            d[f[0]][1].append(int(f[4]))
            d[f[0]][2].append(gid)
            d[f[0]][3].append(bt)
    out = {}
    B = 1000000
    for c, (s, e, g, b) in d.items():
        s2, e2 = np.array(s), np.array(e)
        gg = np.array(g, dtype=object)
        bb = np.array(b, dtype=object)
        bins = defaultdict(list)
        for k in range(len(s2)):
            for bn in range(int(s2[k] // B), int((e2[k] - 1) // B) + 1):
                bins[bn].append(k)
        out[c] = (s2, e2, gg, bb, bins)
    return out


def lift_pos(chain, pos):
    i = bisect_right(chain.block_t, pos) - 1
    if i < 0:
        return None
    bt, bq, size = chain.blocks[i]
    off = pos - bt
    if off < 0 or off >= size:
        return None
    return bq + off


def lift_interval(db, chrom, s, e):
    cands = get_cands(db, chrom)
    if not cands:
        return None
    for c in cands:
        if not (c.t_start <= s and e - 1 < c.t_end):
            continue
        qs, qe = lift_pos(c, s), lift_pos(c, e - 1)
        if qs is None or qe is None:
            continue
        if c.q_strand == "+":
            return dict(chrom=c.q_name, start=qs, end=qe + 1)
        return dict(chrom=c.q_name, start=c.q_size - qe - 1,
                    end=c.q_size - qs)
    return None


def main():
    with step("48_age_synteny") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        genes = pd.read_csv(os.path.join(WORKDIR, "data", "genes",
                                         "genes.tsv"), sep="\t",
                            low_memory=False)
        genes = genes.set_index("gene_id")
        lnc = meta[meta.group == "lnc"].index
        g = genes.reindex(lnc)
        st.note("human lnc genes: {}".format(len(g)))

        res = pd.DataFrame(index=g.index)
        for sp, level, chain_rel, gtf_rel in SPECIES:
            db = ChainDB(open_chain(os.path.join(CMP, chain_rel)))
            tgt = load_genes(os.path.join(CMP, gtf_rel))
            present = []
            for gid, r in g.iterrows():
                if pd.isna(r.chr) or pd.isna(r.start):
                    present.append(0)
                    continue
                lf = lift_interval(db, str(r.chr),
                                   int(r.start) - 1, int(r.end))
                if lf is None:
                    present.append(0)
                    continue
                qc, qs, qe = lf["chrom"], lf["start"], lf["end"]
                if qe - qs < 0.5 * (r.end - r.start):
                    present.append(0)
                    continue
                arr = tgt.get(qc) or tgt.get(qc.replace("chr", ""))
                if arr is None:
                    present.append(0)
                    continue
                starts, ends, ids, bts, bins = arr
                cand = set()
                for bn in range(int(qs // 1000000),
                                int((qe - 1) // 1000000) + 1):
                    cand.update(bins.get(bn, []))
                hit = 0
                for k in cand:
                    ov = min(ends[k], qe) - max(starts[k], qs)
                    if ov <= 0:
                        continue
                    if ov >= 0.5 * (qe - qs) and ov >= 0.5 * (ends[k] -
                                                              starts[k]):
                        hit = 1
                        break
                present.append(hit)
            res[sp] = present
            st.note("{}: syntenic present={}/{} ({:.1%})".format(
                sp, int(np.sum(present)), len(present),
                np.mean(present)))
        res["age_level"] = 0
        for sp, level, _, _ in SPECIES:
            res["age_level"] = np.maximum(res["age_level"],
                                          res[sp] * level)
        lab = {0: "human_specific", 1: "primate", 2: "mammalian",
               3: "eutherian"}
        res["age_class"] = res["age_level"].map(lab)
        res.to_csv(safe_write_path(os.path.join(
            OUT, "E_transcript_age_synteny.tsv")), sep="\t")

        m = meta[["strand", "cls", "group"]].join(res)
        tab = pd.crosstab(m[m.group == "lnc"].cls,
                          m[m.group == "lnc"].age_class)
        st.note("age composition:\n" + tab.to_string())
        # strand bias check
        rows = []
        for sp, _, _, _ in SPECIES:
            a = m.loc[(m.group == "lnc") & (m.strand == "+"), sp]
            b = m.loc[(m.group == "lnc") & (m.strand == "-"), sp]
            chi2, p, _, _ = chi2_contingency([[int(a.sum()),
                                              int(len(a) - a.sum())],
                                             [int(b.sum()),
                                              int(len(b) - b.sum())]])
            rows.append(dict(species=sp, plus=float(a.mean()),
                             minus=float(b.mean()), chi2=chi2, pval=p))
        S = pd.DataFrame(rows)
        S.to_csv(safe_write_path(os.path.join(
            OUT, "age_synteny_strand_check.tsv")), sep="\t", index=False)
        st.note("strand check: " + "; ".join(
            "{} +{:.3f}/-{:.3f} p={:.1e}".format(
                r.species, r.plus, r.minus, r.pval) for r in S.itertuples()))
        st.summary_md("results/comparative/step_summaries/48_age.md",
                      "48 synteny-based age assignment",
                      "## Class x age\n\n" + df_to_md(tab.reset_index(), "d")
                      + "\n\n## Strand check\n\n" + df_to_md(S, ".3f"))


if __name__ == "__main__":
    main()
