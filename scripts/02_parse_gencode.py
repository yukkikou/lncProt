#!/usr/bin/env python
"""02_parse_gencode.py — stream GENCODE v46 GTF into compact tables.

Outputs (data/genes/):
  genes.tsv        gene_id, gene_type, gene_name, chr, start, end, strand, level
  transcripts.tsv  gene_id, tx_id, chr, tss, end, strand, length, tsl, basic, mane
  exons.bed        all exons with gene biotype + gene_id (col5=score placeholder)
  introns.bed      PCG introns (for intronic lncRNA classification)
Writes per-group gene BEDs + transcript-TSS BEDs for bedtools.
Reference GTF is opened READ-ONLY; all outputs under WORKDIR.
"""
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR

GTF = "/media/share/resource/reference/homo_sapiens/GRCh38/annotation/gencode.v46.annotation.gtf"
OUT = os.path.join(WORKDIR, "data", "genes")

GENE_GROUPS = {
    "pcg": {"protein_coding"},
    "lnc": {"lncRNA"},
    "pseudo": {t for t in (
        "unprocessed_pseudogene processed_pseudogene transcribed_unprocessed_pseudogene "
        "transcribed_processed_pseudogene transcribed_unitary_pseudogene unitary_pseudogene "
        "polymorphic_pseudogene pseudogene rRNA_pseudogene tRNA_pseudogene snRNA_pseudogene "
        "scaRNA_pseudogene miRNA_pseudogene scRNA_pseudogene TR_V_pseudogene TR_J_pseudogene "
        "TR_C_pseudogene TR_D_pseudogene IG_V_pseudogene IG_J_pseudogene IG_C_pseudogene "
        "IG_D_pseudogene").split()},
    "small": {"snRNA snoRNA miRNA misc_RNA rRNA Mt_rRNA Mt_tRNA scaRNA ribozyme sRNA vaultRNA "
              "artifact TEC".split()[i] for i in range(12)},
}
PSEUDO_TRANSCRIBED = {"transcribed_unprocessed_pseudogene", "transcribed_processed_pseudogene",
                      "transcribed_unitary_pseudogene"}


def group_of(gtype):
    for g, members in GENE_GROUPS.items():
        if gtype in members:
            return g
    return "other"


def main():
    with step("02_parse_gencode") as st:
        os.makedirs(OUT, exist_ok=True)
        attr_re = re.compile(r'(\S+) "([^"]*)"')

        genes = {}          # gid -> dict
        tx = {}             # txid -> dict
        exons = []          # (chr, s, e, gid, gtype, strand)

        op = gzip.open if GTF.endswith(".gz") else open
        with op(GTF, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                chrom, feat, start, end, strand = f[0], f[2], int(f[3]), int(f[4]), f[6]
                a = dict(attr_re.findall(f[8]))

                if feat == "gene":
                    gid = a.get("gene_id", "")
                    g = genes.get(gid)
                    if g is None:
                        genes[gid] = dict(gid=gid, gtype=a.get("gene_type", "?"),
                                          name=a.get("gene_name", "?"), chr=chrom,
                                          start=start, end=end, strand=strand,
                                          level=a.get("level", "?"))
                elif feat == "transcript":
                    tid = a.get("transcript_id", "")
                    gid = a.get("gene_id", "")
                    gtype = genes.get(gid, {}).get("gtype", a.get("gene_type", "?"))
                    tags = f[8].split('tag "')
                    ttags = [t.split('"')[0] for t in tags[1:]]
                    tx[tid] = dict(tid=tid, gid=gid, chr=chrom, start=start, end=end,
                                   strand=strand, gtype=gtype,
                                   tsl=a.get("transcript_support_level", "NA"),
                                   basic=1 if "basic" in ttags else 0,
                                   mane=1 if any(t.startswith("MANE_Select") for t in ttags) else 0,
                                   n_exons=0)
                elif feat == "exon":
                    tid = a.get("transcript_id", "")
                    t = tx.get(tid)
                    if t is not None:
                        t["n_exons"] += 1
                        exons.append((chrom, start - 1, end, tid, t["gid"], t["gtype"], strand))

        # gene-level table
        genes_path = safe_write_path(os.path.join(OUT, "genes.tsv"))
        with open(genes_path, "w") as fh:
            fh.write("gene_id\tgene_type\tgroup\tgene_name\tchr\tstart\tend\tstrand\tlevel\ttranscribed_pseudo\n")
            for g in sorted(genes.values(), key=lambda x: (x["chr"], x["start"])):
                grp = group_of(g["gtype"])
                fh.write(f"{g['gid']}\t{g['gtype']}\t{grp}\t{g['name']}\t{g['chr']}\t"
                         f"{g['start']}\t{g['end']}\t{g['strand']}\t{g['level']}\t"
                         f"{1 if g['gtype'] in PSEUDO_TRANSCRIBED else 0}\n")

        # transcript table with TSS (5' end: + strand -> start; - strand -> end)
        txs_path = safe_write_path(os.path.join(OUT, "transcripts.tsv"))
        with open(txs_path, "w") as fh:
            fh.write("tx_id\tgene_id\tgene_type\tgroup\tchr\ttss\ttx_start\ttx_end\tstrand\t"
                     "length\tn_exons\ttsl\tbasic\tmane\n")
            for t in tx.values():
                tss = t["start"] if t["strand"] == "+" else t["end"]
                fh.write(f"{t['tid']}\t{t['gid']}\t{t['gtype']}\t{group_of(t['gtype'])}\t"
                         f"{t['chr']}\t{tss}\t{t['start']}\t{t['end']}\t{t['strand']}\t"
                         f"{t['end']-t['start']+1}\t{t['n_exons']}\t{t['tsl']}\t"
                         f"{t['basic']}\t{t['mane']}\n")

        # exon BED (merged per transcript not needed; keep per-tx for antisense exon logic)
        ex_path = safe_write_path(os.path.join(OUT, "exons.bed"))
        with open(ex_path, "w") as fh:
            for c, s, e, tid, gid, gtype, strand in exons:
                fh.write(f"{c}\t{s}\t{e}\t{gid}|{tid}\t0\t{strand}\n")

        # per-group gene BED + TSS BED
        for grp in ("pcg", "lnc", "pseudo", "small"):
            gset = [g for g in genes.values() if group_of(g["gtype"]) == grp]
            with open(safe_write_path(os.path.join(OUT, f"{grp}.bed")), "w") as fh:
                for g in sorted(gset, key=lambda x: (x["chr"], x["start"])):
                    fh.write(f"{g['chr']}\t{g['start']-1}\t{g['end']}\t{g['gid']}\t0\t{g['strand']}\n")
            with open(safe_write_path(os.path.join(OUT, f"{grp}_tss.bed")), "w") as fh:
                for g in sorted(gset, key=lambda x: (x["chr"], x["start"])):
                    tss = g["start"] if g["strand"] == "+" else g["end"]
                    fh.write(f"{g['chr']}\t{tss-1}\t{tss}\t{g['gid']}\t0\t{g['strand']}\n")

        st.note(f"genes={len(genes)} transcripts={len(tx)} exons={len(exons)}")
        st.note("groups: " + ", ".join(
            f"{g}={sum(1 for x in genes.values() if group_of(x['gtype'])==g)}"
            for g in ("pcg", "lnc", "pseudo", "small", "other")))
        assert len(genes) == 63086, f"expected 63086 genes, got {len(genes)}"

        st.summary_md("results/step_summaries/02_parse_gencode.md", "02 GENCODE parsing",
            f"- Genes: **{len(genes)}** (expected 63,086)\n"
            f"- Transcripts: {len(tx):,}; exon rows: {len(exons):,}\n"
            f"- PCG={sum(1 for x in genes.values() if group_of(x['gtype'])=='pcg')}, "
            f"lncRNA={sum(1 for x in genes.values() if group_of(x['gtype'])=='lnc')}, "
            f"pseudogene={sum(1 for x in genes.values() if group_of(x['gtype'])=='pseudo')}, "
            f"small/other={sum(1 for x in genes.values() if group_of(x['gtype']) in ('small','other'))}\n")


if __name__ == "__main__":
    main()
