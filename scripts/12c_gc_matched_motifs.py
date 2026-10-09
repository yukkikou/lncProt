#!/usr/bin/env python
"""12c_gc_matched_motifs.py — are candidate-TF motif enrichments real, or GC
artifacts? Absolute-enrichment control: rescan TOP candidate TF motifs on
dinucleotide-matched synthetic promoters (1st-order Markov sequences generated
from each promoter's own dinucleotide composition), compare hit rates.

Outputs: results/gc_matched_motif_control.tsv + summary
"""
import os
import random
import subprocess
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data")
SIF_MEME = "/media/share/resource/container/meme/meme-5.4.1.sif"
HOCO = ("/media/share/resource/database/MotifDB/HUMAN/"
        "HOCOMOCOv11_core_HUMAN_mono_meme_format.meme")


def markov_sequence(dinuc_counts, length, rng):
    """Generate a sequence from 1st-order Markov chain matching dinucleotide
    composition (mono + di frequencies preserved in expectation)."""
    # transition counts: base -> Counter(next base)
    trans = defaultdict(Counter)
    mono = Counter()
    for (a, b), n in dinuc_counts.items():
        trans[a][b] += n
        mono[a] += n
    first = list(mono.keys())
    w = [mono[f] for f in first]
    cur = rng.choices(first, weights=w)[0]
    out = [cur]
    for _ in range(length - 1):
        nxt = list(trans[cur].keys())
        if not nxt:
            cur = rng.choices(first, weights=w)[0]
        else:
            cur = rng.choices(nxt, weights=[trans[cur][x] for x in nxt])[0]
        out.append(cur)
    return "".join(out)


def parse_meme(path):
    """-> (header_lines, [motif_block_lines])"""
    header, blocks, cur = [], [], None
    for ln in open(path):
        if ln.startswith("MOTIF "):
            cur = [ln]
            blocks.append(cur)
        elif cur is None:
            header.append(ln)
        else:
            cur.append(ln)
    return header, blocks


def main():
    with step("12c_gc_matched_motifs") as st:
        rng = random.Random(42)
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                           index_col=0, low_memory=False)
        # top candidate TFs per class (motif-enriched, FDR<0.05)
        me = pd.read_csv(os.path.join(RES, "motif_enrichment.tsv"), sep="\t")
        top = (me[(me.fdr < 0.05) & (me.log2OR > 0.3)]
               .sort_values("log2OR", ascending=False)
               .groupby("cls").head(15))
        classes = [c for c in ("divergent", "intergenic") if c in set(top.cls)]
        tfs = sorted(set(top[top.cls.isin(classes)].tf))
        st.note(f"classes={classes} n_top_TFs={len(tfs)}")

        # motif blocks for these TFs from HOCOMOCO (symbols match motif names)
        header, blocks = parse_meme(HOCO)
        want = {f"{t}.HUMAN" for t in tfs}
        sel, got = [], set()
        tf_up = {t.upper(): t for t in tfs}
        for b in blocks:
            name = b[0].strip().split()[1].upper()
            for tu, t in tf_up.items():
                if name.startswith(f"{tu}_HUMAN") or name.startswith(f"{tu}.HUMAN"):
                    sel.append(b)
                    got.add(t)
                    break
        st.note(f"motifs found for {len(got)}/{len(tfs)} TFs in HOCOMOCO")
        motif_path = safe_write_path(os.path.join(DATA, "qc", "top_tfs.meme"))
        with open(motif_path, "w") as fh:
            fh.write("MEME version 4\nALPHABET= ACGT\n")
            for b in sel:
                fh.writelines(b)

        # promoter sequences per class + synthetic Markov-matched set
        seqs = {}
        name = None
        buf = []
        for line in open(os.path.join(RES, "promoters_fimo.fa")):
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf)

        results = []
        for cls in classes:
            genes = meta[(meta.cls == cls)].index
            real = [seqs[g] for g in genes if g in seqs and len(seqs[g]) > 100]
            syn = []
            for s in real:
                di = Counter(zip(s[:-1], s[1:]))
                syn.append(markov_sequence(di, len(s), rng))
            for tag, dataset in (("observed", real), ("gc_matched", syn)):
                fa = safe_write_path(os.path.join(DATA, "qc", f"mc_{cls}_{tag}.fa"))
                with open(fa, "w") as fh:
                    for i, s in enumerate(dataset):
                        fh.write(f">s{i}\n{s}\n")
                oc = os.path.join(DATA, "qc", f"fimo_mc_{cls}_{tag}")
                cmd = ["singularity", "exec", "--no-mount", "bind-paths",
                       "--bind", "/media/share/resource:/media/share/resource:ro",
                       "--bind", f"{WORKDIR}:{WORKDIR}", SIF_MEME, "fimo",
                       "--no-qvalue", "--thresh", "1e-4", "--oc", oc,
                       motif_path, fa]
                r = subprocess.run(cmd, capture_output=True, text=True)
                if r.returncode != 0:
                    raise RuntimeError(r.stderr[:800])
                hits = set()
                p = os.path.join(oc, "fimo.tsv")
                with open(p) as fh:
                    fh.readline()
                    for ln in fh:
                        f = ln.rstrip("\n").split("\t")
                        if len(f) > 3 and f[2] and not f[2].startswith("#"):
                            hits.add(f[2])          # sequence_name column
                results.append(dict(cls=cls, dataset=tag, n_promoters=len(dataset),
                                    frac_with_any_hit=len(hits) / max(len(dataset), 1)))
                st.note(f"{cls} {tag}: {len(hits)}/{len(dataset)} promoters with hits")
        t = pd.DataFrame(results)
        piv = t.pivot(index="cls", columns="dataset", values="frac_with_any_hit")
        piv["observed/gc_matched_ratio"] = piv["observed"] / piv["gc_matched"]
        piv.to_csv(safe_write_path(
            os.path.join(RES, "gc_matched_motif_control.tsv")), sep="\t")
        st.note("ratios: " + "; ".join(
            f"{r.cls}: {r._2:.2f}" if False else f"{i}: {v:.2f}"
            for i, v in piv["observed/gc_matched_ratio"].items()))
        st.summary_md("results/step_summaries/12c_gc_matched_motifs.md",
                      "12c GC-matched motif control",
                      "Fraction of promoters with >=1 hit from the top-TF set:\n\n"
                      + df_to_md(t, ".3f")
                      + "\n\nA ratio >>1 means enrichment is NOT explained by "
                        "dinucleotide composition.")


if __name__ == "__main__":
    main()
