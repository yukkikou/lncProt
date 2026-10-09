#!/usr/bin/env python
"""10_consv_chromatin.py — conservation (phastCons470) + ChromHMM states.

phastCons: computed by scripts/pybw_stats.py run INSIDE pyBigWig.sif.
ChromHMM: dominant 15-state (largest overlap) within ±1kb of TSS per cell type.
Outputs: results/conservation.tsv, results/chromhmm_states.tsv
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, WORKDIR, df_to_md

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data")
SIF_BT = "/media/share/resource/container/bedtools/bedtools-2.30.0.sif"
SIF_PBW = "/media/share/resource/container/pyBigWig/pyBigWig.sif"
CHROMHMM_DIR = os.path.join(DATA, "chromhmm_hg38")

# Roadmap 15-state mnemonics (as written by 03b: col4 = mnemonic only)
PROM_LIKE = {"TssA", "TssAFlnk"}
ENH_LIKE = {"EnhG", "Enh"}
BIV_LIKE = {"TssBiv", "BivFlnk", "EnhBiv"}
REPR_LIKE = {"ReprPC", "ReprPCWk"}


def sing(*args):
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}"] + list(args)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(" ".join(cmd) + "\n" + r.stderr[:1000])
    return r.stdout


def main():
    with step("10_consv_chromatin") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                           index_col=0)

        # ---- phastCons via pyBigWig container
        raw = os.path.join(RES, "conservation_raw.tsv")
        if not os.path.exists(raw):
            sing(SIF_PBW, "python",
                 os.path.join(WORKDIR, "scripts", "pybw_stats.py"))
        cons = pd.read_csv(raw, sep="\t", index_col=0)
        out = meta.join(cons, how="left")
        out.to_csv(safe_write_path(os.path.join(RES, "conservation.tsv")), sep="\t")
        st.note(f"phastCons mean (PCG)={out[out.group=='pcg']['phast_mean'].mean():.3f}")

        # ---- ChromHMM dominant state per promoter per cell type
        prom = pd.read_csv(os.path.join(RES, "promoters_feat.bed"), sep="\t",
                           names=["chrom", "start", "end", "gene_id", "x", "strand"])
        prom["w_start"] = (prom["start"] + prom["end"]) // 2 - 1000
        prom["w_end"] = (prom["start"] + prom["end"]) // 2 + 1000
        win = prom[["chrom", "w_start", "w_end", "gene_id"]].copy()
        win.columns = ["chrom", "start", "end", "gene_id"]
        win = win.sort_values(["chrom", "start"])
        win_bed = safe_write_path(os.path.join(DATA, "prom_tss_2kb.bed"))
        win.to_csv(win_bed, sep="\t", index=False, header=False)

        cells = sorted(f[:-4] for f in os.listdir(CHROMHMM_DIR)
                       if f.endswith(".bed")) if os.path.isdir(CHROMHMM_DIR) else []
        st.note(f"chromhmm cells available={cells}")
        state_tbl = pd.DataFrame(index=meta.index)
        for E in cells:
            bed = os.path.join(CHROMHMM_DIR, E + ".bed")
            ovf = os.path.join(DATA, f"chromhmm_ov_{E}.tsv")
            cmd = ["singularity", "exec", "--no-mount", "bind-paths",
                   "--bind", "/media/share/resource:/media/share/resource:ro",
                   "--bind", f"{WORKDIR}:{WORKDIR}",
                   SIF_BT, "bedtools", "intersect", "-a", win_bed, "-b", bed, "-wo"]
            with open(safe_write_path(ovf), "w") as fh:
                subprocess.run(cmd, stdout=fh)
            best = {}
            with open(ovf) as fh:
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if len(f) < 9 or f[4] == ".":
                        continue
                    gid, state, ov = f[3], f[7], int(f[8])
                    if gid not in best or ov > best[gid][1]:
                        best[gid] = (state, ov)
            state_tbl[E] = pd.Series({g: s for g, (s, _) in best.items()})
            st.note(f"{E}: assigned={len(best)} "
                    f"TssA={sum(1 for s,_ in best.values() if s=='TssA')}")
        state_tbl.to_csv(safe_write_path(os.path.join(RES, "chromhmm_states.tsv")),
                         sep="\t")

        # per-class state distribution summary (state 1/2/3/4 vs enh 7/8/9 vs 12/13)
        smrows = []
        for g in ["PCG", "divergent", "antisense", "intronic", "sense_overlapping",
                  "intergenic", "pseudogene", "transcribed_pseudogene"]:
            sub_idx = meta[(meta.group == "pcg")].index if g == "PCG" \
                else meta[meta.cls == g].index
            if not len(sub_idx) or not cells:
                continue
            sub = state_tbl.reindex(sub_idx)
            row = dict(cls=g, n=len(sub_idx))
            for E in cells:
                v = sub[E].dropna().astype(str)
                row[f"{E}_prom_like"] = v.isin(PROM_LIKE).mean()
                row[f"{E}_enh_like"] = v.isin(ENH_LIKE).mean()
                row[f"{E}_biv_like"] = v.isin(BIV_LIKE).mean()
                row[f"{E}_repressed"] = v.isin(REPR_LIKE).mean()
            smrows.append(row)
        sm = pd.DataFrame(smrows)
        sm.to_csv(safe_write_path(os.path.join(RES, "chromhmm_by_class.tsv")),
                  sep="\t", index=False)
        c0 = f"{cells[0]}_prom_like" if cells else None
        if c0 and c0 in sm.columns:
            st.note("chromhmm by class: " + "; ".join(
                f"{r.cls}: prom={getattr(r, c0):.2f}" for r in sm.itertuples()))
        st.summary_md("results/step_summaries/10_consv_chromatin.md",
                      "10 conservation + chromatin",
                      df_to_md(sm, ".2f") if len(sm) else "_(no cells)")


if __name__ == "__main__":
    main()
