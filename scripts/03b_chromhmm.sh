#!/bin/bash
# =============================================================================
# 03b_chromhmm.sh — lift Roadmap ChromHMM hg19 segments to hg38 + two-level QC:
#   A) chainlift validation: hg19 CpG islands -> overlap with native hg38 CpG
#      islands (must be high, validates our liftover implementation)
#   B) module QC: PCG TSS enrichment in TssActive/TssFlnk states per cell type
# =============================================================================
source "$(dirname "$0")/common.sh"
trap 'on_err 03b_chromhmm $BASH_COMMAND' ERR
CURRENT_STEP="03b_chromhmm"
progress_init

echo ">>> $(date) 03b start"

# ---- 0) chain file: Ensembl assembly-mapping chain (reference copy truncated) --
CHAIN="$DATA/qc/grch37_to_grch38.chain"
if [ ! -s "$CHAIN" ]; then
    curl -fsSL --retry 3 -o "$DATA/qc/grch37_to_grch38.chain.gz" \
        "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz"
    gunzip -kf "$DATA/qc/grch37_to_grch38.chain.gz"
fi

# ---- A) chainlift validation via CpG islands ---------------------------------
mkdir -p "$DATA/qc"
curl -fsSL --retry 3 -o "$DATA/qc/hg19_cpgIslandExt.txt.gz" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz"
zcat "$DATA/qc/hg19_cpgIslandExt.txt.gz" | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' > "$DATA/qc/hg19_cpg.bed"
zcat "$DATA/cpgIslandExt.txt.gz"         | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' > "$DATA/qc/hg38_cpg_native.bed"

"$PY" scripts/03a_chain_liftover.py --in "$DATA/qc/hg19_cpg.bed" \
    --out "$DATA/qc/hg19_cpg.hg38.bed" --label 03a_cpg_validation

sort -k1,1 -k2,2n "$DATA/qc/hg19_cpg.hg38.bed" > "$DATA/qc/hg19_cpg.hg38.sorted.bed"
sort -k1,1 -k2,2n "$DATA/qc/hg38_cpg_native.bed" > "$DATA/qc/hg38_cpg_native.sorted.bed"
N_LIFT=$(wc -l < "$DATA/qc/hg19_cpg.hg38.sorted.bed")
N_HIT=$(bt intersect -a "$DATA/qc/hg19_cpg.hg38.sorted.bed" -b "$DATA/qc/hg38_cpg_native.sorted.bed" -sorted -u | wc -l)
RATE=$(awk -v a="$N_HIT" -v b="$N_LIFT" 'BEGIN{if(b>0) printf "%.3f", a/b; else print 0}')
echo "chainlift CpG validation: $N_HIT / $N_LIFT lifted hg19 CpG islands overlap a native hg38 island (rate=$RATE)"
echo "chainlift_cpg_overlap_rate=$RATE" > "$DATA/qc/chainlift_qc.txt"
# fail hard if liftover is broken (<70% concordance)
awk -v r="$RATE" 'BEGIN{exit !(r>=0.70)}' || { echo "FATAL: chainlift validation failed (<0.70) — ChromHMM module degraded"; progress "03b_chromhmm" "FAILED" "chainlift_qc_rate=$RATE"; exit 1; }

# ---- B) lift each Roadmap cell type ------------------------------------------
for E in $CHROMHMM_CELLS; do
    SRC="$DATA/chromhmm_hg19/${E}_15_coreMarks_dense.bed.gz"
    [ -s "$SRC" ] || { echo "WARN: missing $SRC — skip $E"; continue; }
    zcat "$SRC" | awk 'BEGIN{OFS="\t"}{split($4,a,"_"); print $1,$2,$3,a[2]}' > "$DATA/qc/${E}.hg19.bed"
    "$PY" scripts/03a_chain_liftover.py --in "$DATA/qc/${E}.hg19.bed" \
        --out "$DATA/chromhmm_hg38/${E}.bed" --label "03a_lift_$E"
    sort -k1,1 -k2,2n "$DATA/chromhmm_hg38/${E}.bed" > "$DATA/chromhmm_hg38/${E}.sorted.bed"
    mv "$DATA/chromhmm_hg38/${E}.sorted.bed" "$DATA/chromhmm_hg38/${E}.bed"
done

# ---- module QC: PCG TSS in state 1/2 ------------------------------------------
"$PY" - <<'PYEOF'
import os, sys
sys.path.insert(0, "scripts")
from pylib.progress import step
DATA = "data"
with step("03b_chromhmm_qc") as st:
    # 2kb windows around PCG major-TSS (approx: gene 5' ends from pcg.bed gene bodies)
    import subprocess
    lines = [l.split() for l in open(f"{DATA}/genes/pcg.bed")]
    with open(f"{DATA}/qc/pcg_tss2k.bed", "w") as fh:
        for c, s, e, gid, _, strand in lines:
            tss = int(s) + 1 if strand == "+" else int(e)
            fh.write(f"{c}\t{max(0,tss-2000)}\t{tss+2000}\t{gid}\n")
    out = []
    for E in os.environ.get("CHROMHMM_CELLS", "").split():
        p = f"{DATA}/chromhmm_hg38/{E}.bed"
        if not os.path.exists(p):
            continue
        n_hit = n_tot = 0
        # dominant state per TSS = state of largest overlap
        wd = os.path.abspath(".")
        r = subprocess.run(["singularity","exec","--no-mount","bind-paths",
            "--bind","/media/share/resource:/media/share/resource:ro",
            "--bind",f"{wd}:{wd}",
            os.environ["SIF_BEDTOOLS"],"bedtools","intersect",
            "-a",f"{DATA}/qc/pcg_tss2k.bed","-b",p,"-wao"],
            capture_output=True, text=True)
        best = {}
        for ln in r.stdout.splitlines():
            f = ln.split("\t")
            if len(f) < 9 or f[4] == ".":
                continue
            gid, state, ov = f[3], f[7], int(f[8])
            if gid not in best or ov > best[gid][1]:
                best[gid] = (state, ov)
        n_tot = len(lines)
        n_hit = sum(1 for s,_ in best.values() if s in ("TssA","TssAFlnk"))
        out.append((E, n_hit, n_tot, n_hit/max(n_tot,1)))
    body = "| cell | TSS in state 1/2 | n PCG | rate |\n|---|---|---|---|\n" + \
        "\n".join(f"| {e} | {h} | {t} | {r:.2f} |" for e,h,t,r in out)
    for e,h,t,r in out:
        st.note(f"{e}: TssActive/TssFlnk at PCG TSS={r:.2f}")
    st.summary_md("results/step_summaries/03b_chromhmm.md", "03b ChromHMM liftover QC", body)
PYEOF

progress "03b_chromhmm" "DONE" "cells=$(ls $DATA/chromhmm_hg38/*.bed 2>/dev/null | wc -l)"
[ -n "$(ls -A $DATA/chromhmm_hg38/*.bed 2>/dev/null)" ] && touch "$DATA/chromhmm_hg38/.done"
echo ">>> $(date) 03b complete"
