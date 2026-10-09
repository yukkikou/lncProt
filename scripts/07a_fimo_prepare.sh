#!/bin/bash
# =============================================================================
# 07a_fimo_scan.sh — FIMO motif scan over union promoter FASTA (-500/+100).
# Two-phase: 'prepare' merges motif DBs and splits into chunks;
#            'scan' (run per array task via SLURM_ARRAY_TASK_ID) scans one chunk.
# Usage: 07a_fimo_scan.sh prepare | scan
# =============================================================================
source "$(dirname "$0")/common.sh"
CURRENT_STEP="07a_fimo"
progress_init

MODE="${1:-prepare}"
MOTIF_DIR="$DATA/motifs"
CHUNK_DIR="$MOTIF_DIR/chunks"
N_CHUNKS=$MAXCHUNK_NJOBS

merge_motifs() {
    mkdir -p "$CHUNK_DIR"
    # HOCOMOCO: v14 tarball ships raw log-odds .pwm (not MEME) — use the
    # verified local v11 core MEME bundle (401 human motifs) instead.
    echo "HOCOMOCO layer: local v11 core"
    cp "$HOCOMOCO_LOCAL" "$CHUNK_DIR/hoco.meme"
    if [ -s "$MOTIF_DIR/JASPAR2024_vertebrates_meme.zip" ]; then
        ( cd "$MOTIF_DIR" && unzip -oq JASPAR2024_vertebrates_meme.zip )
    fi
    # JASPAR zip unpacks to many single-motif MA*.meme files — parse MOTIF
    # blocks properly and concat under the first file's header.
    "$PY" - "$MOTIF_DIR" "$CHUNK_DIR" <<'PYEOF'
import glob, os, sys
mdir, cdir = sys.argv[1], sys.argv[2]

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

jas = sorted(glob.glob(os.path.join(mdir, "MA*.meme")))
header, allblocks = None, []
for f in jas:
    h, b = parse_meme(f)
    if header is None and h:
        header = h
    allblocks.extend(b)
with open(os.path.join(cdir, "jaspar.meme"), "w") as fh:
    fh.writelines(header or ["MEME version 4\n", "ALPHABET= ACGT\n"])
    for blk in allblocks:
        fh.writelines(blk)
print(f"jaspar motifs: {len(allblocks)}")
PYEOF
    # combined library via the same strict parser
    "$PY" - "$CHUNK_DIR" <<'PYEOF'
import os, sys
cdir = sys.argv[1]

def parse_meme(path):
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

h1, b1 = parse_meme(os.path.join(cdir, "hoco.meme"))
h2, b2 = parse_meme(os.path.join(cdir, "jaspar.meme"))
with open(os.path.join(cdir, "all.meme"), "w") as fh:
    fh.writelines(h1)
    for blk in b1 + b2:
        fh.writelines(blk)
print(f"merged motifs: {len(b1) + len(b2)}")
PYEOF
}

split_chunks() {
    rm -f "$CHUNK_DIR"/chunk_*.meme
    "$PY" - "$CHUNK_DIR/all.meme" "$N_CHUNKS" <<'PYEOF'
import os, sys
path, n = sys.argv[1], int(sys.argv[2])
header, blocks, cur = [], [], None
for ln in open(path):
    if ln.startswith("MOTIF "):
        cur = [ln]
        blocks.append(cur)
    elif cur is None:
        header.append(ln)
    else:
        cur.append(ln)
bufs = [[] for _ in range(n)]
for i, blk in enumerate(blocks):
    bufs[i % n].append(blk)
for k, chunks in enumerate(bufs):
    with open(os.path.join(os.path.dirname(path), f"chunk_{k}.meme"), "w") as fh:
        fh.writelines(header)
        for c in chunks:
            fh.writelines(c)
    print(f"chunk_{k}: {len(chunks)} motifs")
PYEOF
}

scan_chunk() {
    local k="$SLURM_ARRAY_TASK_ID"
    [ -z "$k" ] && k=0
    local out="$DATA/fimo/chunk_$k"
    [ -d "$out" ] && [ -s "$out/fimo.tsv" ] && { echo ">>> fimo chunk $k exists, skip"; return 0; }
    mkdir -p "$out"
    echo ">>> fimo scan chunk $k $(date)"
    singularity exec --no-mount bind-paths \
        --bind /media/share/resource:/media/share/resource:ro --bind "$WORKDIR:$WORKDIR" \
        "$SIF_MEME" fimo --no-qvalue --thresh "$FIMO_P" \
        --oc "$out" "$CHUNK_DIR/chunk_$k.meme" "$RESULTS/promoters_fimo.fa" \
        || { echo "FIMO chunk $k failed" >&2; return 1; }
    # keep only the essential columns, gzipped
    awk -F'\t' 'NR==1 || $1!=""{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' \
        "$out/fimo.tsv" 2>/dev/null | gzip -c > "$out/fimo.hits.tsv.gz" || true
    echo "    chunk $k done: $(zcat "$out/fimo.hits.tsv.gz" | wc -l) rows"
}

case $MODE in
  prepare) merge_motifs; split_chunks; progress "07a_fimo_prepare" "DONE" ;;
  scan)    trap 'on_err 07a_fimo_scan $BASH_COMMAND' ERR
           scan_chunk; progress "07a_fimo_scan_$SLURM_ARRAY_TASK_ID" "DONE" ;;
esac
