#!/bin/bash
# =============================================================================
# common.sh — shared shell helpers: singularity wrappers, progress ledger,
# failure recording, output guards. Sourced by all pipeline shell scripts.
# =============================================================================
# Conventions:
#   * Reference tree (/media/share/resource) mounted READ-ONLY into containers
#   * All outputs must live under $WORKDIR (guard_out enforces)
#   * Progress -> $LOGS/PROGRESS.md ; failures -> $LOGS/failures/<step>.err
# =============================================================================

source "${BASH_SOURCE[0]%/*}/../config.sh"

# ---- Singularity wrappers (system config binds nonexistent /work -> disable) --
SING="singularity exec --no-mount bind-paths"
# Mount resource tree read-only, workdir read-write:
binds=(--bind /media/share/resource:/media/share/resource:ro --bind "$WORKDIR:$WORKDIR")
bt()    { $SING "${binds[@]}" "$SIF_BEDTOOLS"  bedtools "$@"; }
fimo()  { $SING "${binds[@]}" "$SIF_MEME"      fimo     "$@"; }
pbwpy() { $SING "${binds[@]}" "$SIF_PYBIGWIG"  python   "$@"; }

# ---- Output path guard: refuse any write target outside $WORKDIR -------------
guard_out() {
    local p
    for p in "$@"; do
        case "$(readlink -m "$p" 2>/dev/null || echo "$p")/" in
            "$WORKDIR"/*) ;;
            *) echo "REFUSED: write target outside WORKDIR: $p" >&2; return 1 ;;
        esac
    done
}

# ---- Progress ledger -----------------------------------------------------------
# progress <step> <status> <key_numbers>
progress() {
    local step="$1" status="$2" nums="${3:-}"
    local ts; ts=$(date '+%Y-%m-%d %H:%M:%S')
    guard_out "$LOGS/PROGRESS.md" || return 0
    printf '| %s | %s | %s | %s |\n' "$step" "$status" "$ts" "$nums" >> "$LOGS/PROGRESS.md"
}
progress_init() {
    [ -f "$LOGS/PROGRESS.md" ] || { guard_out "$LOGS/PROGRESS.md" && \
        printf '| step | status | time | key numbers |\n|---|---|---|---|\n' > "$LOGS/PROGRESS.md"; }
}

# ---- Failure trap ----------------------------------------------------------------
# fail_trap <step> <line> <command>  (installed via: trap 'on_err STEPNAME $LINENO $BASH_COMMAND' ERR)
CURRENT_STEP="shell"
on_err() {
    local step="$CURRENT_STEP" line="$1" cmd="$2" code=$?
    progress "$step" "FAILED(code=$code)" "line=$line cmd=${cmd:0:120}"
    mkdir -p "$LOGS/failures"
    { echo "step=$step line=$line exit=$code"
      echo "cmd: $cmd"
      echo "--- last stderr (from SLURM log if available) ---"
      local jid="${SLURM_JOB_ID:-none}"
      [ "$jid" != none ] && tail -n 50 "/media/share/home/1610305236/log/"*"$jid"*.err 2>/dev/null
    } > "$LOGS/failures/${step}.err" 2>/dev/null || true
    guard_out "$RESULTS/.${step}.failed" && touch "$RESULTS/.${step}.failed" || true
    return $code
}

# ---- Reference integrity manifest -------------------------------------------------
ref_manifest() {  # ref_manifest build|check
    local mode="$1" mf="$LOGS/00_ref_manifest.txt"
    guard_out "$mf" || return 1
    local files=( "$GENCODE" "$GENOME" "$TE_BED" "$HERV_BED" \
                  "$PHASTCONS_BW" "$REMAP_NR" "$REMAP_ALL" "$HOCOMOCO_LOCAL" "$JASPAR_LOCAL" )
    case $mode in
      build)
        { echo "# reference manifest (size:mtime:sha256-if<100MB)  built $(date)"
          for f in "${files[@]}"; do
            [ -f "$f" ] || { echo "MISSING $f"; continue; }
            local sz mt h
            sz=$(stat -c%s "$f"); mt=$(stat -c%y "$f" | cut -d. -f1)
            if [ "$sz" -lt 100000000 ]; then h=$(sha256sum "$f" | cut -d' ' -f1); else h="(large:skip-hash)"; fi
            echo -e "$sz\t$mt\t$h\t$f"
          done; } > "$mf" ;;
      check)
        local tmp; tmp=$(mktemp "$WORKDIR/tmp_manifest.XXXXXX")
        local files_bak; files_bak=()
        cp "$mf" "$tmp"
        # rebuild to second file and compare
        local mf2="$LOGS/.00_ref_manifest.check.txt"
        { echo "# reference manifest (size:mtime:sha256-if<100MB)  built $(date)"
          for f in "${files[@]}"; do
            [ -f "$f" ] || { echo "MISSING $f"; continue; }
            local sz mt h
            sz=$(stat -c%s "$f"); mt=$(stat -c%y "$f" | cut -d. -f1)
            if [ "$sz" -lt 100000000 ]; then h=$(sha256sum "$f" | cut -d' ' -f1); else h="(large:skip-hash)"; fi
            echo -e "$sz\t$mt\t$h\t$f"
          done; } > "$mf2"
        if diff <(tail -n +2 "$tmp") <(tail -n +2 "$mf2") > /dev/null; then
            echo "REFERENCE INTEGRITY: OK (zero changes)"; rm -f "$tmp" "$mf2"; return 0
        else
            echo "REFERENCE INTEGRITY: VIOLATION — changes detected:" >&2
            diff <(tail -n +2 "$tmp") <(tail -n +2 "$mf2") >&2 || true
            rm -f "$tmp"; return 1
        fi ;;
    esac
}
