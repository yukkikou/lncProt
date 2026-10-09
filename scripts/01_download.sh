#!/bin/bash
# =============================================================================
# 01_download.sh — fetch all external data into $DATA (resumable, verified).
# Reference tree untouched. All URLs verified reachable 2026-09-30.
# =============================================================================
source "$(dirname "$0")/common.sh"
trap 'on_err 01_download $BASH_COMMAND' ERR
CURRENT_STEP="01_download"
progress_init

FANTOM_BASE="https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra"

fetch() { # fetch <url> <dest>  — persistent retry; permanent 4xx fails fast
    local url="$1" dest="$2"
    guard_out "$dest"
    mkdir -p "$(dirname "$dest")"
    if [ -s "$dest" ] && [ ! "${3:-}" = "--force" ]; then
        echo ">>> exists, skip: $dest"
    else
        echo ">>> downloading $(basename "$dest")"
        # permanent-failure probe (some servers lack HEAD; use 1-byte range GET)
        local code
        code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 30 --range 0-0 "$url" || echo 000)
        if [ "$code" = "404" ] || [ "$code" = "403" ] || [ "$code" = "410" ]; then
            echo "PERMANENT FAIL ($code): $url" >&2
            return 1
        fi
        local attempt=0 max=40 ok=0
        while [ $attempt -lt $max ]; do
            attempt=$((attempt+1))
            # --speed-limit 1000 --speed-time 120: abort if <1KB/s for 2min (stall kill)
            if curl -fsSL --retry 2 --retry-delay 15 -C - --connect-timeout 40 \
                 --speed-limit 1000 --speed-time 120 -o "$dest.part" "$url"; then
                ok=1; break
            fi
            local rc=$?
            if [ "$rc" = "22" ]; then   # HTTP error again -> re-probe, maybe permanent
                code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 30 --range 0-0 "$url" || echo 000)
                [ "$code" = "404" ] || [ "$code" = "403" ] || [ "$code" = "410" ] && \
                    { echo "PERMANENT FAIL ($code): $url" >&2; return 1; }
            fi
            echo "    attempt $attempt failed (rc=$rc), retrying in 30s ($(date +%H:%M:%S))"
            sleep 30
        done
        [ $ok -eq 1 ] || { echo "DOWNLOAD FAILED after $max attempts: $url" >&2; return 1; }
        mv "$dest.part" "$dest"
    fi
    local sz; sz=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    [ "$sz" -gt 1000 ] || { echo "DOWNLOAD TOO SMALL: $dest ($sz bytes)" >&2; return 1; }
    echo "    ok: $sz bytes"
}

echo ">>> $(date) 01_download start"

# ---- 1. UCSC CpG islands (hg38) ---------------------------------------------
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz" \
      "$DATA/cpgIslandExt.txt.gz"
# convert to simple bed: chr, start, end, name=isl name (+length info kept in name)
zcat "$DATA/cpgIslandExt.txt.gz" | awk -F'\t' 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}' \
    > "$DATA/cpgIslandExt.full.txt" 2>/dev/null || true

# ---- 2. FANTOM5 hg38 native: CAGE peaks + enhancers (critical, small) --------
# NOTE: annotation table only exists for the "liftover+new" flavor, not
# "fair+new"; peaks (here) and TPM matrix (below) share the fair+new peak IDs,
# so no annot join is needed — locus classes are computed from coordinates.
fetch "$FANTOM_BASE/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz" \
      "$DATA/fantom/hg38_cage_peaks.bed.gz"
fetch "$FANTOM_BASE/enhancer/F5.hg38.enhancers.bed.gz" \
      "$DATA/fantom/F5.hg38.enhancers.bed.gz"

# ---- 3. Roadmap ChromHMM + motifs (different hosts, usually fast) -------------
for E in $CHROMHMM_CELLS; do
    fetch "$CHROMHMM_URLBASE/${E}_15_coreMarks_dense.bed.gz" \
          "$DATA/chromhmm_hg19/${E}_15_coreMarks_dense.bed.gz" || \
      echo "WARNING: ChromHMM $E failed (module will degrade gracefully)"
done
fetch "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.zip" \
      "$DATA/motifs/JASPAR2024_vertebrates_meme.zip"
( cd "$DATA/motifs" && unzip -o JASPAR2024_vertebrates_meme.zip )
fetch "https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE/H14CORE_pwm.tar.gz" \
      "$DATA/motifs/H14CORE_pwm.tar.gz" && \
( cd "$DATA/motifs" && tar xzf H14CORE_pwm.tar.gz ) || \
  { echo "WARNING: HOCOMOCO v14 unpack failed — fallback to local v11 ($HOCOMOCO_LOCAL)"; }

# ---- 4. FANTOM TPM matrix LAST (844MB, slow mirror; non-critical for 04) -------
fetch "$FANTOM_BASE/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz" \
      "$DATA/fantom/hg38_cage_tpm_ann.osc.txt.gz" || \
  { echo "WARNING: TPM matrix failed; trying rel_expr fallback"; \
    fetch "$FANTOM_BASE/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_rel_expr.txt.gz" \
          "$DATA/fantom/hg38_cage_rel_expr.txt.gz" || \
    echo "WARNING: no expression matrix — 09 will run peaks-only (binary activity) mode"; }

# list what we got
echo ">>> inventory:"
find "$DATA" -type f \( -name "*.gz" -o -name "*.meme" -o -name "*.txt" -o -name "*.PWM" \) -printf "%s\t%p\n" | sort -k2

progress "01_download" "DONE" "files=$(find "$DATA" -type f | wc -l)"
echo ">>> $(date) 01_download complete"
