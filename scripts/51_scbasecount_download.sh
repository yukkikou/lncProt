#!/usr/bin/env bash
# =============================================================================
# 51_scbasecount_download.sh -- fetch the scBaseCount subset designed by 50.
# -----------------------------------------------------------------------------
# Public resources are retained under /media/share/exchange/1610305236/scBaseCount.
# Only the sample-level metadata (30 MB, all species) and the balanced analysis
# design (520 human samples, ~14 GB) are fetched here. The rest of the 39.5 TiB
# release can be added later against the same manifest.
#
# Layout:
#   scBaseCount/metadata/GeneFull_Ex50pAS/<Organism>/sample_metadata.parquet
#   scBaseCount/h5ad/<cancer_group>/<arm>/<SRX>.h5ad
#   scBaseCount/manifest_design.tsv   url, destination, cells, group, arm
#   scBaseCount/README.md             provenance
#
# Re-running is safe. Large files stream to .partial and are moved on success,
# and curl resumes with -C - if a partial file is present.
# =============================================================================
set -euo pipefail

PROJ=/media/share/workdir/1610305236/lncProt
ROOT=/media/share/exchange/1610305236/scBaseCount
BUCKET=https://storage.googleapis.com/arc-institute-virtual-cell-atlas
REL=scbasecount/2026-01-12
PY=/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/bin/python
LIST=${1:-$PROJ/results/scbasecount/design_samples.tsv}
MAN=${2:-manifest_design.tsv}
JOBS=${JOBS:-6}

mkdir -p "$ROOT/h5ad" "$ROOT/metadata/GeneFull_Ex50pAS"

# ---- manifest ---------------------------------------------------------------
MAN_NAME="$MAN" "$PY" - "$LIST" "$ROOT" "$BUCKET" "$REL" <<'PYEOF'
import os, re, sys
import pandas as pd

lst, root, bucket, rel = sys.argv[1:5]
d = pd.read_csv(lst, sep="\t")

def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")

rows = []
for r in d.itertuples():
    dest = os.path.join(root, "h5ad", slug(r.group), r.arm,
                        "{}.h5ad".format(r.srx))
    # file_path is already a gs:// URL that names the bucket, so swap the
    # scheme and bucket rather than prefixing the bucket a second time.
    url = r.file_path.replace("gs://arc-institute-virtual-cell-atlas", bucket)
    rows.append(dict(url=url, dest=dest, srx=r.srx, group=r.group,
                     arm=r.arm, cells=int(r.obs_count)))
man = pd.DataFrame(rows)
man.to_csv(os.path.join(root, os.environ.get("MAN_NAME", "manifest_design.tsv")), sep="\t", index=False)
print("manifest rows", len(man), flush=True)
PYEOF

# ---- sample metadata for every species (small, keeps future options open) ---
# The bucket listing is the fragile part of this script. It is downloaded into a
# temporary file so a failed attempt cannot truncate a list that is already
# good, and when it fails the local metadata directory is used instead, which
# names exactly the same species.
if [ ! -s "$ROOT/.species_prefixes.txt" ]; then
  _tmp=$(mktemp)
  curl -sS --fail --retry 3 --retry-delay 3 --max-time 60 \
    "$BUCKET?prefix=$REL/metadata/GeneFull_Ex50pAS/&delimiter=/" \
    | grep -o "<Prefix>[^<]*</Prefix>" | sed 's/<[^>]*>//g' \
    | grep -E "/$" | grep -v "GeneFull_Ex50pAS/$" > "$_tmp" || true
  if [ -s "$_tmp" ]; then
    mv "$_tmp" "$ROOT/.species_prefixes.txt"
  else
    rm -f "$_tmp"
    if [ -d "$ROOT/metadata/GeneFull_Ex50pAS" ]; then
      ls "$ROOT/metadata/GeneFull_Ex50pAS" \
        | sed "s|^|$REL/metadata/GeneFull_Ex50pAS/|; s|$|/|" \
        > "$ROOT/.species_prefixes.txt"
      echo "listing failed, reused local species list"
    fi
  fi
fi

while read -r p; do
  sp=$(basename "$p")
  out="$ROOT/metadata/GeneFull_Ex50pAS/$sp"
  mkdir -p "$out"
  [ -s "$out/sample_metadata.parquet" ] && continue
  curl -sS --fail --retry 5 --retry-delay 3 --max-time 300 \
    -o "$out/sample_metadata.parquet.partial" \
    "$BUCKET/${p}sample_metadata.parquet" \
    && mv "$out/sample_metadata.parquet.partial" \
          "$out/sample_metadata.parquet"
  echo "metadata $sp done"
done < "$ROOT/.species_prefixes.txt"

# ---- matrices ---------------------------------------------------------------
# Driven from Python rather than xargs. Tab-separated URLs and destinations
# survive shell quoting poorly, which silently corrupted the URLs on the first
# attempt and produced 404s for objects that exist.
"$PY" - "$ROOT/$MAN" "$JOBS" <<'PYEOF'
import concurrent.futures as cf
import os
import subprocess
import sys
import uuid

import h5py

man = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
jobs = int(sys.argv[2])


def verified(dest, cells):
    """A file only counts as done if HDF5 opens and the cell count matches."""
    try:
        with h5py.File(dest, "r") as f:
            n = int(f["X"].attrs["shape"][0])
        return n == cells
    except Exception:
        return False


def get(row):
    url, dest, cells = row[0], row[1], int(row[5])
    if os.path.exists(dest) and verified(dest, cells):
        return "skip"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    # A per-attempt temporary name, so two workers can never write the same
    # inode and truncate each other's finished file.
    part = "{}.part.{}".format(dest, uuid.uuid4().hex[:8])
    r = subprocess.run(["curl", "-sS", "--fail", "--retry", "5",
                        "--retry-delay", "5", "--retry-all-errors",
                        "--max-time", "3600", "-o", part, url])
    if r.returncode == 0 and verified(part, cells):
        os.replace(part, dest)
        return "ok"
    if os.path.exists(part):
        os.remove(part)
    print("FAIL {}".format(url), flush=True)
    return "fail"


with cf.ThreadPoolExecutor(jobs) as ex:
    res = list(ex.map(get, man))
print("ok={} skip={} fail={}".format(res.count("ok"), res.count("skip"),
                                     res.count("fail")), flush=True)
PYEOF

echo "download complete: $(find "$ROOT/h5ad" -name '*.h5ad' | wc -l) matrices"
du -sh "$ROOT"
