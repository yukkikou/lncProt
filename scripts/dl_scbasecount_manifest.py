#!/usr/bin/env python
"""dl_scbasecount_manifest.py — verified matrix download from a manifest.

This is the matrix half of scripts/51_scbasecount_download.sh, split out so a
design can be fetched without repeating the per-species metadata step, whose
bucket listing is the fragile part of that script. Behaviour is unchanged. A
file counts as done only when HDF5 opens and the cell count matches the
manifest, so re-running is always safe.

Usage
  python scripts/dl_scbasecount_manifest.py <manifest.tsv> [jobs]
"""
import concurrent.futures as cf
import os
import subprocess
import sys
import uuid

import h5py


def verified(dest, cells):
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
    # per-attempt temporary name so two workers can never share an inode
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


def main():
    man_path = sys.argv[1]
    jobs = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    rows = [l.rstrip("\n").split("\t") for l in open(man_path)][1:]
    rows = [r for r in rows if len(r) >= 6]
    print("manifest {} rows {}".format(os.path.basename(man_path), len(rows)),
          flush=True)
    with cf.ThreadPoolExecutor(jobs) as ex:
        res = list(ex.map(get, rows))
    print("ok={} skip={} fail={}".format(res.count("ok"), res.count("skip"),
                                         res.count("fail")), flush=True)


if __name__ == "__main__":
    main()
