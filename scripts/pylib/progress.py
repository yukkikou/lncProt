#!/usr/bin/env python
"""Progress ledger, write-guard, failure recording, and step context manager.

Every pipeline script uses:
    from pylib.progress import step, safe_open, safe_write_path
    with step("04_classify_loci", "genes=63086 classes=5") as st:
        ... st.note("divergent=1234") ...
LEDGER: $LOGS/PROGRESS.md (appended; markdown table)
Failure -> FAILED row + stderr tail to $LOGS/failures/<step>.err + marker file.
"""
import os
import sys
import time
import traceback
import datetime
from contextlib import contextmanager

WORKDIR = "/media/share/workdir/1610305236/lncProt"
LOGS = os.path.join(WORKDIR, "logs")
RESULTS = os.path.join(WORKDIR, "results")
LEDGER = os.path.join(LOGS, "PROGRESS.md")
FAILDIR = os.path.join(LOGS, "failures")

_FORBIDDEN_PREFIXES = ("/media/share/resource",)  # never writable


def _now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def safe_write_path(path):
    """Assert that a write target is inside WORKDIR and not inside the read-only
    reference tree. Hard guarantee that the reference is never modified."""
    ap = os.path.abspath(path)
    if ap.startswith(_FORBIDDEN_PREFIXES):
        raise PermissionError(f"write into read-only reference refused: {ap}")
    if not ap.startswith(WORKDIR + "/"):
        raise PermissionError(f"write outside WORKDIR refused: {ap} (WORKDIR={WORKDIR})")
    os.makedirs(os.path.dirname(ap), exist_ok=True)
    return ap


def _ledger_row(step, status, nums=""):
    os.makedirs(LOGS, exist_ok=True)
    new = not os.path.exists(LEDGER)
    with open(LEDGER, "a") as fh:
        if new:
            fh.write("| step | status | time | key numbers |\n|---|---|---|---|\n")
        fh.write(f"| {step} | {status} | {_now()} | {nums} |\n")


def df_to_md(df, floatfmt=".3f", max_rows=40):
    """Plain-text markdown table (no `tabulate` dependency)."""
    import pandas as pd
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else format(v, floatfmt))
    d = d.head(max_rows)
    cols = [str(c) for c in d.columns]
    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, row in d.iterrows():
        lines.append("| " + " | ".join(str(v) for v in row.values) + " |")
    return "\n".join(lines)


class _Step:
    def __init__(self, name):
        self.name = name
        self.notes = []
        self.t0 = time.time()

    def note(self, text):
        self.notes.append(str(text))
        print(f"[{self.name}] {text}", flush=True)

    def summary_md(self, path, title, body):
        """Write per-step result summary (markdown) under results/step_summaries/."""
        path = safe_write_path(path)
        with open(path, "w") as fh:
            fh.write(f"# {title}\n\n_Generated {_now()} by {self.name}_\n\n{body}\n")
        self.note(f"summary -> {path}")


@contextmanager
def step(name, done_nums=""):
    """Context manager: RUNNING row on enter, DONE/FAILED row on exit."""
    st = _Step(name)
    _ledger_row(name, "RUNNING")
    print(f">>> {name} start {_now()}", flush=True)
    try:
        yield st
    except Exception:
        etype, evalue, tb = sys.exc_info()
        tail = "\n".join(traceback.format_exception(etype, evalue, tb)[-8:])
        _ledger_row(name, "FAILED", tail.replace("|", "/").replace("\n", " ")[:300])
        os.makedirs(FAILDIR, exist_ok=True)
        with open(os.path.join(FAILDIR, f"{name}.err"), "w") as fh:
            fh.write(tail + "\n")
        marker = safe_write_path(os.path.join(RESULTS, f".{name}.failed"))
        open(marker, "w").close()
        raise
    dt = time.time() - st.t0
    nums = done_nums + ("" if not st.notes else " | " + "; ".join(st.notes))
    _ledger_row(name, f"DONE({dt:.0f}s)", nums)
    marker = os.path.join(RESULTS, f".{name}.failed")
    if os.path.exists(marker):
        os.remove(marker)
    print(f"<<< {name} done in {dt:.0f}s {_now()}", flush=True)
