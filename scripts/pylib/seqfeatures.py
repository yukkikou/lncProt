#!/usr/bin/env python
"""seqfeatures.py — promoter sequence features: GC, CpG, HCP/LCP, core elements.

Core promoter elements are scanned as degenerate motifs with <=1 mismatch at
canonical offsets relative to the TSS (position 0 in the provided sequence =
-500 upstream; sequence is 5'->3' on the gene strand).
"""
import re

# IUPAC expansion to regex
_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def iupac_regex(pattern):
    return "".join(_IUPAC[c] for c in pattern.upper())


def degenerate_match(seq, pattern, max_mismatch=1):
    """True if `pattern` (IUPAC) occurs in seq with <= max_mismatch mismatches.
    Efficient enough for ~700bp promoters x 55k genes."""
    L = len(pattern)
    rx = iupac_regex(pattern)
    core = re.compile(rx)
    # exact hit first (cheap)
    if core.search(seq):
        return True
    if max_mismatch == 0:
        return False
    # <=1 mismatch: for each position, allow any base at one slot ('.' is the
    # regex wildcard — 'N' would match a literal N!)
    for i in range(L):
        p = iupac_regex(pattern[:i]) + "." + iupac_regex(pattern[i + 1:])
        if re.search(p, seq):
            return True
    return False


def gc_content(seq):
    n = len(seq)
    if n == 0:
        return 0.0
    return (seq.count("G") + seq.count("C")) / n


def cpg_oe(seq):
    """CpG observed/expected (Gardiner-Garden): obs CpG / (fC * fG * n)."""
    n = len(seq)
    if n == 0:
        return 0.0
    c, g = seq.count("C"), seq.count("G")
    cpg = seq.count("CG")
    if c == 0 or g == 0:
        return 0.0
    return cpg / (c * g / n)


def cpg_density(seq):
    """CpG count per bp."""
    return seq.count("CG") / max(len(seq), 1)


def hcp_icp_lcp(gc, oe, cpg_count, length=700):
    """Sharp-style promoter classes on the -500..+200 window:
    HCP: (GC>=0.55 and CpG o/e>=0.60 and CpG count>=15)  [lenient modern cut]
    LCP: (GC<0.42 and CpG o/e<0.40) or CpG count<5
    ICP: everything in between."""
    if gc >= 0.55 and oe >= 0.60 and cpg_count >= 15:
        return "HCP"
    if (gc < 0.42 and oe < 0.40) or cpg_count < 5:
        return "LCP"
    return "ICP"


# ---- core promoter element scan within offset windows ----------------------
# offsets relative to TSS; seq[i] corresponds to offset i - 500 (for 700bp seq)
CORE_ELEMENTS = {
    # name: (IUPAC motif, lo, hi)  window in bp relative to TSS
    "TATA":   ("TATAWAGR", -40, -20),
    "Inr":    ("YYANWYY",  -3, 4),
    "CCAAT":  ("RRCCAATC", -100, -40),
    "GCbox":  ("GGGCGG", -120, -20),
    "BREu":   ("SSCGGCC", -45, -18),
    "DPE":    ("RGWYV", 25, 35),
    "TCT":    ("YYCTTYTY", -2, 6),
    "XCPE1":  ("GSCGCCG", -3, 6),
}


def scan_core_elements(seq, up=500):
    """Return dict element -> 0/1 presence in its canonical window."""
    out = {}
    for name, (motif, lo, hi) in CORE_ELEMENTS.items():
        a = max(0, lo + up)
        b = min(len(seq), hi + up)
        sub = seq[a:b] if b > a else ""
        out[name] = 1 if (sub and degenerate_match(sub, motif, 1)) else 0
    return out
