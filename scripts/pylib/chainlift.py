#!/usr/bin/env python
"""chainlift.py — minimal, strict UCSC chain-file interval liftover (pure stdlib).

Algorithm (Kent-style endpoint mapping):
  * chain records: target(t)=hg19, query(q)=hg38 for hg19ToHg38.over.chain
  * blocks give aligned (tStart, qStart, size) runs; q coords are in qStrand frame
  * lift interval [s, e): map s and e-1 through blocks of the SAME chain
      - qStrand '+': out = [q(s), q(e-1)+1)
      - qStrand '-': out = [qSize - q(e-1) - 1, qSize - q(s))
  * strict: endpoints must fall inside aligned blocks of one chain, else drop
  * intervals spanning internal block gaps are kept (endpoint mapping) but
    the number of crossed gaps is reported for QC
Reference chain file is read READ-ONLY.
"""
from bisect import bisect_right
from collections import defaultdict


class Chain:
    __slots__ = ("t_start", "t_end", "q_name", "q_size", "q_strand",
                 "block_t", "blocks", "index")

    def __init__(self, t_start, t_end, q_name, q_size, q_strand):
        self.t_start, self.t_end = t_start, t_end
        self.q_name, self.q_size, self.q_strand = q_name, q_size, q_strand
        self.blocks = []        # list of (t_start, q_start, size)
        self.block_t = []       # sorted t starts for bisect


class ChainDB:
    def __init__(self, chain_path):
        self.by_tname = defaultdict(list)
        self._parse(chain_path)
        for lst in self.by_tname.values():
            lst.sort(key=lambda c: c.t_start)

    @staticmethod
    def _canon(name):
        """UCSC 'chr1' <-> Ensembl '1' normalization."""
        return name[3:] if name.lower().startswith("chr") else name

    def _parse(self, path):
        cur = None
        t_cursor = q_cursor = 0
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    if cur is not None:
                        cur.block_t = [b[0] for b in cur.blocks]
                    cur = None
                    continue
                if line.startswith("#"):
                    continue
                if line.startswith("chain"):
                    p = line.split()
                    # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                    assert p[4] == "+", "target strand must be + in chain"
                    cur = Chain(int(p[5]), int(p[6]), p[7], int(p[8]), p[9])
                    self.by_tname[p[2]].append(cur)
                    t_cursor, q_cursor = int(p[5]), int(p[10])
                    continue
                p = line.split()
                if len(p) == 3 and cur is not None:
                    size, dt, dq = int(p[0]), int(p[1]), int(p[2])
                    cur.blocks.append((t_cursor, q_cursor, size))
                    t_cursor += size + dt
                    q_cursor += size + dq

    @staticmethod
    def _lift_pos(chain, pos):
        """position (target frame) -> q-frame position, or None if in a gap."""
        i = bisect_right(chain.block_t, pos) - 1
        if i < 0:
            return None
        bt, bq, size = chain.blocks[i]
        off = pos - bt
        if off >= size:
            return None
        return bq + off

    @staticmethod
    def _gap_count(chain, s, e):
        """number of internal block boundaries crossed by [s, e)."""
        i = bisect_right(chain.block_t, s) - 1
        j = bisect_right(chain.block_t, e - 1) - 1
        return max(0, j - i)

    def lift(self, chrom, start, end):
        """Lift 0-based half-open interval. Returns dict or None.
        Dict: {chrom, start, end, strand, gaps}. Input chrom may carry a 'chr'
        prefix (UCSC) or not (Ensembl); output keeps the input's convention."""
        if end <= start:
            return None
        had_chr = chrom.lower().startswith("chr")
        cands = self.by_tname.get(self._canon(chrom))
        if not cands:
            return None
        # candidate chains whose t-range covers [start, end-1]
        import bisect
        lo = bisect.bisect_right(cands, start - 1, key=lambda c: c.t_end)  # first chain with t_end > start-1
        mapped, gaps = None, 0
        for c in cands[lo:]:
            if c.t_start > end - 1:
                break
            if not (c.t_start <= start and end - 1 < c.t_end):
                continue
            qs = self._lift_pos(c, start)
            qe = self._lift_pos(c, end - 1)
            if qs is None or qe is None:
                continue
            gaps = self._gap_count(c, start, end)
            qname = c.q_name
            if had_chr and not qname.lower().startswith("chr"):
                qname = "chr" + qname
            if c.q_strand == "+":
                mapped = (qname, qs, qe + 1, "+")
            else:
                ps = c.q_size - qe - 1
                pe = c.q_size - qs
                mapped = (qname, ps, pe, "-")
            break
        if mapped is None:
            return None
        return dict(chrom=mapped[0], start=mapped[1], end=mapped[2],
                    strand=mapped[3], gaps=gaps)
