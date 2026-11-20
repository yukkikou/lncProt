#!/usr/bin/env python
"""34_comparative.py — comparative-genome analyses (Sun et al. inspiration).

Modes:
  E  Per-lncRNA transcript age by orthologous TSS presence in macaque /
     mouse / rat (chain liftover + target GTF transcript TSS within 1 kb).
  G  Human-vs-macaque promoter substitution rate from chain-aligned blocks
     (whole promoter and TSS +/-50 core), by birth class.
Outputs -> results/comparative/*.tsv + step summaries.
"""
import gzip
import os
import sys
from bisect import bisect_right
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR
from pylib.chainlift import ChainDB

RES = os.path.join(WORKDIR, "results")
CMP = os.path.join(WORKDIR, "data", "comparative")
OUT = os.path.join(RES, "comparative")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]
HUMAN_FA = ("/media/share/resource/reference/homo_sapiens/GRCh38/"
            "genome/GRCh38.CHR.fa")
SPECIES = [("macaque", 1, "chains/hg38ToRheMac10.over.chain.gz",
            "gtf/Macaca_mulatta.Mmul_10.116.gtf.gz"),
           ("mouse", 2, "chains/hg38ToMm39.over.chain.gz",
            "gtf/Mus_musculus.GRCm39.116.gtf.gz"),
           ("rat", 3, "chains/hg38ToRn7.over.chain.gz",
            "gtf/Rattus_norvegicus.mRatBN7.2.110.gtf.gz")]
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def mwu(y, x):
    u, p = mannwhitneyu(y, x, alternative="two-sided")
    return float(2 * u / (len(y) * len(x)) - 1), float(p)


def open_chain(path):
    if path.endswith(".gz"):
        plain = path[:-3]
        if not os.path.exists(plain):
            with gzip.open(path, "rt") as fi, open(plain, "w") as fo:
                for line in fi:
                    fo.write(line)
        return plain
    return path


def load_tss(gtf_path):
    d = defaultdict(lambda: {"+": [], "-": []})
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            chrom = f[0].replace("chr", "")
            start, end, strand = int(f[3]), int(f[4]), f[6]
            d[chrom][strand].append(start if strand == "+" else end)
    for c in d:
        for s in ("+", "-"):
            d[c][s] = np.array(sorted(d[c][s]))
    return d


def present(db, chrom, pos, strand=None, win=1000, either=True):
    arrs = []
    for s in (("+", "-") if either else (strand,)):
        a = db.get(chrom.replace("chr", ""), {}).get(s)
        if a is not None and len(a):
            arrs.append(a)
    for arr in arrs:
        i = np.searchsorted(arr, pos)
        for j in (i - 1, i):
            if 0 <= j < len(arr) and abs(int(arr[j]) - pos) <= win:
                return True
    return False


def get_cands(db, chrom):
    for key in (chrom, db._canon(chrom), "chr" + db._canon(chrom)):
        if key in db.by_tname:
            return db.by_tname[key]
    return None


def lift_point(db, chrom, pos):
    """Point liftover (target->query) using a direct chain scan."""
    cands = get_cands(db, chrom)
    if not cands:
        return None
    for c in cands:
        if not (c.t_start <= pos < c.t_end):
            continue
        i = bisect_right(c.block_t, pos) - 1
        if i < 0:
            continue
        bt, bq, size = c.blocks[i]
        off = pos - bt
        if off < 0 or off >= size:
            continue
        if c.q_strand == "+":
            qpos = bq + off
        else:
            qpos = c.q_size - (bq + off) - 1
        return c.q_name, qpos, c.q_strand
    return None


# ------------------------------------------------------------------ mode E --
def mode_e():
    with step("34_cmp_E") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        res = pd.DataFrame(index=meta.index)
        for sp, level, chain_rel, gtf_rel in SPECIES:
            cpath = os.path.join(CMP, chain_rel)
            gpath = os.path.join(CMP, gtf_rel)
            if not (os.path.exists(cpath) and os.path.exists(gpath)):
                st.note("{} missing; skipped".format(sp))
                continue
            db = ChainDB(open_chain(cpath))
            tss_db = load_tss(gpath)
            vals = []
            for chrom, tss, strand in zip(meta.chrom, meta.tss,
                                          meta.strand):
                hit = False
                if pd.notna(tss) and pd.notna(chrom):
                    pt = int(tss) - 1
                    lf = lift_point(db, str(chrom), pt)
                    if lf is not None:
                        qname, qpos, qstrand = lf
                        tstrand = "+" if qstrand == "+" else "-"
                        hit = present(tss_db, qname, qpos, tstrand, either=False)
                vals.append(int(hit))
            res[sp] = vals
            st.note("{}: present={}/{} ({:.1%})".format(
                sp, int(np.sum(vals)), len(vals),
                np.mean(vals)))
        levels = [c for c in ["macaque", "mouse", "rat"] if c in res.columns]
        res["age_level"] = 0
        for c in levels:
            res["age_level"] = np.maximum(res["age_level"],
                                          res[c] * dict(
                                              (s[0], s[1]) for s in SPECIES)[c])
        lab = {0: "human_specific", 1: "primate", 2: "mammalian",
               3: "eutherian"}
        res["age_class"] = res["age_level"].map(lab)
        res.to_csv(safe_write_path(os.path.join(
            OUT, "E_transcript_age.tsv")), sep="\t")
        m = meta[["cls", "group"]].join(res)
        rows = []
        for c in ["human_specific", "primate", "mammalian", "eutherian"]:
            rows.append(dict(age_class=c, n_lnc=int(
                ((m.group == "lnc") & (m.age_class == c)).sum()),
                frac_lnc=float(((m.group == "lnc") &
                                (m.age_class == c)).sum()
                               / max((m.group == "lnc").sum(), 1)),
                n_PCG=int(((m.group == "pcg") & (m.age_class == c)).sum()),
                frac_PCG=float(((m.group == "pcg") &
                                (m.age_class == c)).sum()
                               / max((m.group == "pcg").sum(), 1))))
        A = pd.DataFrame(rows)
        # class x age
        tab = pd.crosstab(m[m.group == "lnc"].cls,
                          m[m.group == "lnc"].age_class)
        tab.to_csv(safe_write_path(os.path.join(
            OUT, "E_class_x_age.tsv")), sep="\t")
        from scipy.stats import chi2_contingency
        chi2, pchi, _, _ = chi2_contingency(tab.values)
        ph = pd.read_csv(os.path.join(RES, "phylop", "phylop_per_gene.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["phylop_mean"]]
        mm = m.join(ph)
        r, pr = spearmanr(mm.age_level, mm.phylop_mean, nan_policy="omit")
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["tss_in_ltr"]]
        mm = mm.join(te)
        ig = mm[(mm.cls == "intergenic") & mm.tss_in_ltr.notna()]
        a = ig.loc[ig.tss_in_ltr == 1, "age_level"]
        b = ig.loc[ig.tss_in_ltr != 1, "age_level"]
        dlt, pl = mwu(a.values, b.values)
        st.note("E age: class x age chi2 p={:.1e}; age~phyloP rho={:.3f}; "
                "intergenic LTR-vs-nonLTR age d={:.3f} p={:.1e}".format(
                    pchi, r, dlt, pl))
        st.summary_md("results/comparative/step_summaries/34_e.md",
                      "34E transcript age by orthologous TSS",
                      "## Lnc vs PCG age distribution\n\n"
                      + df_to_md(A, ".3f")
                      + "\n\n## Class x age (lnc)\n\n"
                      + df_to_md(tab.reset_index(), "d")
                      + "\n\nage~phyloP rho={:.3f}, P={:.2e}".format(r, pr))


# ------------------------------------------------------------------ mode G --
def find_chain(db, chrom, s, e):
    cands = get_cands(db, chrom)
    if not cands:
        return None
    for c in cands:
        if c.t_start <= s and e - 1 < c.t_end:
            return c
    return None


def blocks_for(chain, s, e):
    out = []
    i0 = max(0, bisect_right(chain.block_t, s) - 1)
    for bt, bq, size in chain.blocks[i0:]:
        if bt >= e:
            break
        bs, be = max(bt, s), min(bt + size, e)
        if be <= bs:
            continue
        off, L = bs - bt, be - bs
        if chain.q_strand == "+":
            qs, qe, strand = bq + off, bq + off + L, "+"
        else:
            qs = chain.q_size - (bq + off + L)
            qe = chain.q_size - (bq + off)
            strand = "-"
        out.append((bs, be, qs, qe, strand))
    return out


def fetch(fa, chrom, s, e):
    for name in (chrom, chrom.replace("chr", ""),
                 "chr" + chrom.replace("chr", "")):
        try:
            if name in fa.references:
                return fa.fetch(name, s, e)
        except Exception:
            pass
    return None


def mode_g():
    with step("34_cmp_G") as st:
        prom = pd.read_csv(os.path.join(RES, "promoters_feat.bed"), sep="\t",
                           names=["chrom", "start", "end", "gene_id", "x",
                                  "strand"])
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", index_col=0, low_memory=False)
        cpath = os.path.join(CMP, "chains",
                             "hg38ToRheMac10.over.chain.gz")
        cand = [os.path.join(CMP, "genome", "rheMac10.fa.gz"),
                os.path.join(CMP, "genome",
                             "Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz")]
        gpath = next((p for p in cand if os.path.exists(p)), None)
        if gpath is None:
            raise RuntimeError("macaque genome missing: " + "; ".join(cand))
        fa_path = gpath[:-3]
        if not os.path.exists(fa_path):
            with gzip.open(gpath, "rb") as fi, open(fa_path, "wb") as fo:
                while True:
                    b = fi.read(1 << 24)
                    if not b:
                        break
                    fo.write(b)
        if not os.path.exists(fa_path + ".fai"):
            pysam.faidx(fa_path)
        db = ChainDB(open_chain(cpath))
        hfa = pysam.FastaFile(HUMAN_FA)
        mfa = pysam.FastaFile(fa_path)
        tss = meta["tss"].reindex(prom.gene_id).values
        rows = []
        for i, (chrom, s, e, gid) in enumerate(zip(
                prom.chrom, prom["start"], prom["end"], prom.gene_id)):
            ch = find_chain(db, chrom, s, e)
            if ch is None:
                rows.append((gid, np.nan, np.nan, np.nan, np.nan, 0, 0))
                continue
            hs = hfa.fetch(chrom, s, e) if chrom in hfa.references else None
            if hs is None:
                rows.append((gid, np.nan, np.nan, np.nan, np.nan, 0, 0))
                continue
            subs = aln = csubs = caln = 0
            scpg = lcpg = snon = lnon = 0
            core_s = core_e = None
            if pd.notna(tss[i]):
                core_s, core_e = int(tss[i]) - 51, int(tss[i]) + 49
            for bs, be, qs, qe, strand in blocks_for(ch, s, e):
                ms = fetch(mfa, ch.q_name, qs, qe)
                if ms is None:
                    continue
                if strand == "-":
                    ms = ms.translate(RC)[::-1]
                hb = hs[bs - s:be - s].upper()
                mb = ms.upper()
                if len(hb) != len(mb):
                    continue
                a_mask = np.array([x in "ACGT" for x in hb]) & \
                    np.array([x in "ACGT" for x in mb])
                if a_mask.sum() == 0:
                    continue
                aln += int(a_mask.sum())
                hh = np.frombuffer(hb.encode(), dtype="S1")
                mm_ = np.frombuffer(mb.encode(), dtype="S1")
                diff = (hh != mm_) & a_mask
                subs += int(diff.sum())
                # CpG-aware split (dinucleotide in the human sequence)
                cpg = np.zeros(len(hh), dtype=bool)
                if len(hh) > 1:
                    is_c = hh[:-1] == b"C"
                    is_g = hh[1:] == b"G"
                    cpg[:-1] |= is_c & is_g
                    cpg[1:] |= is_c & is_g
                cm = cpg & a_mask
                scpg += int((diff & cm).sum())
                lcpg += int(cm.sum())
                snon += int((diff & ~cm).sum())
                lnon += int((~cpg & a_mask).sum())
                if core_s is not None:
                    lo, hi = max(core_s, bs), min(core_e, be)
                    if hi > lo:
                        cs = lo - bs
                        ce = hi - bs
                        cmask = a_mask[cs:ce]
                        caln += int(cmask.sum())
                        csubs += int(((hh[cs:ce] != mm_[cs:ce]) &
                                      cmask).sum())
            rows.append((gid, subs, aln, csubs, caln,
                         scpg, lcpg, snon, lnon,
                         1 if aln > 0 else 0,
                         len(blocks_for(ch, s, e))))
        T = pd.DataFrame(rows, columns=["gene_id", "subs", "aln", "csubs",
                                        "caln", "subs_cpg", "len_cpg",
                                        "subs_non", "len_non", "mapped",
                                        "n_blocks"])
        T["rate_all"] = T.subs / T.aln.replace(0, np.nan)
        T["rate_core"] = T.csubs / T.caln.replace(0, np.nan)
        T["rate_cpg"] = T.subs_cpg / T.len_cpg.replace(0, np.nan)
        T["rate_non"] = T.subs_non / T.len_non.replace(0, np.nan)
        T = T.set_index("gene_id")
        T.to_csv(safe_write_path(os.path.join(
            OUT, "G_promoter_substitution.tsv")), sep="\t")
        m = meta[["cls", "group"]].join(T)
        rows = []
        pcg = m[m.group == "pcg"]
        for feat in ["rate_all", "rate_core", "rate_cpg", "rate_non"]:
            x = pcg[feat].dropna().values
            for cls in CLASSES + ["pseudogene", "transcribed_pseudogene"]:
                y = m.loc[m.cls == cls, feat].dropna().values
                if len(y) < 20:
                    continue
                dlt, p = mwu(y, x)
                rows.append(dict(feature=feat, cls=cls, n=len(y),
                                 median_class=float(np.median(y)),
                                 median_PCG=float(np.median(x)),
                                 cliffs_delta=dlt, pval=p))
        ig = m[m.cls == "intergenic"]
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"),
                         sep="\t", index_col=0, low_memory=False)[
            ["tss_in_ltr"]]
        ig = ig.join(te)
        for feat in ["rate_all", "rate_core", "rate_cpg", "rate_non"]:
            a = ig.loc[ig.tss_in_ltr == 1, feat].dropna().values
            b = ig.loc[ig.tss_in_ltr != 1, feat].dropna().values
            if len(a) >= 20 and len(b) >= 20:
                dlt, p = mwu(a, b)
                rows.append(dict(feature=feat, cls="intergenic_LTR_vs_nonLTR",
                                 n=len(a) + len(b),
                                 median_class=float(np.median(a)),
                                 median_PCG=float(np.median(b)),
                                 cliffs_delta=dlt, pval=p))
        R = pd.DataFrame(rows)
        R["fdr"] = bh_fdr(R.pval.values)
        R.to_csv(safe_write_path(os.path.join(
            OUT, "G_substitution_tests.tsv")), sep="\t", index=False)
        st.note("G mapped={}/{}; median rate all={:.4f}; LTR-vs-nonLTR "
                "d={:.3f} p={:.1e}".format(
                    int(T.mapped.sum()), len(T),
                    float(T.rate_all.median()),
                    R.iloc[-1].cliffs_delta, R.iloc[-1].pval))
        st.summary_md("results/comparative/step_summaries/34_g.md",
                      "34G human-macaque promoter substitution",
                      "## Class vs PCG\n\n" + df_to_md(R, ".4f")
                      + "\n\nMapped promoters: {}/{}".format(
                          int(T.mapped.sum()), len(T)))


if __name__ == "__main__":
    mode = (sys.argv[1] if len(sys.argv) > 1 else "").upper()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "step_summaries"), exist_ok=True)
    if mode == "E":
        mode_e()
    elif mode == "G":
        mode_g()
    else:
        raise SystemExit("usage: 34_comparative.py E|G")
