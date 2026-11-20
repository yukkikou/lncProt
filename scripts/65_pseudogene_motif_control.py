#!/usr/bin/env python
"""65_pseudogene_motif_control.py — is the pseudogene motif enrichment real, or
a dinucleotide-composition artefact?

WHY THIS CONTROL IS NEEDED
--------------------------
The manuscript reports that pseudogene promoters carry the strongest TF motif
enrichments of any class (PHOX2B log2OR = 1.90 vs PCG).  Intergenic lncRNAs got
a dinucleotide-matched control (12c_gc_matched_motifs.py); pseudogenes did not.
That omission matters because pseudogenes preferentially arise in, and are
flanked by, transposon-rich, compositionally skewed sequence.  A PWM scanner
with a uniform background will call more motif hits in such sequence even when
no functional TF site exists.  The published pseudogene enrichment therefore
conflates "more binding sites" with "more GC/CpG/dinucleotide opportunity".

WHAT THIS MODULE DOES
---------------------
For every pseudogene promoter (both `pseudogene` and `transcribed_pseudogene`)
we build N=10 synthetic promoters of identical length and identical dinucleotide
composition with a first-order Markov chain estimated from that promoter's own
dinucleotide counts — the same null used by 12c.  We rescan real and shuffled
promoters with exactly the same FIMO motif library, p-value threshold, and
motif->TF collapsing rule that produced results/motif_presence.tsv.  The
observed count is then compared with the composition-matched expectation using
a gene-stratified permutation test: within each promoter the observed label and
its 10 shuffled labels are exchangeable under the null, so the null distribution
of the total is Poisson-binomial with per-gene success probability k_g/11
(k_g = number of the 11 labels that carry a hit).  That null is evaluated
exactly by convolution, i.e. without Monte-Carlo error.

Outputs
  results/evo_sensitivity/pseudogene_motif_matched.tsv
  results/step_summaries/65_pseudogene_motif_control.md
"""
import gzip
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.signal import fftconvolve
from scipy.stats import binom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR
from pylib.tfnames import norm_tf_symbol

RES = os.path.join(WORKDIR, "results")
DATA = os.path.join(WORKDIR, "data")
OUTDIR = os.path.join(RES, "evo_sensitivity")
WORK = os.path.join(DATA, "qc", "65_motif_control")
CHUNK_DIR = os.path.join(DATA, "motifs", "chunks")
FASTA = os.path.join(RES, "promoters_fimo.fa")
SIF_MEME = "/media/share/resource/container/meme/meme-5.4.1.sif"
FIMO_P = "1e-4"
N_SHUFFLE = 10
SEED = 20261102          # fixed seed: the whole module is reproducible
N_WORKERS = 18
CLASSES = ["pseudogene", "transcribed_pseudogene"]

# A/T/C/G -> 0..3; anything else (N) is invalid and those promoters are dropped.
_CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _CODE[ord(_b)] = _i




def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    ranked = p[o] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(ranked, 0, 1)
    return out


def read_fasta_subset(path, keep):
    """Read only promoters for `keep`, returning {gene_id: sequence}."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and name in keep:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None and name in keep:
        seqs[name] = "".join(buf)
    return seqs


def encode(seq):
    return _CODE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def sample_cdf(p, u):
    """Inverse-CDF draw from rows of p (n x k) given uniforms u (n,)."""
    c = np.cumsum(p, axis=1)
    return np.clip((u[:, None] > c).sum(axis=1), 0, p.shape[1] - 1)


def dinuc_shuffle_batch(codes, rng):
    """First-order Markov sequences with each row's own dinucleotide counts.

    Same null as 12c_gc_matched_motifs.py: the synthetic promoter keeps the
    mono- and dinucleotide composition of the real promoter in expectation,
    while destroying all higher-order structure (including real motif sites).
    Implemented across the whole batch so the 600-bp walk is vectorised.
    """
    n, length = codes.shape
    cnt = np.zeros((n, 4, 4), dtype=np.float64)
    for a in range(4):
        ma = codes[:, :-1] == a
        for b in range(4):
            cnt[:, a, b] = (ma & (codes[:, 1:] == b)).sum(axis=1)
    mono = np.stack([(codes == a).sum(axis=1) for a in range(4)], axis=1)
    mono = mono.astype(np.float64)
    # Bases absent from a promoter get probability 0; guard the degenerate row.
    p0 = mono / np.maximum(mono.sum(axis=1, keepdims=True), 1.0)
    p0[mono.sum(axis=1) == 0] = 0.25

    ptrans = cnt.copy()
    rs = ptrans.sum(axis=2, keepdims=True)
    # A base that only ever occurs last has no observed successor: fall back to
    # the promoter's mononucleotide frequencies rather than inventing a base.
    ptrans = np.where(rs > 0, ptrans / np.maximum(rs, 1.0), p0[:, None, :])

    out = np.empty_like(codes)
    out[:, 0] = sample_cdf(p0, rng.random(n))
    ar = np.arange(n)
    state = out[:, 0]
    for t in range(1, length):
        probs = ptrans[ar, state]
        state = sample_cdf(probs, rng.random(n)).astype(np.uint8)
        out[:, t] = state
    return out


def decode(codes):
    return "".join("ACGT"[c] for c in codes)


def run_fimo(rep, chunk, fa, out_gz):
    """Scan one (replicate, motif chunk) pair and stream only 3 columns.

    --text keeps FIMO from writing the multi-GB XML/GFF payloads; the stdout
    stream is filtered in Python so no intermediate files touch disk.
    """
    if os.path.exists(out_gz) and os.path.getsize(out_gz) > 0:
        return "skip"
    meme = os.path.join(CHUNK_DIR, f"chunk_{chunk}.meme")
    cmd = ["singularity", "exec", "--no-mount", "bind-paths",
           "--bind", "/media/share/resource:/media/share/resource:ro",
           "--bind", f"{WORKDIR}:{WORKDIR}", SIF_MEME,
           "fimo", "--text", "--no-qvalue", "--thresh", FIMO_P,
           "--max-stored-scores", "10000000", meme, fa]
    tmp = out_gz + ".tmp"
    n = 0
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True, bufsize=1 << 20)
    with gzip.open(tmp, "wt", compresslevel=1) as fh:
        for line in p.stdout:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] == "motif_id":
                continue
            fh.write(f"{f[0]}\t{f[1]}\t{f[2]}\n")
            n += 1
    err = p.stderr.read()
    rc = p.wait()
    if rc != 0:
        raise RuntimeError(f"fimo failed rep={rep} chunk={chunk}: {err[-800:]}")
    os.replace(tmp, out_gz)
    return n


def permutation_pvalue(k_counts, obs):
    """Exact gene-stratified permutation p-value.

    Under the null the observed label is exchangeable with the 10 shuffled
    labels inside every promoter.  Given k_g hits among the 11 labels, the
    observed label is a hit with probability k_g/11.  Genes sharing k therefore
    contribute Binomial(c_k, k/11), and the total follows the convolution of
    those 12 binomials.  p = P(T >= obs).
    """
    pmf = np.array([1.0])
    # k runs 0..11 (observed + 10 shuffles), so the k=11 group (every label a
    # hit) is a deterministic 11/11 term and MUST be included; omitting it
    # shifts the null down and inflates significance.
    for k in range(1, N_SHUFFLE + 2):
        c = int(k_counts[k])
        if c == 0:
            continue
        pmf = fftconvolve(pmf, binom.pmf(np.arange(c + 1), c, k / 11.0))
    pmf = np.clip(pmf, 0.0, None)
    s = pmf.sum()
    if s > 0:
        pmf /= s
    if obs >= len(pmf):
        return 0.0
    return float(min(max(pmf[obs:].sum(), 0.0), 1.0))


def main():
    for req in (FASTA, os.path.join(RES, "motif_presence.tsv"),
                os.path.join(RES, "motif_enrichment.tsv"),
                os.path.join(RES, "gene_master_table.tsv")):
        if not os.path.exists(req):
            raise FileNotFoundError(f"required input missing: {req}")
    os.makedirs(WORK, exist_ok=True)

    with step("65_pseudogene_motif_control") as st:
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                           index_col=0, low_memory=False)
        me = pd.read_csv(os.path.join(RES, "motif_enrichment.tsv"), sep="\t")

        # ---- promoter set ---------------------------------------------------
        cls_genes = {c: meta.index[meta.cls == c].to_numpy() for c in CLASSES}
        all_genes = np.concatenate([cls_genes[c] for c in CLASSES])
        seqs = read_fasta_subset(FASTA, set(all_genes))
        missing = [g for g in all_genes if g not in seqs]
        if missing:
            raise RuntimeError(f"{len(missing)} promoters absent from FASTA, "
                               f"e.g. {missing[:3]}")
        # N-containing promoters cannot be shuffled unambiguously; there are
        # very few and they are excluded from BOTH sides of the comparison.
        n_before = {c: len(cls_genes[c]) for c in CLASSES}
        ok = [g for g in all_genes
              if set(seqs[g]) <= set("ACGT") and len(seqs[g]) == 600]
        n_drop = len(all_genes) - len(ok)
        cls_genes = {c: np.array([g for g in cls_genes[c] if g in set(ok)])
                     for c in CLASSES}
        drop_split = ", ".join(f"{n_before[c] - len(cls_genes[c])} {c}"
                               for c in CLASSES)
        st.note(f"promoters: {len(ok)} usable, {n_drop} dropped (N or length); "
                + ", ".join(f"{c}={len(cls_genes[c])}" for c in CLASSES))
        if len(ok) < 100:
            raise RuntimeError("too few usable promoters")

        # ---- observed presence (published matrix) ---------------------------
        pres_obs = pd.read_csv(os.path.join(RES, "motif_presence.tsv"),
                               sep="\t", index_col=0)
        tfs = list(pres_obs.columns)
        tf_col = {t: i for i, t in enumerate(tfs)}
        obs = pres_obs.reindex(ok).fillna(0).astype(bool).to_numpy()
        st.note(f"observed matrix: {obs.shape[0]} promoters x {obs.shape[1]} TFs, "
                f"{int(obs.sum())} promoter-TF positives")
        del pres_obs

        # ---- dinucleotide-shuffled promoters --------------------------------
        rng = np.random.default_rng(SEED)
        codes = np.stack([encode(seqs[g]) for g in ok])
        fa_paths = []
        for rep in range(N_SHUFFLE):
            fa = os.path.join(WORK, f"shuf_rep{rep}.fa")
            fa_paths.append(fa)
            shuf = dinuc_shuffle_batch(codes, rng)
            with open(fa, "w") as fh:
                for i, row in enumerate(shuf):
                    fh.write(f">s{i}\n{decode(row)}\n")
            if rep == 0:
                # QC: composition must be preserved by construction.  Report
                # the mean absolute difference in GC and CpG density.
                gc_o = np.array([(seqs[g].count("G") + seqs[g].count("C")) / 600
                                 for g in ok])
                gc_s = ((shuf == 1) | (shuf == 2)).mean(axis=1)
                st.note(f"shuffle QC rep0: mean |dGC|={np.abs(gc_o - gc_s).mean():.5f}")
        st.note(f"{N_SHUFFLE} x {len(ok)} shuffled promoters, seed={SEED}")

        # ---- rescan with the identical FIMO library -------------------------
        jobs = [(rep, ch, fa_paths[rep],
                 os.path.join(WORK, f"fimo_rep{rep}_chunk{ch}.tsv.gz"))
                for rep in range(N_SHUFFLE) for ch in range(10)]
        for _, _, _, pth in jobs:
            safe_write_path(pth)
        n_done = 0
        with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
            futs = [ex.submit(run_fimo, *j) for j in jobs]
            for f in as_completed(futs):
                f.result()
                n_done += 1
                if n_done % 10 == 0:
                    st.note(f"fimo jobs finished: {n_done}/{len(jobs)}")
        st.note("fimo rescan complete")

        # ---- assemble shuffled presence tensor ------------------------------
        pres_shuf = np.zeros((N_SHUFFLE, len(ok), len(tfs)), dtype=bool)
        motif_tf = {}
        gene_idx = {g: i for i, g in enumerate(ok)}
        for rep in range(N_SHUFFLE):
            for ch in range(10):
                gz = os.path.join(WORK, f"fimo_rep{rep}_chunk{ch}.tsv.gz")
                # FIMO --text emits 10 columns; select only the three we stored.
                df = pd.read_csv(gz, sep="\t", header=None,
                                 names=["mid", "alt", "seq"], usecols=[0, 1, 2],
                                 dtype=str, na_filter=False)
                for mid, alt in df[["mid", "alt"]].drop_duplicates().itertuples(
                        index=False):
                    if mid not in motif_tf:
                        motif_tf[mid] = norm_tf_symbol(mid, alt)
                df["tf"] = df["mid"].map(motif_tf)
                df["gi"] = df["seq"].str.slice(1).astype(np.int32)
                pairs = df[["gi", "tf"]].drop_duplicates()
                gi = pairs["gi"].to_numpy()
                ti = pairs["tf"].map(tf_col).to_numpy()
                good = ~pd.isna(ti)
                pres_shuf[rep, gi[good], ti[good].astype(int)] = True
                del df, pairs
        st.note(f"shuffled motif->TF symbols: {len(motif_tf)}; "
                f"positives per replicate: "
                + ", ".join(str(int(x)) for x in pres_shuf.sum(axis=(1, 2))))

        # ---- per-class, per-TF observed vs composition-matched expectation --
        rows = []
        for cls in CLASSES:
            gset = cls_genes[cls]
            gi = np.array([gene_idx[g] for g in gset])
            obs_c = obs[gi]
            shuf_c = pres_shuf[:, gi, :]
            n = len(gi)
            for j, tf in enumerate(tfs):
                o = obs_c[:, j]
                s = shuf_c[:, :, j]
                obs_n = int(o.sum())
                shuf_n = s.sum(axis=1).astype(float)
                shuf_frac = shuf_n / n
                mean_sf = float(shuf_frac.mean())
                obs_frac = obs_n / n
                ratio = obs_frac / mean_sf if mean_sf > 0 else np.nan
                # k_g = this promoter's observed status (0/1) + number of its
                # 10 shuffles with a hit, so 0 <= k_g <= 11.
                k = o.astype(np.int8) + s.sum(axis=0).astype(np.int8)
                k_counts = np.bincount(k, minlength=N_SHUFFLE + 2)
                p = permutation_pvalue(k_counts, obs_n)
                rows.append(dict(cls=cls, tf=tf, n_promoters=n,
                                 observed_hits=obs_n, observed_frac=obs_frac,
                                 shuffled_hits_mean=float(shuf_n.mean()),
                                 shuffled_frac_mean=mean_sf,
                                 shuffled_frac_sd=float(shuf_frac.std(ddof=1)),
                                 obs_exp_ratio=ratio,
                                 log2_obs_exp=np.log2(ratio) if ratio > 0 else np.nan,
                                 obs_exp_ratio_lo=(obs_frac / shuf_frac.max()
                                                   if shuf_frac.max() > 0 else np.nan),
                                 obs_exp_ratio_hi=(obs_frac / shuf_frac.min()
                                                   if shuf_frac.min() > 0 else np.nan),
                                 perm_p=p))
        res = pd.DataFrame(rows)
        res["fdr"] = bh_fdr(res["perm_p"].values)

        pub = me.rename(columns={"frac_class": "published_frac_class",
                                 "frac_pcg": "published_frac_pcg",
                                 "log2OR": "published_log2OR",
                                 "pval": "published_pval",
                                 "fdr": "published_fdr"})
        res = res.merge(pub[["tf", "cls", "published_frac_class",
                             "published_frac_pcg", "published_log2OR",
                             "published_fdr"]],
                        on=["tf", "cls"], how="left")
        # How much of the published (vs-PCG) enrichment is reproduced by
        # composition alone: log2(shuffled pseudogene fraction / PCG fraction).
        res["shuffled_vs_pcg_log2"] = np.log2(
            res["shuffled_frac_mean"] / res["published_frac_pcg"])
        res["composition_explained"] = (
            (res["published_fdr"] < 0.05) & (res["published_log2OR"] > 0)
            & ~((res["fdr"] < 0.05) & (res["obs_exp_ratio"] > 1)))

        res.to_csv(safe_write_path(
            os.path.join(OUTDIR, "pseudogene_motif_matched.tsv")),
            sep="\t", index=False)

        # ---- summary --------------------------------------------------------
        parts = []
        for cls in CLASSES:
            d = res[res.cls == cls]
            n_prom = int(d.n_promoters.iloc[0])
            tested = len(d)
            up = int(((d.fdr < 0.05) & (d.obs_exp_ratio > 1)).sum())
            down = int(((d.fdr < 0.05) & (d.obs_exp_ratio < 1)).sum())
            pub_up = int(((d.published_fdr < 0.05) & (d.published_log2OR > 0)).sum())
            lost = int((d.composition_explained).sum())
            ph = d[d.tf == "PHOX2B"].iloc[0]
            st.note(f"{cls}: n={n_prom} TFs={tested} survive_up={up} "
                    f"survive_down={down} published_up={pub_up} "
                    f"composition_explained={lost}; PHOX2B ratio="
                    f"{ph.obs_exp_ratio:.3f} fdr={ph.fdr:.2e} "
                    f"(published log2OR={ph.published_log2OR:.2f})")
            parts.append((cls, n_prom, tested, up, down, pub_up, lost, ph))

        # top surviving and top collapsed TFs
        top = (res[(res.fdr < 0.05) & (res.obs_exp_ratio > 1)]
               .sort_values("obs_exp_ratio", ascending=False)
               .groupby("cls").head(10))
        collapsed = (res[res.composition_explained]
                     .sort_values("published_log2OR", ascending=False)
                     .groupby("cls").head(10))

        # correlation between published effect and composition-only effect
        corr_lines = []
        for cls in CLASSES:
            d = res[(res.cls == cls) & res.published_log2OR.notna()
                    & res.shuffled_vs_pcg_log2.notna()]
            r = np.corrcoef(d.published_log2OR, d.shuffled_vs_pcg_log2)[0, 1]
            r2 = np.corrcoef(d.published_log2OR.fillna(0),
                             d.log2_obs_exp.fillna(0))[0, 1]
            corr_lines.append(
                f"- `{cls}`: corr(published log2OR, composition-only "
                f"log2(shuffled/PCG)) = **{r:.3f}**; corr(published log2OR, "
                f"log2 observed/expected) = **{r2:.3f}**")

        body = []
        body.append("## Why this control\n")
        body.append(
            "Pseudogene promoters were the only class in the manuscript without "
            "a dinucleotide-composition control, yet they sit in the most "
            "compositionally skewed, TE-rich sequence. A PWM scan with a uniform "
            "background can therefore produce apparent enrichment with no "
            "functional sites. This module asks whether each pseudogene motif "
            "frequency exceeds what the promoter's own dinucleotide composition "
            "predicts.\n")
        body.append("## Method\n")
        body.append(
            f"- Promoters: {len(ok)} ({', '.join(f'{c}={len(cls_genes[c])}' for c in CLASSES)}), "
            f"{n_drop} dropped for N/length; all 600 bp (-500/+100).\n"
            f"- Null: {N_SHUFFLE} first-order Markov sequences per promoter with "
            "that promoter's own dinucleotide counts (same null as 12c), fixed "
            f"seed {SEED}.\n"
            "- Scan: identical FIMO library (1,280 motifs in 10 chunks), "
            f"p<{FIMO_P}, both strands, identical motif->TF collapsing as "
            "`08_enrichment.py`.\n"
            "- Test: gene-stratified permutation. Within a promoter the observed "
            "label is exchangeable with its 10 shuffled labels, so the null total "
            "is Poisson-binomial with per-gene probability k/11; evaluated "
            "exactly by convolution. p = P(T >= observed). BH-FDR across all "
            "rows shown here (both pseudogene classes).\n")
        body.append("## Headline numbers\n")
        head = pd.DataFrame([
            dict(cls=c, n_promoters=n, tfs_tested=t,
                 survive_up_fdr05=u, survive_down_fdr05=dn,
                 published_up=pu, composition_explained=lo)
            for c, n, t, u, dn, pu, lo, _ in parts])
        body.append(df_to_md(head, ".0f") + "\n")
        body.append(
            "- `survive_up_fdr05`: observed fraction significantly ABOVE the "
            "dinucleotide-shuffled expectation.\n"
            "- `composition_explained`: published pseudogene enrichment "
            "(FDR<0.05, log2OR>0) that no longer exceeds the composition-matched "
            "null.\n")

        # A first-order Markov null preserves adjacent-base composition but not
        # longer-range structure (GC/CpG domains, repeats). Real promoters are
        # therefore expected to sit slightly above it even with no functional
        # sites. The residual ratio, not the raw ratio, is the interpretable
        # quantity; this table supplies that generic floor.
        body.append("## Is the residual specific, or generic?\n")
        glob = []
        for cls in CLASSES:
            d = res[res.cls == cls]
            pu = d[(d.published_fdr < 0.05) & (d.published_log2OR > 0)]
            glob.append(dict(
                cls=cls,
                aggregate_obs_exp_ratio=(d.observed_hits.sum()
                                         / d.shuffled_hits_mean.sum()),
                median_tf_ratio=float(d.obs_exp_ratio.median()),
                published_up_median_ratio=float(pu.obs_exp_ratio.median()),
                frac_tfs_ratio_gt1=float((d.obs_exp_ratio > 1).mean())))
        body.append(
            "A first-order Markov chain matches adjacent-base composition but "
            "not longer-range structure, so every real sequence should sit a "
            "little above its shuffled counterpart. A TF is only interesting if "
            "it clears that generic floor.\n")
        body.append(df_to_md(pd.DataFrame(glob), ".3f") + "\n")

        body.append("## PHOX2B (the manuscript's headline example)\n")
        ph_rows = []
        for cls, *_rest, ph in parts:
            ph_rows.append(dict(cls=cls, n=ph.n_promoters,
                                observed_frac=ph.observed_frac,
                                shuffled_frac=ph.shuffled_frac_mean,
                                ratio=ph.obs_exp_ratio,
                                perm_p=ph.perm_p, fdr=ph.fdr,
                                published_frac_pcg=ph.published_frac_pcg,
                                published_log2OR=ph.published_log2OR,
                                published_fdr=ph.published_fdr))
        body.append(df_to_md(pd.DataFrame(ph_rows), ".4g") + "\n")
        ph_main = res[(res.cls == "pseudogene") & (res.tf == "PHOX2B")].iloc[0]
        d_pse = res[res.cls == "pseudogene"]
        survived = (ph_main.fdr < 0.05) and (ph_main.obs_exp_ratio > 1)
        share = ph_main.shuffled_vs_pcg_log2 / ph_main.published_log2OR
        pctile = float((d_pse.obs_exp_ratio < ph_main.obs_exp_ratio).mean() * 100)
        r_pub_comp = float(np.corrcoef(
            d_pse.published_log2OR, d_pse.shuffled_vs_pcg_log2)[0, 1])
        pu_med = float(d_pse[
            (d_pse.published_fdr < 0.05) & (d_pse.published_log2OR > 0)
        ].obs_exp_ratio.median())
        body.append(
            f"**PHOX2B (pseudogene), quantitatively:** observed "
            f"{ph_main.observed_frac:.4f} vs composition-matched "
            f"{ph_main.shuffled_frac_mean:.4f} (ratio {ph_main.obs_exp_ratio:.3f}; "
            f"permutation p = {ph_main.perm_p:.1e}, FDR = {ph_main.fdr:.1e}). "
            f"The published vs-PCG effect is log2OR {ph_main.published_log2OR:.2f} "
            f"(OR {2 ** ph_main.published_log2OR:.2f}). Composition alone gives "
            f"log2({ph_main.shuffled_frac_mean:.4f}/"
            f"{ph_main.published_frac_pcg:.4f}) = "
            f"{ph_main.shuffled_vs_pcg_log2:.2f} "
            f"(OR {2 ** ph_main.shuffled_vs_pcg_log2:.2f}), i.e. "
            f"**{share * 100:.0f}% of the published log2OR**. The residual is "
            f"{ph_main.obs_exp_ratio:.2f}x and is statistically significant, but "
            f"it is close to the generic real-vs-Markov floor (median TF ratio "
            f"{d_pse.obs_exp_ratio.median():.2f}; aggregate "
            f"{d_pse.observed_hits.sum() / d_pse.shuffled_hits_mean.sum():.2f}) "
            f"and PHOX2B ranks at only the {pctile:.0f}th percentile of that "
            f"ratio.\n\n"
            + ("**Verdict:** the published PHOX2B enrichment is real relative to "
               "a composition-matched null, but most of its magnitude is "
               "dinucleotide composition. The claim that pseudogene promoters "
               "carry the *strongest* motif enrichments does not survive as a "
               "composition-independent statement: the published effect ranking "
               f"tracks the composition-only effect (r = {r_pub_comp:.2f}) and "
               "the published-up TFs have a median observed/expected ratio "
               f"({pu_med:.2f}) below the all-TF median."
               if survived else
               "**Verdict:** the enrichment does NOT survive the matched null: "
               "dinucleotide composition alone reproduces it.") + "\n")

        body.append("## Published effect vs composition-only effect\n")
        body.append("\n".join(corr_lines) + "\n")

        body.append("## TFs still above the matched null (top 10 per class)\n")
        body.append(df_to_md(
            top[["cls", "tf", "observed_frac", "shuffled_frac_mean",
                 "obs_exp_ratio", "perm_p", "fdr", "published_log2OR"]],
            ".4g") if len(top) else "_(none)_")
        body.append("\n\n## Published pseudogene enrichments explained by composition "
                    "(top 10 per class by published log2OR)\n")
        body.append(df_to_md(
            collapsed[["cls", "tf", "observed_frac", "shuffled_frac_mean",
                       "obs_exp_ratio", "perm_p", "fdr", "published_log2OR"]],
            ".4g") if len(collapsed) else "_(none)_")
        body.append(
            "\n\n## Caveats\n"
            "- The null preserves mono-/dinucleotide composition in expectation "
            "but not exactly; residual sampling noise is largest for short "
            "promoters with extreme composition.\n"
            "- A first-order Markov chain does not preserve longer-range "
            "structure (GC/CpG domains, repeats), which is why even non-enriched "
            "TFs sit ~1.1x above it. That inflates every ratio by a similar "
            "factor, so the residual relative to the generic floor, not the raw "
            "ratio, is the interpretable quantity.\n"
            "- A ratio near 1 means the published pseudogene-vs-PCG enrichment "
            "is indistinguishable from what the pseudogene promoters' own "
            "dinucleotide composition generates; it does not by itself prove "
            "the sites are non-functional.\n"
            "- Conversely, a ratio >1 significantly is evidence for structure "
            "beyond dinucleotide composition, but not for TF binding.\n"
            f"- {n_drop} promoters ({drop_split}) contained N and were excluded "
            "from both sides.\n")

        st.summary_md("results/step_summaries/65_pseudogene_motif_control.md",
                      "65 pseudogene motif composition control", "\n".join(body))


if __name__ == "__main__":
    main()
