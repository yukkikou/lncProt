#!/usr/bin/env python
"""87_rampage_matched.py — cell-type-matched RAMPAGE support.

The unmatched RAMPAGE comparison in module 84 asks whether any of 225 peak sets
supports a TSS. A lncRNA expressed in one tissue can only be capped in samples
of that tissue, so the unmatched fraction partly measures expression breadth
rather than promoter use. This module matches every peak set to an organ group
and counts support only in peak sets whose organ matches a GTEx organ in which
the gene is expressed (mean TPM >= 1). The RAMPAGE test then asks, within
samples that should be competent for the gene, whether the TSS is capped.

A second test conditions on the partner. For divergent pairs, the peak sets in
which the protein-coding partner is supported are a set of samples in which the
locus is demonstrably transcribed. If the lncRNA borrows that promoter, its TSS
should also be capped in most of those sets. The same containment is computed
for non-divergent opposite-strand pairs at 1-2 kb spacing, which have no
borrowing claim.

Outputs
  results/rampage_matched/A_matched_support.tsv
  results/rampage_matched/B_pair_containment.tsv
  results/rampage_matched/C_independent_initiation.tsv
  results/rampage_matched/D_organ_map.tsv
"""
import glob
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

RES = os.path.join(WORKDIR, "results")
OUT = os.path.join(RES, "rampage_matched")
BED = os.path.join(WORKDIR, "data", "encode_rampage", "bed")
META = os.path.join(OUT, "rampage_biosample_meta.tsv")
HALF = 250
CLASSES = ["PCG", "divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic", "pseudogene", "transcribed_pseudogene"]

# Biosample term names to organ groups. Cell lines are assigned to the organ
# of origin, because the question is whether the gene can be expressed there.
ORGAN_KEYWORDS = [
    ("thyroid", "thyroid"),
    ("adrenal", "adrenal"),
    ("stomach", "stomach"),
    ("gastroesophageal", "esophagus"),
    ("esophagus", "esophagus"),
    ("colon", "colon"),
    ("peyer", "colon"),
    ("intestine", "colon"),
    ("spleen", "spleen"),
    ("lung", "lung"),
    ("alveolus", "lung"),
    ("bronchus", "lung"),
    ("airway", "lung"),
    ("skin", "skin"),
    ("melanoma", "skin"),
    ("keratinocyte", "skin"),
    ("hair follicle", "skin"),
    ("adipose", "adipose"),
    ("fat pad", "adipose"),
    ("heart", "heart"),
    ("cardiac", "heart"),
    ("atrium", "heart"),
    ("ventricle", "heart"),
    ("liver", "liver"),
    ("hepatocyte", "liver"),
    ("muscle", "muscle"),
    ("myoblast", "muscle"),
    ("myotube", "muscle"),
    ("myometrial", "uterus"),
    ("nerve", "nerve"),
    ("spinal cord", "spinal_cord"),
    ("brain", "brain"),
    ("cortex", "brain"),
    ("cerebell", "brain"),
    ("lobe", "brain"),
    ("diencephalon", "brain"),
    ("neural", "brain"),
    ("neuroblastoma", "brain"),
    ("daoy", "brain"),
    ("a172", "brain"),
    ("m059j", "brain"),
    ("sjcrh30", "muscle"),
    ("uterus", "uterus"),
    ("cervix", "uterus"),
    ("vagina", "uterus"),
    ("fallopian", "uterus"),
    ("ovary", "ovary"),
    ("prostate", "prostate"),
    ("pc-3", "prostate"),
    ("testis", "testis"),
    ("kidney", "kidney"),
    ("metanephros", "kidney"),
    ("glomerular", "kidney"),
    ("mesangial", "kidney"),
    ("proximal tubule", "kidney"),
    ("bladder", "bladder"),
    ("pancreas", "pancreas"),
    ("breast", "breast"),
    ("mammary", "breast"),
    ("artery", "artery"),
    ("aorta", "artery"),
    ("blood", "blood"),
    ("lymphocyte", "blood"),
    ("k562", "blood"),
    ("gm12878", "blood"),
    ("gm23", "blood"),
    ("oci-ly", "blood"),
    ("karpas", "blood"),
    ("fibroblast", "fibroblast"),
    ("ht1080", "fibroblast"),
    ("pituitary", "pituitary"),
    ("salivary", "salivary"),
    ("tongue", "esophagus"),
    ("eye", "eye"),
    ("umbilical cord", "blood"),
]

GTEX_ORGAN = {
    "thyroid": ["Thyroid"],
    "adrenal": ["Adrenal Gland"],
    "stomach": ["Stomach"],
    "esophagus": ["Esophagus - Mucosa", "Esophagus - Muscularis",
                  "Esophagus - Gastroesophageal Junction"],
    "colon": ["Colon - Transverse", "Colon - Sigmoid",
              "Small Intestine - Terminal Ileum"],
    "spleen": ["Spleen"],
    "lung": ["Lung"],
    "skin": ["Skin - Sun Exposed (Lower leg)",
             "Skin - Not Sun Exposed (Suprapubic)"],
    "adipose": ["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)"],
    "heart": ["Heart - Left Ventricle", "Heart - Atrial Appendage"],
    "liver": ["Liver"],
    "muscle": ["Muscle - Skeletal"],
    "nerve": ["Nerve - Tibial"],
    "spinal_cord": ["Brain - Spinal cord (cervical c-1)"],
    "brain": ["Brain - Cerebellum", "Brain - Cortex",
              "Brain - Frontal Cortex (BA9)", "Brain - Hippocampus",
              "Brain - Hypothalamus", "Brain - Putamen (basal ganglia)",
              "Brain - Caudate (basal ganglia)",
              "Brain - Nucleus accumbens (basal ganglia)",
              "Brain - Amygdala", "Brain - Substantia nigra",
              "Brain - Anterior cingulate cortex (BA24)",
              "Brain - Cerebellar Hemisphere"],
    "uterus": ["Uterus", "Cervix - Ectocervix", "Cervix - Endocervix",
               "Fallopian Tube", "Vagina"],
    "ovary": ["Ovary"],
    "prostate": ["Prostate"],
    "testis": ["Testis"],
    "kidney": ["Kidney - Cortex"],
    "bladder": ["Bladder"],
    "pancreas": ["Pancreas"],
    "breast": ["Breast - Mammary Tissue"],
    "artery": ["Artery - Tibial", "Artery - Aorta", "Artery - Coronary"],
    "blood": ["Whole Blood", "Cells - EBV-transformed lymphocytes"],
    "fibroblast": ["Cells - Cultured fibroblasts"],
    "pituitary": ["Pituitary"],
    "salivary": ["Minor Salivary Gland"],
}


def to_organ(name):
    n = str(name).lower()
    for key, org in ORGAN_KEYWORDS:
        if key in n:
            return org
    return None


def load_bed(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            c = f[0].replace("chr", "")
            out.setdefault(c, []).append((int(f[1]), int(f[2])))
    merged = {}
    for c, v in out.items():
        v.sort()
        s, e = [], []
        for a, b in v:
            if s and a <= e[-1]:
                e[-1] = max(e[-1], b)
            else:
                s.append(a)
                e.append(b)
        merged[c] = (np.array(s), np.array(e))
    return merged


def gtex_organ_means(st):
    g = pd.read_csv(os.path.join(WORKDIR, "data", "tissue_disease",
                                 "toil_genes.tsv"), sep="\t")
    s = pd.read_csv(os.path.join(WORKDIR, "data", "tissue_disease",
                                 "toil_samples.tsv"), sep="\t")
    z = np.load(os.path.join(WORKDIR, "data", "tissue_disease",
                             "toil_class_expr.npz"), mmap_mode="r")
    M = np.asarray(z["M"])
    ok = (s.cohort == "GTEx").values
    cols = np.where(ok)[0]
    tis = s.tissue.values[cols]
    ex = {}
    for org, tissues in GTEX_ORGAN.items():
        sel = np.isin(tis, tissues)
        if sel.sum() < 5:
            st.note("organ {} has {} GTEx samples".format(org, sel.sum()))
            continue
        ex[org] = np.nanmean(M[:, cols[sel]], axis=1)
    del M
    E = pd.DataFrame(ex, index=g.gene_id.values)
    E.index = E.index.astype(str).str.split(".").str[0]
    E = E[~E.index.duplicated()]
    return E


def main():
    with step("87_rampage_matched") as st:
        os.makedirs(OUT, exist_ok=True)
        meta = pd.read_csv(META, sep="\t")
        meta["organ"] = meta.biosample.map(to_organ)
        meta.to_csv(safe_write_path(os.path.join(
            OUT, "D_organ_map.tsv")), sep="\t", index=False)
        st.note("peak sets mapped to an organ: {}/{}".format(
            meta.organ.notna().sum(), len(meta)))
        mapped = meta.dropna(subset=["organ"]).reset_index(drop=True)

        m = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"), sep="\t",
                        index_col=0, low_memory=False)
        m["tss"] = pd.to_numeric(m.tss, errors="coerce")
        m["chrom_c"] = m.chrom.astype(str).str.replace("chr", "", regex=False)
        m = m.dropna(subset=["tss"]).reset_index()
        m["base"] = m.gene_id.astype(str).str.split(".").str[0]
        by_chr = {c: np.where(m.chrom_c == c)[0] for c in set(m.chrom_c)}

        E = gtex_organ_means(st)
        st.note("GTEx organ matrices: {}".format(", ".join(E.columns)))
        base_ok = pd.Index(m.base)
        Ecols = list(E.columns)
        E = E.reindex(base_ok).values  # genes x organs

        # hit matrix, genes x peak sets
        H = np.zeros((len(m), len(mapped)), dtype=bool)
        for k, path in enumerate(glob.glob(os.path.join(BED, "*.bed.gz"))):
            acc = os.path.basename(path).split(".")[0]
            if acc not in set(mapped.acc):
                continue
            j = int(mapped.index[mapped.acc == acc][0])
            iv = load_bed(path)
            for c, idx in by_chr.items():
                t = iv.get(c)
                if t is None:
                    continue
                s, e = t
                q = m.tss.values[idx].astype(np.int64)
                jj = np.searchsorted(s, q + HALF, "right") - 1
                valid = jj >= 0
                ok = np.zeros(len(idx), dtype=bool)
                ok[valid] = e[jj[valid]] >= (q[valid] - HALF)
                H[idx[ok], j] = True
        st.note("hit matrix {} x {}".format(*H.shape))

        organs = sorted(set(mapped.organ))
        oidx = {o: i for i, o in enumerate(Ecols)}
        set_org = np.array([oidx.get(o, -1) for o in mapped.organ])

        # expression for each gene in the organ of each peak set
        M = np.zeros_like(H)
        for j, o in enumerate(set_org):
            if o < 0:
                continue
            M[:, j] = E[:, o] >= 1.0
        HM = H & M

        rows = []
        for cls in CLASSES:
            sel = np.where(m.cls.values == cls)[0]
            if len(sel) < 20:
                continue
            ns = M[sel].sum(1)
            nr = HM[sel].sum(1)
            have = ns >= 5
            frac = np.divide(nr[have], ns[have], out=np.zeros(have.sum()),
                             where=ns[have] > 0)
            rows.append(dict(cls=cls, n=len(sel),
                             frac_any_unmatched=float(H[sel].any(1).mean()),
                             n_matched_ge5=int(have.sum()),
                             frac_any_matched=float((nr[have] > 0).mean()),
                             mean_sets_gt100=float(
                                 (frac >= 0.1).mean())))
        A = pd.DataFrame(rows)
        A.to_csv(safe_write_path(os.path.join(
            OUT, "A_matched_support.tsv")), sep="\t", index=False)

        # breadth-matched class contrast. Divergent lncRNAs are broadly
        # expressed, so they see more matched peak sets. Within strata of the
        # number of matched sets, compare each class with intergenic lncRNAs.
        ns_all = M.sum(1)
        any_hit = (HM.sum(1) > 0).astype(int)
        strata = np.digitize(ns_all, [5, 10, 20, 40, 80, 160])
        inter = np.where((m.cls.values == "intergenic") & (ns_all >= 5))[0]
        rows = []
        for cls in ("divergent", "antisense", "intronic",
                    "sense_overlapping", "pseudogene",
                    "transcribed_pseudogene"):
            cas = np.where((m.cls.values == cls) & (ns_all >= 5))[0]
            if len(cas) < 100:
                continue
            num = den = 0.0
            for s in np.unique(strata[cas]):
                a = cas[strata[cas] == s]
                b = inter[strata[inter] == s]
                if len(a) < 10 or len(b) < 10:
                    continue
                fa, fb = any_hit[a].mean(), any_hit[b].mean()
                n1, n2 = len(a), len(b)
                num += fa * (1 - fb) * n1 * n2 / (n1 + n2)
                den += fb * (1 - fa) * n1 * n2 / (n1 + n2)
            rows.append(dict(cls=cls, n_class=len(cas), n_intergenic=len(inter),
                             crude_frac_class=float(any_hit[cas].mean()),
                             crude_frac_intergenic=float(any_hit[inter].mean()),
                             mh_or=float(num / den) if den > 0 else np.nan))
        F = pd.DataFrame(rows)
        F.to_csv(safe_write_path(os.path.join(
            OUT, "F_breadth_matched.tsv")), sep="\t", index=False)

        # pair containment
        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        gcol = pd.Series(np.arange(len(m)), index=m.gene_id.values)
        dp["li"] = dp.lnc_gid.map(gcol)
        dp["pi"] = dp.pcg_gid.map(gcol)
        dp = dp.dropna(subset=["li", "pi"]).astype({"li": int, "pi": int})

        def contain(a, b, a_expr, lab):
            """P(a supported | b supported and a expressed) over peak sets."""
            den = b & a_expr
            num = a & den
            n_b = den.sum(1).astype(float)
            ok = n_b >= 5
            c = num[ok].sum(1) / n_b[ok]
            return dict(subset=lab, n=int(ok.sum()),
                        median_partner_sets=float(np.median(n_b[ok])),
                        mean_containment=float(c.mean()),
                        median_containment=float(np.median(c)),
                        frac_containment_gt_half=float((c > 0.5).mean()))

        AA = np.ones(H.shape, dtype=bool)
        B = [contain(H[dp.li.values], H[dp.pi.values], AA[dp.li.values],
                     "all_sets"),
             contain(H[dp.li.values], H[dp.pi.values], M[dp.li.values],
                     "matched_sets")]

        # orientation split, as in module 85
        genes = pd.read_csv(os.path.join(WORKDIR, "data", "genes",
                                         "genes.tsv"), sep="\t").set_index(
                                             "gene_id")
        d2 = dp.copy()
        tssser = pd.Series(m.tss.values, index=m.gene_id.values)
        tssser = tssser[~tssser.index.duplicated()]
        d2["ltss"] = d2.lnc_gid.map(tssser)
        d2["ptss"] = d2.pcg_gid.map(tssser)
        for c, col in (("lstart", genes.start), ("lend", genes.end),
                       ("pstart", genes.start), ("pend", genes.end)):
            key = c
            src = d2.lnc_gid if c.startswith("l") else d2.pcg_gid
            d2[key] = src.map(col)
        d2 = d2.dropna(subset=["ltss", "ptss", "lstart", "lend", "pstart",
                               "pend"])
        ps = d2.pcg_gid.map(m.set_index("gene_id").strand)
        h2h = np.where(ps == "+", d2.ltss < d2.ptss, d2.ltss > d2.ptss)
        for lab, sub in (("matched_h2h", d2[h2h]), ("matched_t2t", d2[~h2h])):
            B.append(contain(H[sub.li.values], H[sub.pi.values],
                             M[sub.li.values], lab))
        B = pd.DataFrame(B)

        # non-divergent opposite-strand control pairs at 1-2 kb spacing. These
        # have the same geometry but sit outside the 1 kb class window, so any
        # borrowing signal that survives here is not specific to the class.
        pg = m[m.group == "pcg"]
        cl, cp, cdist = [], [], []
        for c in set(m.chrom_c):
            L = m[(m.group == "lnc") & (m.chrom_c == c)]
            P = pg[pg.chrom_c == c]
            if not len(L) or not len(P):
                continue
            for s in ("+", "-"):
                Psub = P[P.strand == s]
                Lsub = L[L.strand != s]
                if not len(Psub) or not len(Lsub):
                    continue
                o = np.argsort(Psub.tss.values)
                pv = Psub.tss.values[o]
                pidx = Psub.index.values[o]
                q = Lsub.tss.values
                j = np.searchsorted(pv, q)
                for k, pos in enumerate(q):
                    best, bi = np.inf, -1
                    for jj in (j[k] - 1, j[k]):
                        if 0 <= jj < len(pv):
                            d = abs(pv[jj] - pos)
                            if d < best:
                                best, bi = d, pidx[jj]
                    if 1001 <= best <= 2000:
                        cl.append(int(Lsub.index[k]))
                        cp.append(int(bi))
                        cdist.append(best)
        cl = np.array(cl, dtype=int)
        cp = np.array(cp, dtype=int)
        st.note("control pairs at 1-2 kb: {}".format(len(cl)))
        if len(cl) >= 50:
            B = pd.concat([B, pd.DataFrame([contain(
                H[cl], H[cp], M[cl], "control_1_2kb")])],
                ignore_index=True)
            cdist = np.array(cdist)
            for lo, hi in ((1001, 1250), (1251, 1500), (1501, 2000)):
                k = (cdist >= lo) & (cdist <= hi)
                if k.sum() >= 20:
                    B = pd.concat([B, pd.DataFrame([contain(
                        H[cl[k]], H[cp[k]], M[cl[k]],
                        "control_{}_{}bp".format(lo, hi))])],
                        ignore_index=True)
            # per-pair containment, divergent versus 1-2 kb control
            def per_pair(a, b, a_expr):
                den = b & a_expr
                n = den.sum(1)
                num = (a & den).sum(1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    return np.where(n >= 5, num / np.maximum(n, 1), np.nan)
            cd = per_pair(H[dp.li.values], H[dp.pi.values], M[dp.li.values])
            cc = per_pair(H[cl], H[cp], M[cl])
            cd, cc = cd[~np.isnan(cd)], cc[~np.isnan(cc)]
            from scipy.stats import mannwhitneyu
            u, p = mannwhitneyu(cd, cc, alternative="greater")
            E = pd.DataFrame([dict(
                n_divergent=len(cd), n_control=len(cc),
                mean_divergent=float(cd.mean()), mean_control=float(cc.mean()),
                mwu_U=float(u), mwu_p_greater=float(p),
                cliffs_delta=float(2 * u / (len(cd) * len(cc)) - 1))])
            E.to_csv(safe_write_path(os.path.join(
                OUT, "E_containment_test.tsv")), sep="\t", index=False)
            st.note("containment {} v {}: mean {:.2f} v {:.2f}, "
                    "MWU P={:.2g}, Cliff delta={:.2f}".format(
                        "divergent", "control", cd.mean(), cc.mean(), p,
                        E.cliffs_delta.iloc[0]))

        # independent initiation: peak sets in which one TSS is supported and
        # the other is not, and matched sets in which the partner is not
        # expressed but the lncRNA is
        li, pi = dp.li.values, dp.pi.values
        both = H[li] & H[pi]
        lnc_only = H[li] & ~H[pi]
        pcg_only = ~H[li] & H[pi]
        # sets in which the lncRNA is expressed and the partner is not, so a
        # capped lncRNA TSS cannot be host transcription
        no_host = M[li] & ~M[pi]
        indep = H[li] & no_host
        n_nh = no_host.sum(1)
        ok = n_nh >= 5
        C = []
        C.append(dict(subset="divergent",
                      n_pairs=len(dp),
                      lnc_only_sets=int(lnc_only.sum()),
                      pcg_only_sets=int(pcg_only.sum()),
                      both_sets=int(both.sum()),
                      lnc_only_pairs=int((lnc_only.sum(1) > 0).sum()),
                      pcg_only_pairs=int((pcg_only.sum(1) > 0).sum()),
                      n_pairs_host_silent_ge5=int(ok.sum()),
                      frac_lnc_capped_host_silent=float(
                          indep[ok].sum(1).mean() / n_nh[ok].mean())
                      if ok.any() else np.nan))
        Ct = pd.DataFrame(C)
        Ct.to_csv(safe_write_path(os.path.join(
            OUT, "C_independent_initiation.tsv")), sep="\t", index=False)
        B.to_csv(safe_write_path(os.path.join(
            OUT, "B_pair_containment.tsv")), sep="\t", index=False)

        np.save(os.path.join(OUT, "hit_matrix.npy"), H)
        np.save(os.path.join(OUT, "expression_matrix.npy"), M)
        np.save(os.path.join(OUT, "matched_matrix.npy"), HM)

        for r in A.itertuples():
            st.note("{}: unmatched {:.1%}, matched (>=5 sets) {:.1%} any, "
                    "{:.1%} in >=10% sets, n={}".format(
                        r.cls, r.frac_any_unmatched, r.frac_any_matched,
                        r.mean_sets_gt100, r.n_matched_ge5))
        for r in F.itertuples():
            st.note("{}: crude {:.1%} v intergenic {:.1%}, MH OR {:.2f}"
                    .format(r.cls, r.crude_frac_class,
                            r.crude_frac_intergenic, r.mh_or))
        for r in B.itertuples():
            st.note("{}: n={} median partner sets={:.0f} containment={:.1%}"
                    .format(r.subset, r.n, r.median_partner_sets,
                            r.mean_containment))
        st.note("independent initiation: lnc-only {}/{} pairs, pcg-only {}/{}"
                .format(int(Ct.lnc_only_pairs.iloc[0]), len(dp),
                        int(Ct.pcg_only_pairs.iloc[0]), len(dp)))
        st.summary_md(os.path.join(
            RES, "step_summaries", "87_rampage_matched.md"),
            "87 cell-type-matched RAMPAGE support",
            df_to_md(A, ".4f") + "\n\n" + df_to_md(B, ".4f") + "\n\n"
            + df_to_md(Ct, ".0f"))


if __name__ == "__main__":
    main()
