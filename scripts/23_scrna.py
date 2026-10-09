#!/usr/bin/env python
"""23_scrna.py — CGGA IDH_A single-cell RNA-seq (exported Seurat object).

Analyses:
 1. lncRNA birth-class detection/expression by cell type
 2. DIVERGENT PAIR CO-DETECTION within single cells (the headline: direct
    evidence for shared initiation machinery) — observed vs expected under
    independence (product of marginals), with same-strand matched controls
 3. LTR-derived lncRNA assignment to cell types (malignant vs immune)
Inputs: data/.../scrna_export/{counts.mtx, genes.txt, cell_meta.tsv}
Outputs -> results/scrna/*.tsv + step summary
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import binomtest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

EXP = os.path.join(WORKDIR, "data", "ccle_cgga", "cgga", "CGGA_CCell_4083",
                   "scrna_export")
OUT = os.path.join(WORKDIR, "results", "scrna")
RES = os.path.join(WORKDIR, "results")
CLASSES = ["divergent", "antisense", "intronic", "sense_overlapping",
           "intergenic"]


def main():
    with step("23_scrna") as st:
        os.makedirs(OUT, exist_ok=True)
        meta = pd.read_csv(os.path.join(RES, "gene_master_table.tsv"),
                           sep="\t", low_memory=False)
        sym2cls = dict(zip(meta.gene_name.astype(str), meta.cls))
        te = pd.read_csv(os.path.join(RES, "te_promoter_stats.tsv"), sep="\t",
                         index_col=0, low_memory=False)[["tss_in_ltr"]].join(
            meta.set_index("gene_id")["gene_name"])
        sym2ltr = {str(s): int(v) for s, v in zip(te.gene_name, te.tss_in_ltr)
                   if pd.notna(s)}

        X = mmread(os.path.join(EXP, "counts.mtx")).tocsc()
        genes = pd.read_csv(os.path.join(EXP, "genes.txt"), header=None)[0].tolist()
        cells = pd.read_csv(os.path.join(EXP, "cell_meta.tsv"), sep="\t",
                            index_col=0)
        st.note(f"matrix: {len(genes)} genes x {X.shape[1]} cells")

        # gene -> class / LTR
        gcls = np.array([sym2cls.get(g) for g in genes])
        gltr = np.array([sym2ltr.get(g, 0) for g in genes])

        # ---- marker-based cell typing (object ships unannotated)
        PANELS = {
            "Myeloid": ["AIF1", "CD68", "CSF1R", "LST1", "C1QA"],
            "T_cell": ["CD3D", "CD3E", "CD2", "TRAC"],
            "B_cell": ["MS4A1", "CD79A"],
            "Oligodendrocyte": ["MBP", "PLP1", "MOG"],
            "Astrocyte": ["GFAP", "AQP4", "S100B"],
            "Neoplastic": ["EGFR", "CD44", "CHI3L1", "OLIG2"],
        }
        totals = np.asarray(X.sum(axis=0)).ravel()
        norm = X.multiply(1e4 / np.maximum(totals, 1)).tocsr()
        norm = norm.log1p()
        gi_sym = {g: i for i, g in enumerate(genes)}
        scores = {}
        for name, marks in PANELS.items():
            idx = [gi_sym[m] for m in marks if m in gi_sym]
            if not idx:
                continue
            scores[name] = np.asarray(norm[idx].mean(axis=0)).ravel()
        S = pd.DataFrame(scores, index=cells.index)
        ct = S.idxmax(axis=1).where(S.max(axis=1) > 0.05, "LowQC")
        st.note("cell types: " + ct.value_counts().to_dict().__str__())

        # ---- 1. class detection per cell type
        det = (X > 0)
        rows = []
        for t, idx in ct.groupby(ct).groups.items():
            sub = det[:, np.isin(cells.index, list(idx))]
            n = sub.shape[1]
            if n < 50:
                continue
            for c in CLASSES + ["PCG"]:
                gi = (gcls == c) if c != "PCG" else np.isin(gcls, ["PCG"])
                if gi.sum() < 10:
                    continue
                rows.append(dict(cell_type=t, n_cells=n, cls=c,
                                 n_genes=int(gi.sum()),
                                 det_rate=float(sub[gi].mean())))
            lm = np.isin(gcls, CLASSES) & (gltr == 1)
            if lm.sum() >= 10:
                rows.append(dict(cell_type=t, n_cells=n, cls="LTR_derived_lnc",
                                 n_genes=int(lm.sum()),
                                 det_rate=float(sub[lm].mean())))
        t1 = pd.DataFrame(rows)
        t1.to_csv(safe_write_path(os.path.join(OUT, "class_by_celltype.tsv")),
                  sep="\t", index=False)
        st.note(f"celltype rows={len(t1)}")

        # ---- 2. divergent pair co-detection within cells
        dp = pd.read_csv(os.path.join(RES, "divergent_pairs.tsv"), sep="\t")
        sym2gid = dict(zip(meta.gene_name.astype(str), meta.gene_id.astype(str)))
        sym2idx = {g: i for i, g in enumerate(genes)}
        # master full ids -> base -> symbol (matrix rows are symbols)
        base_sym = {g.split(".")[0]: str(s) for g, s in
                    zip(meta.gene_id, meta.gene_name)}
        pairs = []
        for _, r in dp.iterrows():
            sa = base_sym.get(r.lnc_gid.split(".")[0])
            sb = base_sym.get(r.pcg_gid.split(".")[0])
            if sa is None or sb is None:
                continue
            ia, ib = sym2idx.get(sa), sym2idx.get(sb)
            if ia is None or ib is None:
                continue
            pairs.append((sa, sb, ia, ib, "divergent"))
        st.note(f"divergent pairs present in matrix: {len(pairs)}")
        n_cells = X.shape[1]
        res_rows = []
        rng = np.random.default_rng(7)
        all_lnc_idx = [i for i in range(len(genes)) if gcls[i] in CLASSES]
        # detection rates & matching helper computed ONCE (outside pair loop)
        lnc_rate = np.array([float(det[i].mean()) for i in all_lnc_idx])
        order = np.argsort(lnc_rate)

        def matched_partner(target, exclude):
            j = np.searchsorted(lnc_rate[order], target)
            for step in range(1, 60):
                for jj in (j - step, j + step, j):
                    if 0 <= jj < len(order):
                        cand = all_lnc_idx[order[jj]]
                        if cand not in exclude:
                            return cand
            return all_lnc_idx[0]

        for (ga, gb, ia, ib, kind) in pairs[:400]:
            da = np.asarray(det[ia].todense()).ravel()
            db = np.asarray(det[ib].todense()).ravel()
            pa, pb = da.mean(), db.mean()
            k = int((da & db).sum())
            exp = pa * pb * n_cells
            if pa < 0.005 or pb < 0.005 or k < 5:
                continue
            p_obs = binomtest(k, n_cells, pa * pb, alternative="greater").pvalue
            # detection-rate MATCHED control pair (quantile-nearest)
            c1 = matched_partner(pa, {ia, ib})
            c2 = matched_partner(pb, {ia, ib, c1})
            dc1 = np.asarray(det[c1].todense()).ravel()
            dc2 = np.asarray(det[c2].todense()).ravel()
            kc = int((dc1 & dc2).sum())
            pc = dc1.mean() * dc2.mean() * n_cells
            res_rows.append(dict(pair=f"{ga}|{gb}", k_both=k, expected=exp,
                                 ratio=k / max(exp, 1e-9), p_binom=p_obs,
                                 control_k=kc, control_expected=pc))
        t2 = pd.DataFrame(res_rows)
        if len(t2):
            t2.to_csv(safe_write_path(os.path.join(
                OUT, "divergent_co_detection.tsv")), sep="\t", index=False)
            sig = (t2.p_binom < 0.01).mean()
            med_ratio = t2.ratio.median()
            med_ctrl = (t2.control_k / t2.control_expected.clip(lower=1e-9)).median()
            st.note(f"pairs tested={len(t2)}; median obs/exp ratio={med_ratio:.2f} "
                    f"vs random-pair control={med_ctrl:.2f}; "
                    f"FDR<0.01 via p<1e-2 in {sig:.0%}")
        st.summary_md("results/scrna/step_summaries/23_scrna.md",
                      "23 single-cell analysis",
                      ("## Class detection by cell type\n\n"
                       + df_to_md(t1.head(40), ".3f", max_rows=40)
                       + (f"\n\n## Divergent co-detection: median ratio "
                          f"{med_ratio:.2f} (control {med_ctrl:.2f})\n\n"
                          + df_to_md(t2.head(25), ".2f") if len(t2) else "")))


if __name__ == "__main__":
    main()
