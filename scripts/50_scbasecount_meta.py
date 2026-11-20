#!/usr/bin/env python
"""50_scbasecount_meta.py — screen the scBaseCount sample index before download.

Reads the human sample-level metadata that the atlas publishes next to the
matrices (3 MB, 35,263 samples) and turns it into a candidate list. Tumour and
normal samples are separated by the free-text disease field, grounded in the
MONDO identifier where one is present, and normals are matched to tumours by
tissue ontology so that a paired contrast is possible per cancer.

No matrices are touched here. Sizes are estimated at 2.1 kB per cell, measured
from a downloaded 3,562-cell sample (7.4 MB).

Run inside a container that has pyarrow, for example
  singularity exec --bind /media/share <scimilarity.sif> python \
      scripts/50_scbasecount_meta.py
"""
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pylib.progress import step, safe_write_path, df_to_md, WORKDIR

IN = os.path.join(WORKDIR, "data", "scbasecount_meta",
                  "human_sample_metadata.parquet")
OUT = os.path.join(WORKDIR, "results", "scbasecount")
KB_PER_CELL = 0.0021  # MB per cell

NORMAL_PAT = re.compile(
    r"\b(normal|healthy|control|adjacent|uninvolved|non[- ]?tumou?r)\b", re.I)
TUMOUR_PAT = re.compile(
    r"(cancer|tumou?r|neoplas|carcinoma|leuka?emia|lymphoma|myeloma|melanoma"
    r"|gliom|glioblastoma|sarcoma|adenoma|adenocarcinoma|squamous|mesothelioma"
    r"|neuroblastoma|medulloblastoma|cholangiocarcinoma|astrocytoma"
    r"|meningioma|thymoma|teratoma)", re.I)
# "none" alone is ambiguous, so it is kept as its own class rather than being
# counted as normal tissue.
UNSPEC_PAT = re.compile(
    r"^(none|not specified|unsure|other|unknown|none reported|)$", re.I)


def classify(txt):
    if txt is None or (isinstance(txt, float) and pd.isna(txt)):
        txt = ""
    t = str(txt).strip()
    if NORMAL_PAT.search(t):
        return "normal"
    if TUMOUR_PAT.search(t):
        return "tumour"
    if UNSPEC_PAT.match(t):
        return "unspecified"
    return "other"


def main():
    with step("50_scbasecount_meta") as st:
        os.makedirs(OUT, exist_ok=True)
        d = pd.read_parquet(IN)
        st.note("samples={} cells={:,}".format(
            len(d), int(d.obs_count.sum())))

        d["disease_class"] = d.disease.map(classify)
        st.note("classes: " + ", ".join(
            "{}={}".format(k, v) for k, v in
            d.disease_class.value_counts().items()))

        # Group tumours by ontology identifier where present, else by cleaned
        # text. The modal text supplies a readable label for each group.
        d["grp"] = d.disease_ontology_term_id.fillna("").astype(str)
        blank = d.grp == ""
        d.loc[blank, "grp"] = (d.loc[blank, "disease"].fillna("")
                               .str.lower().str.strip().str[:40])
        lab = (d[d.disease_class == "tumour"]
               .groupby("grp").disease
               .agg(lambda s: s.value_counts().index[0]))
        d["cancer_group"] = d.grp.map(lab).fillna("")

        tum = d[d.disease_class == "tumour"].copy()
        nor = d[d.disease_class == "normal"].copy()

        rows = []
        for g, sub in tum.groupby("grp"):
            ids = sub.tissue_ontology_term_id.dropna().astype(str)
            ids = ids[ids.str.strip() != ""]
            tissue = ids.value_counts().index[0] if len(ids) else ""
            # A blank tissue identifier must never match the normal pool, or
            # every unannotated tumour borrows the same normal samples.
            if tissue.strip() == "":
                same = nor.iloc[0:0]
            else:
                same = nor[nor.tissue_ontology_term_id.astype(str) == tissue]
            rows.append(dict(
                cancer_group=sub.disease.value_counts().index[0][:60],
                n_tumour=len(sub),
                cells_tumour=int(sub.obs_count.sum()),
                tissue_ontology=tissue,
                tissue_label=str(sub.tissue.value_counts().index[0])[:45],
                n_normal_tissue=len(same),
                cells_normal_tissue=int(same.obs_count.sum()),
                n_tumour_high_conf=int((sub.single_disease_confidence ==
                                        "high").sum())))
        T = pd.DataFrame(rows)
        T["est_GB_two_arms"] = ((T.cells_tumour + T.cells_normal_tissue)
                                * KB_PER_CELL / 1000).round(2)
        T = T.sort_values("n_tumour", ascending=False)
        T.to_csv(safe_write_path(os.path.join(OUT, "cancer_summary.tsv")),
                 sep="\t", index=False)

        cand = d[d.disease_class.isin(["tumour", "normal"])].copy()
        cand["est_MB"] = (cand.obs_count * KB_PER_CELL).round(1)
        keep = ["srx_accession", "file_path", "organism", "tissue",
                "tissue_ontology_term_id", "disease",
                "disease_ontology_term_id", "disease_class", "cancer_group",
                "obs_count", "est_MB", "tech_10x", "cell_prep",
                "single_disease_confidence"]
        cand[keep].to_csv(
            safe_write_path(os.path.join(OUT, "candidate_samples.tsv")),
            sep="\t", index=False)

        paired = T[(T.n_tumour >= 10) & (T.n_normal_tissue >= 10)]
        st.note("tumour groups={}; with >=10 samples each arm={}".format(
            len(T), len(paired)))
        st.note("candidate samples={} (tumour {} / normal {}), cells={:,}, "
                "est download={:.1f} GB".format(
                    len(cand), int((cand.disease_class == "tumour").sum()),
                    int((cand.disease_class == "normal").sum()),
                    int(cand.obs_count.sum()),
                    cand.obs_count.sum() * KB_PER_CELL / 1000))
        st.summary_md(os.path.join(OUT, "step_summaries", "50.md"),
                      "50 scBaseCount sample screen",
                      "Human samples in the 2026-01-12 release. Size assumes "
                      "2.1 kB per cell.\n\n"
                      + df_to_md(T.head(30), ".3f", max_rows=30))

        # ---- a balanced, download-sized design -----------------------------
        # Solid tumours that carry genuine same-tissue normals and are large
        # enough to cap at PER_ARM samples per arm. Normals are de-duplicated
        # across groups because several breast groups share one control pool.
        design = design_samples(d, T)
        st.note("design: {} samples, {:,} cells, est {:.1f} GB".format(
            len(design), int(design.obs_count.sum()),
            design.est_MB.sum() / 1000))


DESIGN_GROUPS = ["glioblastoma", "breast cancer", "hepatocellular carcinoma",
                 "lung adenocarcinoma", "non-small cell lung cancer (NSCLC)",
                 "clear cell renal cell carcinoma (ccRCC)",
                 "colorectal cancer"]
PER_ARM = 40


def design_samples(d, T):
    rows = []
    used_normals = set()
    for grp in DESIGN_GROUPS:
        sub = T[T.cancer_group == grp]
        if not len(sub):
            continue
        tissue = str(sub.tissue_ontology.iloc[0])
        tum = d[(d.disease_class == "tumour") & (d.cancer_group == grp)]
        tum = tum.sort_values("obs_count", ascending=False).head(PER_ARM)
        nor = d[(d.disease_class == "normal") &
                (d.tissue_ontology_term_id.astype(str) == tissue) &
                (~d.srx_accession.isin(used_normals))]
        nor = nor.sort_values("obs_count", ascending=False).head(PER_ARM)
        used_normals.update(nor.srx_accession)
        for _, r in pd.concat([tum, nor]).iterrows():
            rows.append(dict(group=grp, arm=r.disease_class,
                             srx=r.srx_accession, file_path=r.file_path,
                             tissue=str(r.tissue)[:40],
                             obs_count=int(r.obs_count),
                             est_MB=round(r.obs_count * KB_PER_CELL, 1)))
    out = pd.DataFrame(rows)
    out.to_csv(safe_write_path(os.path.join(OUT, "design_samples.tsv")),
               sep="\t", index=False)
    return out


if __name__ == "__main__":
    main()
