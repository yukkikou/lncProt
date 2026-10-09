#!/usr/bin/env Rscript
# 21_scrna_export.R — open CGGA IDH_A Seurat object, export counts (MatrixMarket)
# + cell metadata for downstream Python analysis. Runs inside seurat.sif.
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
W <- "/media/share/workdir/1610305236/lncProt"
OUT <- file.path(W, "data", "ccle_cgga", "cgga", "CGGA_CCell_4083", "scrna_export")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

e <- new.env()
load(file.path(W, "data", "ccle_cgga", "cgga", "CGGA_CCell_4083",
               "CGGA_IDH_A_scRNA_SeuratObj_20250915.rdata"), envir = e)
objs <- ls(envir = e)
message("loaded objects: ", paste(objs, collapse = ", "))
obj <- get(objs[1], envir = e)
if (!inherits(obj, "Seurat")) {
  for (o in objs) {
    if (inherits(get(o, envir = e), "Seurat")) { obj <- get(o, envir = e); break }
  }
}
message("class: ", class(obj)[1], " | dims: ", nrow(obj), " genes x ", ncol(obj), " cells")

counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
meta <- obj@meta.data
message("meta cols: ", paste(colnames(meta), collapse = ", "))

writeMM(counts, file.path(OUT, "counts.mtx"))
writeLines(rownames(counts), file.path(OUT, "genes.txt"))
write.table(meta, file.path(OUT, "cell_meta.tsv"), sep = "\t",
            quote = FALSE, row.names = TRUE)
um <- tryCatch(Embeddings(obj, reduction = "umap"), error = function(e) NULL)
if (!is.null(um)) {
  write.table(um, file.path(OUT, "umap.tsv"), sep = "\t", quote = FALSE,
              row.names = TRUE)
}
message("export complete: ", ncol(counts), " cells, ", nrow(counts), " genes")
