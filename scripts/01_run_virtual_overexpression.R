#!/usr/bin/env Rscript
## ============================================================================
##  scripts/01_run_virtual_overexpression.R
##
##  In-silico 2-fold gain-of-function perturbation of a target gene in the
##  macrophages of the control group, followed by differential regulation.
##
##  Usage
##    Rscript scripts/01_run_virtual_overexpression.R
##  Edit the settings block below; all outputs go to results/.
##
##  The analysis depends only on the public scTenifoldNet package. See
##  R/pipeline.R for what each step does and README.md for the method notes.
## ============================================================================

## ---- settings --------------------------------------------------------------

DATA_FILE       <- "seurat.Rdata"   # Seurat object with a counts assay
OBJECT_NAME     <- NULL             # name of the Seurat object; NULL = auto

TARGET_GENE     <- "HLA-DOA"
CELL_TYPE       <- "macrophages"
CELLTYPE_COLUMN <- "cellType"       # metadata column holding the cell type
GROUP_COLUMN    <- "Type"           # metadata column holding the group
GROUP_KEEP      <- "Control"        # the perturbation starts from the WT group

N_FEATURES      <- 10000            # highly variable genes (vst)
MIN_LIBRARY     <- 1000             # cells with fewer counts are dropped
MT_THRESHOLD    <- 0.1              # maximum mitochondrial fraction
MIN_CELLS       <- 0                # genes detected in fewer cells are dropped

N_NET           <- 10               # networks built from cell subsamples
N_CELLS         <- 500              # cells per network
N_COMP          <- 3                # principal components per regression
Q_EDGES         <- 0.9              # keep the strongest 10% of edges
TD_K            <- 3                # rank of the CP decomposition
TD_DECIMALS     <- 3
NC_LAMBDA       <- 0                # directional pruning (0 = off)
MA_DIM          <- 2                # dimensions of the aligned embedding

FOLD            <- 2                # 2-fold gain-of-function
FDR_CUTOFF      <- 0.05
SEED            <- 123
N_CORES         <- parallel::detectCores()

OUT_DIR         <- "results"
PREVIOUS_RESULT <- NULL             # optional sigDiff.txt from an earlier run

## ---- setup -----------------------------------------------------------------

required <- c("Seurat", "scTenifoldNet", "Matrix", "MASS", "ggplot2", "ggrepel")
missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("missing packages: ", paste(missing, collapse = ", "))
}
suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir  <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()
root_dir    <- normalizePath(file.path(script_dir, ".."))

source(file.path(root_dir, "R", "pipeline.R"))
source(file.path(root_dir, "R", "plots.R"))

setwd(root_dir)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- input -----------------------------------------------------------------

if (!file.exists(DATA_FILE)) {
  stop("cannot find ", DATA_FILE, " in ", getwd(),
       "; set DATA_FILE to its full path")
}

env <- new.env()
load(DATA_FILE, envir = env)
objects <- ls(env)
seurat_objects <- objects[vapply(objects,
                                 function(x) inherits(get(x, envir = env), "Seurat"),
                                 logical(1))]

if (is.null(OBJECT_NAME)) {
  if (length(seurat_objects) != 1) {
    stop("found ", length(seurat_objects), " Seurat objects (",
         paste(seurat_objects, collapse = ", "), "); set OBJECT_NAME")
  }
  OBJECT_NAME <- seurat_objects
}
obj <- get(OBJECT_NAME, envir = env)
rm(env); invisible(gc())

meta <- obj[[]]
for (col in c(GROUP_COLUMN, CELLTYPE_COLUMN)) {
  if (!col %in% colnames(meta)) {
    stop("metadata column '", col, "' not found; available: ",
         paste(colnames(meta), collapse = ", "))
  }
}
cells_keep <- rownames(meta)[meta[[GROUP_COLUMN]] == GROUP_KEEP &
                             meta[[CELLTYPE_COLUMN]] == CELL_TYPE]
if (length(cells_keep) == 0) {
  stop("no cells with ", GROUP_COLUMN, " == '", GROUP_KEEP, "' and ",
       CELLTYPE_COLUMN, " == '", CELL_TYPE, "'")
}
obj <- subset(obj, cells = cells_keep)

obj <- FindVariableFeatures(obj, selection.method = "vst",
                           nfeatures = N_FEATURES)
features <- unique(c(TARGET_GENE, VariableFeatures(obj)))

counts <- tryCatch(GetAssayData(obj, layer = "counts"),
                   error = function(e) GetAssayData(obj, slot = "counts"))

if (!TARGET_GENE %in% rownames(counts)) {
  stop("target gene ", TARGET_GENE, " is not in the count matrix")
}
counts <- counts[features, , drop = FALSE]

rm(obj); invisible(gc())

message(sprintf("input: %d genes x %d cells (%s, %s)",
                nrow(counts), ncol(counts), GROUP_KEEP, CELL_TYPE))

## ---- quality control -------------------------------------------------------

counts <- qc_cells(counts, mt_threshold = MT_THRESHOLD,
                   min_library = MIN_LIBRARY, min_cells = MIN_CELLS)

## ---- perturbation ----------------------------------------------------------

result <- run_perturbation(counts,
                           target_gene   = TARGET_GENE,
                           mode          = "gain",
                           fold          = FOLD,
                           n_net         = N_NET,
                           n_cells       = N_CELLS,
                           n_comp        = N_COMP,
                           q_edges       = Q_EDGES,
                           td_k          = TD_K,
                           td_decimals   = TD_DECIMALS,
                           nc_lambda     = NC_LAMBDA,
                           manifold_dim  = MA_DIM,
                           n_cores       = N_CORES,
                           seed          = SEED)

## ---- output ----------------------------------------------------------------

diff <- result$diff
diff$is_target <- diff$gene == TARGET_GENE

write.table(diff,
            file = file.path(OUT_DIR, "diffRegulation_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sig <- diff[!diff$is_target & diff$p.adj < FDR_CUTOFF, , drop = FALSE]
write.table(sig,
            file = file.path(OUT_DIR, "diffRegulation_FDR05.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

save_figure(plot_top_genes(diff, top_n = 20),
            file.path(OUT_DIR, "fig_barplot_top20"), 6, 5)
save_figure(plot_volcano(diff, fdr = FDR_CUTOFF),
            file.path(OUT_DIR, "fig_volcano"), 6, 5)
save_figure(plot_proportion(diff, fdr = FDR_CUTOFF),
            file.path(OUT_DIR, "fig_proportion"), 6, 5)

## ---- run record ------------------------------------------------------------

seurat_ver <- tryCatch(as.character(packageVersion("Seurat")),
                       error = function(e) "not installed")

params <- c(
  sprintf("run_date\t%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("R\t%s", R.version.string),
  sprintf("scTenifoldNet\t%s", as.character(packageVersion("scTenifoldNet"))),
  sprintf("Seurat\t%s", seurat_ver),
  sprintf("data_file\t%s", DATA_FILE),
  sprintf("target_gene\t%s", TARGET_GENE),
  sprintf("group\t%s", GROUP_KEEP),
  sprintf("cell_type\t%s", CELL_TYPE),
  sprintf("n_features\t%d", N_FEATURES),
  sprintf("qc_min_library\t%d", MIN_LIBRARY),
  sprintf("qc_mt_threshold\t%s", MT_THRESHOLD),
  sprintf("qc_min_cells\t%d", MIN_CELLS),
  sprintf("n_net\t%d", N_NET),
  sprintf("n_cells\t%d", N_CELLS),
  sprintf("n_comp\t%d", N_COMP),
  sprintf("q_edges\t%s", Q_EDGES),
  sprintf("td_k\t%d", TD_K),
  sprintf("td_decimals\t%d", TD_DECIMALS),
  sprintf("nc_lambda\t%s", NC_LAMBDA),
  sprintf("manifold_dim\t%d", MA_DIM),
  sprintf("fold\t%s", FOLD),
  sprintf("seed\t%d", SEED),
  sprintf("n_cores\t%d", N_CORES),
  sprintf("genes_in_network\t%d", nrow(counts)),
  sprintf("cells_in_network\t%d", ncol(counts)),
  sprintf("significant_genes\t%d", nrow(sig))
)
writeLines(params, file.path(OUT_DIR, "run_parameters.txt"))
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))

## ---- sanity check against an earlier run -----------------------------------

if (!is.null(PREVIOUS_RESULT) && file.exists(PREVIOUS_RESULT)) {
  old <- utils::read.delim(PREVIOUS_RESULT)
  old_genes <- unique(old$gene)
  new_genes <- unique(sig$gene)
  overlap <- length(intersect(old_genes, new_genes))
  union_n <- length(union(old_genes, new_genes))
  message(sprintf("overlap with %s: %d shared genes (Jaccard %.3f)",
                  PREVIOUS_RESULT, overlap, overlap / union_n))
}

## ---- summary ---------------------------------------------------------------

message(sprintf("%d of %d genes below FDR %s (perturbed gene excluded; it is
kept in diffRegulation_full.tsv)", nrow(sig), nrow(diff), FDR_CUTOFF))
message("top genes by network perturbation score:")
print(utils::head(sig[order(-sig$FC), c("gene", "distance", "Z", "FC", "p.adj")], 10))
message("outputs written to ", normalizePath(OUT_DIR))
