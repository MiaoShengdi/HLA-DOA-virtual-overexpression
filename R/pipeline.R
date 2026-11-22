## ============================================================================
##  R/pipeline.R
##
##  In-silico perturbation of one gene in a single-cell gene regulatory network
##  (scGRN), built on the public scTenifoldNet package.
##
##  Network construction, CP tensor decomposition and manifold alignment are
##  delegated to scTenifoldNet (Osorio et al., Patterns 2020, 1:100139).
##  Cell-level QC, the perturbation of the adjacency matrix and the
##  differential-regulation statistic are implemented here, so the whole
##  analysis runs from this repository without any third-party modified build.
##
##  Perturbation semantics (after the network is transposed so that rows hold
##  the outgoing, i.e. regulatory, edges of a gene):
##
##      loss-of-function / knockout   A[g, ] <- 0
##      gain-of-function  / fold = 2  A[g, ] <- 2 * A[g, ]
##
##  The expression matrix itself is never modified: what is scaled is the
##  inferred regulatory output of the gene, not its mRNA count.
##
##  Method notes
##    * The displacement statistic follows the virtual-knockout framework:
##      Euclidean distance between the WT and perturbed embeddings of the same
##      gene, Box-Cox normalised, chi-square test on the squared-distance ratio,
##      Benjamini-Hochberg correction.
##    * FC is a network perturbation score, not an expression fold change.
## ============================================================================


## ---- cell-level quality control --------------------------------------------

## Prediction band of y given x, used to discard cells whose number of detected
## genes (or mitochondrial reads) is outside the expected range for their
## library size.
.prediction_band <- function(y, x) {
  fit <- stats::lm(y ~ x)
  as.data.frame(stats::predict(fit, data.frame(x = x), interval = "prediction"))
}

#' Filter cells of a count matrix
#'
#' @param counts genes x cells, matrix, data.frame or dgCMatrix.
#' @param mt_threshold maximum fraction of mitochondrial reads.
#' @param min_library minimum library size (total counts per cell).
#' @param min_cells minimum number of cells in which a gene must be detected
#'   to be kept. 0 keeps every gene, which is what the original analysis did;
#'   genes that are constant across cells receive no edges downstream anyway.
#' @return filtered count matrix (dgCMatrix), genes in rows, cells in columns.
qc_cells <- function(counts,
                     mt_threshold = 0.1,
                     min_library  = 1000,
                     min_cells    = 0,
                     verbose      = TRUE) {

  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!inherits(counts, "dgCMatrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)
  stopifnot(!is.null(rownames(counts)), !is.null(colnames(counts)))

  n_cells_in <- ncol(counts)

  lib <- Matrix::colSums(counts)
  counts <- counts[, lib >= min_library, drop = FALSE]
  if (ncol(counts) == 0) {
    stop("no cell has a library size >= ", min_library,
         " (largest observed: ", max(lib), "); check the units of the input")
  }
  lib <- Matrix::colSums(counts)

  n_expr <- Matrix::colSums(counts != 0)
  band_expr <- .prediction_band(n_expr, lib)

  mt_idx <- grep("^MT-", toupper(rownames(counts)))

  if (length(mt_idx) > 0) {
    mt_reads <- Matrix::colSums(counts[mt_idx, , drop = FALSE])
    band_mt  <- .prediction_band(mt_reads, lib)
    keep <- (n_expr > band_expr$lwr) & (n_expr < band_expr$upr) &
            (mt_reads > band_mt$lwr) & (mt_reads < band_mt$upr) &
            (mt_reads / lib <= mt_threshold) &
            (lib < 2 * mean(lib))
  } else {
    keep <- (n_expr > band_expr$lwr) & (n_expr < band_expr$upr) &
            (lib < 2 * mean(lib))
  }

  if (min_cells > 0) {
    counts <- counts[Matrix::rowSums(counts != 0) >= min_cells, , drop = FALSE]
  }

  counts <- counts[, keep, drop = FALSE]

  if (verbose) {
    message(sprintf("QC: %d -> %d cells, %d genes kept",
                    n_cells_in, ncol(counts), nrow(counts)))
  }
  counts
}


## ---- adjacency matrix helpers ----------------------------------------------

#' Remove the weaker of the two reciprocal edges, keeping one direction per pair
#'
#' lambda = 0 leaves the matrix unchanged.
strict_direction <- function(adj, lambda = 0) {
  adj <- as.matrix(adj)
  stronger <- adj
  stronger[abs(stronger) < abs(t(stronger))] <- 0
  Matrix::Matrix((1 - lambda) * adj + lambda * stronger)
}

#' Perturb the outgoing edges of a gene
#'
#' @param adj genes x genes adjacency matrix, rows = regulators.
#' @param gene gene whose outgoing edges are perturbed (must be a row name).
#' @param mode "gain" scales the edges by `fold`, "loss" sets them to zero.
#' @param fold scaling factor for the gain-of-function mode.
perturb_network <- function(adj, gene, mode = c("gain", "loss"), fold = 2) {

  mode <- match.arg(mode)
  adj <- as.matrix(adj)
  if (!gene %in% rownames(adj)) {
    stop("gene ", gene, " is not in the network")
  }

  out <- adj
  if (mode == "loss") {
    out[gene, ] <- 0
  } else {
    out[gene, ] <- out[gene, ] * fold
  }
  out
}


## ---- differential regulation -----------------------------------------------

#' Displacement of each gene between the WT and the perturbed embedding
#'
#' @param aligned output of scTenifoldNet::manifoldAlignment.
#' @param target_gene the perturbed gene.
#' @param exclude_top whether the expected displacement is computed over all
#'   genes except the most displaced one. TRUE reproduces the reference
#'   implementation, in which the perturbed gene is always the most displaced.
#' @return data.frame with gene, distance, Z, FC, p.value, p.adj.
diff_regulation <- function(aligned, target_gene, exclude_top = TRUE) {

  x_rows <- grep("^X_", rownames(aligned))
  y_rows <- grep("^Y_", rownames(aligned))
  genes  <- sub("^X_", "", rownames(aligned)[x_rows])

  if (length(x_rows) != length(y_rows) ||
      !identical(genes, sub("^Y_", "", rownames(aligned)[y_rows]))) {
    stop("X_ and Y_ blocks are not aligned; check manifoldAlignment output")
  }

  d <- sqrt(rowSums((aligned[x_rows, , drop = FALSE] -
                     aligned[y_rows, , drop = FALSE])^2))

  ## order as in the reference implementation: largest displacement first
  ord   <- order(d, decreasing = TRUE)
  d     <- unname(d[ord])
  genes <- genes[ord]

  ## Box-Cox transformation of the displacements
  lambda_grid <- seq(-2, 2, length.out = 1000)
  lambda_grid <- lambda_grid[lambda_grid != 0]
  bc <- try(MASS::boxcox(d[d > 0] ~ 1, plot = FALSE, lambda = lambda_grid),
            silent = TRUE)
  if (inherits(bc, "try-error")) {
    dn <- d
  } else {
    lambda_hat <- bc$x[which.max(bc$y)]
    dn <- if (lambda_hat < 0) 1 / (d^lambda_hat) else d^lambda_hat
  }

  z <- if (stats::sd(dn) > 0) as.numeric(scale(dn)) else rep(0, length(dn))

  ## squared-displacement ratio against the expected displacement
  fc  <- d^2 / mean(d[-seq_len(if (exclude_top) 1L else 0L)]^2)
  p   <- stats::pchisq(fc, df = 1, lower.tail = FALSE)
  padj <- stats::p.adjust(p, method = "fdr")

  if (exclude_top && genes[1] != target_gene) {
    warning("the most displaced gene is ", genes[1],
            ", not ", target_gene,
            "; the expected displacement excludes that gene, as in the reference")
  }

  data.frame(gene = genes, distance = d, Z = z, FC = fc,
             p.value = p, p.adj = padj, stringsAsFactors = FALSE)
}


## ---- full pipeline ---------------------------------------------------------

#' Run one in-silico perturbation experiment
#'
#' @param counts filtered genes x cells count matrix.
#' @param target_gene gene to perturb.
#' @param mode "gain" (over-expression) or "loss" (knockout).
#' @param fold scaling factor for "gain".
#' @param n_net number of networks built from cell subsamples.
#' @param n_cells cells sampled per network.
#' @param n_comp principal components used in each regression.
#' @param q_edges quantile at which weak edges are dropped (0.9 keeps the
#'   strongest 10%).
#' @param td_k rank of the CP decomposition.
#' @param td_decimals rounding of the reconstructed network.
#' @param nc_lambda strength of the directional pruning (0 = none).
#' @param manifold_dim dimensions of the aligned embedding.
#' @param n_cores cores handed to scTenifoldNet.
#' @param seed seed for network construction and alignment.
#' @return list with the WT network, the perturbed network, the aligned
#'   embedding and the differential-regulation table.
run_perturbation <- function(counts,
                             target_gene,
                             mode = c("gain", "loss"),
                             fold = 2,
                             n_net = 10,
                             n_cells = 500,
                             n_comp = 3,
                             q_edges = 0.9,
                             td_k = 3,
                             td_max_iter = 1000,
                             td_max_error = 1e-5,
                             td_decimals = 3,
                             nc_lambda = 0,
                             manifold_dim = 2,
                             n_cores = 1,
                             seed = 123,
                             verbose = TRUE) {

  mode <- match.arg(mode)
  if (!target_gene %in% rownames(counts)) {
    stop("target gene ", target_gene, " is not in the count matrix")
  }

  if (verbose) {
    est_gb <- nrow(counts)^2 * n_net * 8 / 1024^3
    message(sprintf("genes: %d | networks: %d | estimated tensor: %.1f GB",
                    nrow(counts), n_net, est_gb))
  }

  set.seed(seed)

  nets <- scTenifoldNet::makeNetworks(X = counts, q = q_edges, nNet = n_net,
                                      nCells = n_cells, nComp = n_comp,
                                      scaleScores = TRUE, symmetric = FALSE,
                                      nCores = n_cores)

  tens <- scTenifoldNet::tensorDecomposition(xList = nets, K = td_k,
                                             maxError = td_max_error,
                                             maxIter = td_max_iter,
                                             nDecimal = td_decimals)

  wt <- strict_direction(as.matrix(tens$X), lambda = nc_lambda)
  wt <- as.matrix(wt)
  diag(wt) <- 0
  wt <- t(wt)                       # rows = regulators, columns = targets

  pert <- perturb_network(wt, target_gene, mode = mode, fold = fold)

  aligned <- scTenifoldNet::manifoldAlignment(wt, pert, d = manifold_dim,
                                              nCores = n_cores)

  dr <- diff_regulation(aligned, target_gene)

  list(wt = wt, perturbed = pert, aligned = aligned, diff = dr)
}
