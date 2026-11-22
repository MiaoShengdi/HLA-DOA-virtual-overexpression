## Smoke test on a small synthetic matrix.
## Checks that the pipeline runs end to end, that the output has the expected
## shape, and that two runs with the same seed give identical results.
##
## Usage: Rscript tests/test_pipeline_synthetic.R

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(normalizePath(script_file)), ".."))

source(file.path(root, "R", "pipeline.R"))

ok <- function(label, condition) {
  if (!isTRUE(condition)) stop("FAILED: ", label)
  cat("ok  ", label, "\n")
}

## ---- synthetic counts ------------------------------------------------------

set.seed(1)
n_genes <- 200
n_cells <- 600
counts <- matrix(
  rnbinom(n_genes * n_cells, mu = rep(rgamma(n_genes, 2, 1), n_cells), size = 3),
  nrow = n_genes
)
rownames(counts) <- c(paste0("MT-", 1:5), paste0("G", 1:195))
colnames(counts) <- paste0("C", 1:n_cells)

## ---- quality control -------------------------------------------------------

qc <- qc_cells(counts, mt_threshold = 0.1, min_library = 150, verbose = FALSE)
ok("QC returns a sparse matrix", inherits(qc, "dgCMatrix"))
ok("QC keeps genes and cells", nrow(qc) == nrow(counts) && ncol(qc) > 0)

## ---- perturbation helper ---------------------------------------------------

wt <- matrix(rnorm(16), 4, 4, dimnames = list(paste0("G", 1:4), paste0("G", 1:4)))
gain <- perturb_network(wt, "G2", mode = "gain", fold = 2)
loss <- perturb_network(wt, "G2", mode = "loss")

ok("gain scales the outgoing edges by 2",
   isTRUE(all.equal(gain["G2", ], 2 * wt["G2", ])))
ok("gain leaves other genes untouched", isTRUE(all.equal(gain["G1", ], wt["G1", ])))
ok("loss sets the outgoing edges to zero", all(loss["G2", ] == 0))

## ---- pipeline --------------------------------------------------------------

run <- function() {
  run_perturbation(counts = qc, target_gene = "G10", mode = "gain", fold = 2,
                   n_net = 3, n_cells = 200, n_comp = 3, q_edges = 0.9,
                   td_k = 3, td_decimals = 2, manifold_dim = 2, n_cores = 1,
                   seed = 11, verbose = FALSE)
}

res1 <- run()
res2 <- run()

ok("network is square and labelled",
   all(dim(res1$wt) == nrow(qc)) && identical(rownames(res1$wt), colnames(res1$wt)))
ok("perturbed network differs only in the target row",
   isTRUE(all.equal(res1$perturbed[-which(rownames(res1$wt) == "G10"), ],
                    res1$wt[-which(rownames(res1$wt) == "G10"), ])))
ok("diff table covers every gene", nrow(res1$diff) == nrow(qc))
ok("diff table has the expected columns",
   identical(names(res1$diff),
             c("gene", "distance", "Z", "FC", "p.value", "p.adj")))
ok("statistics are finite and in range",
   all(is.finite(res1$diff$FC)) &&
   all(res1$diff$p.adj >= 0 & res1$diff$p.adj <= 1))
ok("target gene is the most displaced gene", res1$diff$gene[1] == "G10")
ok("two runs with the same seed agree", isTRUE(all.equal(res1$diff, res2$diff)))

cat("\nall tests passed\n")
