# In-silico 2-fold gain-of-function perturbation of a single gene

Reproducible analysis code for a virtual over-expression experiment on a
single-cell gene regulatory network (scGRN). The target gene of the study is
`HLA-DOA` in control macrophages; any gene and any fold change can be set in
`scripts/01_run_virtual_overexpression.R`.

The perturbation is applied to the **inferred regulatory network**, not to the
expression matrix: after the WT network is built and transposed so that rows
hold the outgoing (regulatory) edges of a gene, the outgoing edge weights of the
target gene are scaled by the fold factor. The count matrix is left untouched.

    loss-of-function (knockout)   A[g, ] <- 0
    gain-of-function  (fold = 2)  A[g, ] <- 2 * A[g, ]

## Pipeline

1. Cells of the WT group are filtered (library size, mitochondrial fraction,
   prediction bands of detected genes and mitochondrial reads).
2. Highly variable genes are selected (`vst`), the target gene is added.
3. `nNet` networks are built from `nCells`-cell subsamples by principal
   component regression (top 10% of edges kept).
4. The networks are combined by a rank-`K` CP tensor decomposition, which
   denoises and averages the edge weights.
5. The consensus network is transposed and the target gene is perturbed.
6. WT and perturbed networks are embedded by non-linear manifold alignment.
7. Differentially regulated genes are ranked by the Euclidean displacement of
   each gene between the two embeddings: Box-Cox normalisation, chi-square test
   on the squared-distance ratio, Benjamini-Hochberg correction.

Steps 3, 4 and 6 are calls to the public **scTenifoldNet** package
(Osorio et al., *Patterns* 2020). Steps 1, 5 and 7 are implemented in
`R/pipeline.R`.

## Why the code is self-contained

The repository does not include, redistribute or depend on any third-party
modified build of `scTenifoldKnk`. Only the public `scTenifoldNet` package is
used; the cell QC, the perturbation and the differential-regulation statistic are
implemented here, from the published method, so that the analysis can be audited
line by line. The perturbation is one assignment (`perturb_network()`), which
makes the difference between a knockout and a gain-of-function run explicit.

The statistic reproduces the published formulation exactly: verified column by
column against the reference implementation (`distance`, `Z`, `FC`, `p.value`,
`p.adj` all differ by 0), and the cell filter selects the same cells. The
statistic is a closed-form function of the two embeddings, so it does not depend
on which scTenifoldNet version supplied them.

## Requirements

- R >= 4.0
- scTenifoldNet **1.3** (the version used for the reported results; the same
  script also runs on 1.4.3, see below), Seurat, Matrix, MASS, ggplot2, ggrepel,
  RSpectra, RhpcBLASctl

```r
# install scTenifoldNet 1.3 from a source tarball
install.packages("scTenifoldNet_1.3.tar.gz", repos = NULL, type = "source")
```

### Package version

All reported numbers were produced with **scTenifoldNet 1.3**. Network
construction was rewritten in 1.4.3 (the principal-component regression moved to
a different eigensolver), so the inferred network — and therefore the ranking of
perturbed genes — differs between the two versions. Results obtained under 1.3
and 1.4.3 must not be mixed in one table or figure.

The version used for a run is written automatically to
`results/run_parameters.txt`, and the full session to `results/sessionInfo.txt`,
so the reported version can be checked against the actual output.

The three calls to the package are version independent: `makeNetworks()`,
`tensorDecomposition()` and `manifoldAlignment()` accept the same arguments in
1.3 and 1.4.3, so the script runs unchanged under either version and only the
numbers change.

## Usage

Edit the settings block at the top of `scripts/01_run_virtual_overexpression.R`
(input file, target gene, cell type, fold, number of genes) and run from the
repository root:

```sh
Rscript scripts/01_run_virtual_overexpression.R
```

Smoke test on a synthetic matrix (no data needed):

```sh
Rscript tests/test_pipeline_synthetic.R
```

## Settings used for the reported analysis

| Setting | Value |
| --- | --- |
| Cells | control-group macrophages |
| Genes | 10,000 highly variable genes (vst) plus the target gene |
| Cell filter | library size >= 1,000 counts, mitochondrial fraction <= 0.1 |
| Gene filter | none (`MIN_CELLS = 0`) |
| Networks | 10, from 500-cell subsamples |
| Regression | 3 principal components, `q = 0.9` (top 10% of edges kept) |
| Tensor | `K = 3`, max error 1e-5, 1,000 iterations, 3 decimals |
| Directional pruning | off (`lambda = 0`) |
| Embedding | 2 dimensions |
| Perturbation | outgoing edges of HLA-DOA multiplied by 2 |
| Test | chi-square on squared displacement, Benjamini-Hochberg FDR < 0.05 |
| Seed | 123 |
| scTenifoldNet | 1.3 |

The gene filter is deliberately set to zero: the reference implementation
computes a per-gene detection filter but never assigns the result, so the
analysis it performs applies no gene filter. Writing `MIN_CELLS = 0` records
that behaviour instead of hiding it behind a default.

## Outputs (`results/`)

| File | Content |
| --- | --- |
| `diffRegulation_full.tsv` | every gene in the network, unfiltered, with an `is_target` flag |
| `diffRegulation_FDR05.tsv` | genes below the FDR cutoff, perturbed gene excluded |
| `fig_barplot_top20.png/.pdf` | top genes by network perturbation score |
| `fig_volcano.png/.pdf` | displacement (Z) against -log10(FDR) |
| `fig_proportion.png/.pdf` | share of significantly displaced genes |
| `run_parameters.txt` | all settings, seed, core count, package versions |
| `sessionInfo.txt` | full R session record |

## Interpretation notes

- `FC` is a **network perturbation score** (squared-displacement ratio against
  the expected displacement), not an expression fold change. It must not be
  reported as a fold change of mRNA levels.
- The method ranks genes by how far they move in the network embedding. It does
  **not** assign a direction, so up- or down-regulation of a target should not be
  claimed from these results alone.
- The perturbed gene is the most displaced gene by construction and is therefore
  excluded from the reported hit list; it stays in the full table.
- With `nGenes = 10,000` and `nNet = 10` the tensor alone needs about 7.5 GB
  (`nGenes^2 * nNet * 8` bytes). Reduce the number of genes or run on a machine
  with enough memory; the script prints the estimate before it starts.

## Reproducibility

`SEED` fixes network subsampling and the embedding. Results are reproducible on
the same scTenifoldNet version with the same number of cores; both are recorded
in `results/run_parameters.txt`. Two consecutive runs of the smoke test give
identical output tables.

## References

- Osorio D, Zhong Y, Li G, Xu Q, Yang Y, Tian Y, Chapkin RS, Huang J, Cai JJ.
  scTenifoldNet: a machine learning workflow for constructing and comparing
  single-cell gene regulatory networks. *Patterns* 2020;1:100139.
- scTenifoldKnk, the virtual-knockout framework this perturbation follows:
  https://github.com/cailab-tamu/scTenifoldKnk

## Licence

GPL-3.0. See `LICENSE`.

## Contact

Miao Shengdi, miaoshengdi@hunnu.edu.cn
