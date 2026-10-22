# BayesMAPJoint

**Joint-indicator Bayesian multi-level association mapping**

`BayesMAPJoint` is an experimental R package for hierarchical Bayesian
association mapping across three connected biological levels:

- SNPs (`delta`)
- genes (`Delta`)
- cell types (`gamma`)

The package supports individual-level genotype and phenotype data, dense GWAS
summary-statistic analyses, and memory-efficient blockwise summary-statistic
analyses. Its joint samplers supplement the standard conditional updates with
local blocked updates of connected SNP-gene-cell indicator triplets.

> **Development status:** research software under active development. Results
> should be checked using multiple chains, trace plots, convergence diagnostics,
> and sensitivity analyses before scientific interpretation.

## Model overview

For SNP `j`, gene `g`, and cell type `c`, the model uses binary inclusion
indicators

```text
delta_j  = SNP inclusion indicator
Delta_g  = gene inclusion indicator
gamma_c  = cell-type inclusion indicator
```

The hierarchy connects the layers through:

- `L` (`m x G`): SNP-to-gene annotations
- `A` (`m x C`): SNP-to-cell annotations
- `B` (`G x C`): gene-to-cell annotations

The SNP- and gene-level prior inclusion probabilities are modelled using probit
links:

```text
Pr(delta_j = 1) = Phi(mu_pi + alpha1 * (L Delta)_j
                            + alpha2 * (A gamma)_j
                            + baseline_j)

Pr(Delta_g = 1) = Phi(mu_Pi + alpha3 * (B gamma)_g + pathways_g)
```

where `Phi()` is the standard-normal cumulative distribution function. The
cell indicators have inclusion probability `rho`.

The prior-based functions additionally use

```text
mu_Pi ~ Normal(qnorm(priorPiGene), sigmaMuPiSq)
rho   ~ Beta(priorRho * rhoPriorStrength,
             (1 - priorRho) * rhoPriorStrength)
```

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")

remotes::install_github(
  "mommy003/BayesMAPJoint",
  upgrade = "never"
)
```

Load and verify the package:

```r
library(BayesMAPJoint)
packageVersion("BayesMAPJoint")
```

Compilation requires a working C++ toolchain because the computational kernels
are implemented with Rcpp.

## Main functions

The prior-based joint samplers are recommended for the main analyses.

| Data type | Recommended function |
|---|---|
| Individual-level data | `bayesmap_joint_prior()` |
| Summary statistics with dense LD | `bayesmap_joint_summary_prior()` |
| Summary statistics with LD blocks | `bayesmap_joint_summary_blocks_prior()` |

The following functions are retained as reference or experimental versions:

| Data type | Separate/reference sampler | Joint sampler |
|---|---|---|
| Individual-level data | `bayesmap()` | `bayesmap_joint()` |
| Summary statistics with dense LD | `bayesmap_summary()` | `bayesmap_joint_summary()` |
| Summary statistics with LD blocks | `bayesmap_summary_blocks()` | `bayesmap_joint_summary_blocks()` |

## Input dimensions

Let `n` be the number of individuals, `m` the number of SNPs, `G` the number of
genes, and `C` the number of cell types.

| Object | Meaning | Required dimensions |
|---|---|---|
| `X` | Individual-level genotype matrix | `n x m` |
| `y` | Phenotype vector | length `n` |
| `bhat` | Marginal SNP-effect estimates | length `m` |
| `LD` | Dense SNP LD matrix | `m x m` |
| `LD_blocks` | List of within-block LD matrices | one matrix per block |
| `block_indices` | SNP indices represented by each LD block | list covering `1:m` exactly once |
| `L` | SNP-to-gene annotation matrix | `m x G` |
| `A` | SNP-to-cell annotation matrix | `m x C` |
| `B` | Gene-to-cell annotation matrix | `G x C` |
| `baseline` | Optional SNP-level covariate | length `m` |
| `pathways` | Optional gene-level covariate | length `G` |

The SNP order must be identical across `X` or `bhat`, the LD representation,
`L`, `A`, `baseline`, and `block_indices`.

## Individual-level example

```r
library(BayesMAPJoint)

stopifnot(
  nrow(L) == ncol(X),
  nrow(A) == ncol(X),
  nrow(B) == ncol(L),
  ncol(B) == ncol(A),
  length(y) == nrow(X)
)

set.seed(1)

fit_individual <- bayesmap_joint_prior(
  X = X,
  y = y,
  L = L,
  A = A,
  B = B,
  niter = 4000,
  burnin = 1200,
  thin = 1,
  joint_blocks_per_iter = 5L,
  startPiSnp = 0.01,
  startPiGene = 0.10,
  startRho = 0.10,
  startH2 = 0.20,
  priorPiGene = 0.10,
  priorRho = 0.05,
  rhoPriorStrength = 20,
  sigmaMuPiSq = 1,
  verbose = TRUE
)
```

## Dense summary-statistic example

```r
stopifnot(
  length(bhat) == nrow(L),
  nrow(LD) == length(bhat),
  ncol(LD) == length(bhat)
)

set.seed(1)

fit_summary <- bayesmap_joint_summary_prior(
  bhat = bhat,
  LD = LD,
  N = N,
  L = L,
  A = A,
  B = B,
  phenotype_variance = 1,
  niter = 4000,
  burnin = 1200,
  thin = 1,
  joint_blocks_per_iter = 5L,
  priorPiGene = 0.10,
  priorRho = 0.05,
  rhoPriorStrength = 20,
  sigmaMuPiSq = 1,
  verbose = TRUE
)
```

## Blockwise summary-statistic example

Use the blockwise sampler when a dense `m x m` LD matrix is too large to store.
Every SNP must occur in exactly one element of `block_indices`, and each element
of `LD_blocks` must follow the corresponding SNP order.

```r
stopifnot(
  length(LD_blocks) == length(block_indices),
  identical(sort(unlist(block_indices, use.names = FALSE)), seq_along(bhat)),
  all(vapply(
    seq_along(LD_blocks),
    function(k) {
      nrow(LD_blocks[[k]]) == length(block_indices[[k]]) &&
        ncol(LD_blocks[[k]]) == length(block_indices[[k]])
    },
    logical(1)
  ))
)

set.seed(1)

fit_blocks <- bayesmap_joint_summary_blocks_prior(
  bhat = bhat,
  LD_blocks = LD_blocks,
  block_indices = block_indices,
  N = N,
  L = L,
  A = A,
  B = B,
  phenotype_variance = 1,
  niter = 4000,
  burnin = 1200,
  thin = 1,
  joint_blocks_per_iter = 5L,
  priorPiGene = 0.10,
  priorRho = 0.05,
  rhoPriorStrength = 20,
  sigmaMuPiSq = 1,
  verbose = TRUE
)
```

The blockwise method approximates the full LD matrix by discarding cross-block
LD. Blocks should therefore be constructed so that cross-block correlations are
small.

## Output

The fitted object contains posterior samples and posterior inclusion
probabilities, including:

| Component | Description |
|---|---|
| `par` | MCMC samples of model parameters |
| `beta` | Stored SNP-effect samples, or `NULL` when `store_beta = FALSE` |
| `delta` | Stored SNP-indicator samples, or `NULL` when `store_delta = FALSE` |
| `Delta` | Gene-indicator samples |
| `gamma` | Cell-indicator samples |
| `pip_snp` | SNP posterior inclusion probabilities |
| `pip_gene` | Gene posterior inclusion probabilities |
| `pip_cell` | Cell-type posterior inclusion probabilities |
| `nkeep` | Number of retained MCMC iterations |

The `par` matrix includes quantities such as `SigmaBetaSq`, `Vare`, `Varg`,
`h2`, `Alpha1`, `Alpha2`, `Alpha3`, `Rho`, `Nnz`, `mu_pi`, and `mu_Pi`.

```r
colMeans(fit_summary$par)
head(fit_summary$pip_snp)
head(fit_summary$pip_gene)
fit_summary$pip_cell
```

## MCMC diagnostics

Inspect convergence before reporting posterior summaries. In particular, check
`Nnz`, `h2`, `mu_pi`, `mu_Pi`, and `Rho`.

```r
diagnostic_parameters <- c("Nnz", "h2", "mu_pi", "mu_Pi", "Rho")

old_par <- par(mfrow = c(3, 2))
on.exit(par(old_par), add = TRUE)

for (parameter in diagnostic_parameters) {
  plot(
    fit_summary$par[, parameter],
    type = "l",
    xlab = "Retained iteration",
    ylab = parameter,
    main = parameter
  )
}
```

For final analyses, use longer chains than the short examples, compare multiple
random seeds, and evaluate posterior stability after burn-in.

## Joint versus separate indicator updates

The reference samplers update the SNP, gene, and cell indicators conditionally
in separate steps. The joint samplers retain those updates and additionally
sample a small number of connected `(SNP, gene, cell)` triplets jointly during
each MCMC iteration. This can improve movement between strongly dependent
indicator configurations, but it does not remove the need for convergence
assessment.

`joint_blocks_per_iter` controls how many connected triplets receive a joint
update per iteration. Values around `1` to `5` are useful starting points;
runtime increases as this number increases.

## Development checks

Run the package tests and full package check from the repository root:

```r
devtools::test()
devtools::check()
```

The current test suite contains individual-level, dense-summary, block-summary,
joint-indicator, prior, and LD-block tests.

## Citation

A formal citation will be added when the associated methods manuscript is
available. Until then, please cite the GitHub repository and the version or
commit used in the analysis.

## Authors

- Jian Zeng
- Moksedul Momin

## License

GPL (>= 3)
