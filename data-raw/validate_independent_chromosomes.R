## ================================================================
## Validate dense versus blockwise summary BayesMAP
## using two approximately independent chromosomes
## ================================================================

set.seed(2026)

## ------------------------------------------------
## 1. Simulation dimensions
## ------------------------------------------------

N <- 10000L

m1 <- 100L
m2 <- 100L
m <- m1 + m2

G <- 40L
C <- 6L

## Within-chromosome LD decay parameters
rho1 <- 0.80
rho2 <- 0.70


## ------------------------------------------------
## 2. Create chromosome-specific LD structures
## ------------------------------------------------

make_ar1_correlation <- function(size, rho) {
  index <- seq_len(size)

  rho^abs(
    outer(index, index, FUN = "-")
  )
}

Sigma1 <- make_ar1_correlation(
  size = m1,
  rho = rho1
)

Sigma2 <- make_ar1_correlation(
  size = m2,
  rho = rho2
)


## ------------------------------------------------
## 3. Simulate independent chromosome genotypes
## ------------------------------------------------

Z1 <- matrix(
  rnorm(N * m1),
  nrow = N,
  ncol = m1
)

Z2 <- matrix(
  rnorm(N * m2),
  nrow = N,
  ncol = m2
)

X1 <- Z1 %*% chol(Sigma1)
X2 <- Z2 %*% chol(Sigma2)

X1 <- scale(
  X1,
  center = TRUE,
  scale = TRUE
)

X2 <- scale(
  X2,
  center = TRUE,
  scale = TRUE
)

X <- cbind(X1, X2)

storage.mode(X) <- "double"


## ------------------------------------------------
## 4. Chromosome and position information
## ------------------------------------------------

chr <- c(
  rep("1", m1),
  rep("2", m2)
)

position <- c(
  seq(
    from = 1e6,
    by = 10000,
    length.out = m1
  ),
  seq(
    from = 1e6,
    by = 10000,
    length.out = m2
  )
)

block_table <- data.frame(
  chr = c("1", "2"),
  start = c(
    min(position[chr == "1"]),
    min(position[chr == "2"])
  ),
  end = c(
    max(position[chr == "1"]),
    max(position[chr == "2"])
  )
)


## ------------------------------------------------
## 5. Simulate causal effects
## ------------------------------------------------

causal_chr1 <- sample(
  seq_len(m1),
  size = 8
)

causal_chr2 <- m1 + sample(
  seq_len(m2),
  size = 8
)

causal <- c(
  causal_chr1,
  causal_chr2
)

beta_true <- numeric(m)

beta_true[causal] <- rnorm(
  length(causal),
  mean = 0,
  sd = 0.12
)


## ------------------------------------------------
## 6. Simulate standardized phenotype
## ------------------------------------------------

genetic_value <- as.numeric(
  X %*% beta_true
)

target_h2 <- 0.50

genetic_variance <- var(genetic_value)

environmental_variance <-
  genetic_variance *
  (1 - target_h2) /
  target_h2

y <- genetic_value +
  rnorm(
    N,
    mean = 0,
    sd = sqrt(environmental_variance)
  )

y <- as.numeric(
  scale(
    y,
    center = TRUE,
    scale = TRUE
  )
)

phenotype_variance <- as.numeric(
  crossprod(y) / N
)


## ------------------------------------------------
## 7. Generate summary statistics and full LD
## ------------------------------------------------

bhat <- as.numeric(
  crossprod(X, y) / N
)

LD_full <- crossprod(X) / N
LD_full <- (LD_full + t(LD_full)) / 2


## ------------------------------------------------
## 8. Build LD blocks using chromosome and position
## ------------------------------------------------

block_info <- make_ld_blocks(
  chr = chr,
  position = position,
  block_table = block_table
)

block_indices <- block_info$block_indices

if (length(block_info$unassigned) > 0L) {
  stop("Some SNPs were not assigned to an LD block.")
}

LD_blocks <- extract_ld_blocks(
  LD = LD_full,
  block_indices = block_indices
)

block_info$block_indices
lengths(block_indices)
block_info$unassigned


## ------------------------------------------------
## 9. Assess discarded cross-chromosome LD
## ------------------------------------------------

ld_diagnostics <- assess_ld_blocks(
  LD = LD_full,
  block_indices = block_indices
)

print(ld_diagnostics)

check_ld_block_quality(
  LD = LD_full,
  block_indices = block_indices,
  frobenius_warning = 0.10,
  max_cross_r_warning = 0.15
)


## ------------------------------------------------
## 10. Simulate SNP-gene-cell annotation matrices
## ------------------------------------------------

L <- matrix(
  rbinom(
    m * G,
    size = 1L,
    prob = 0.04
  ),
  nrow = m,
  ncol = G
)

A <- matrix(
  rbinom(
    m * C,
    size = 1L,
    prob = 0.08
  ),
  nrow = m,
  ncol = C
)

B <- matrix(
  rbinom(
    G * C,
    size = 1L,
    prob = 0.15
  ),
  nrow = G,
  ncol = C
)

baseline <- numeric(m)
pathways <- numeric(G)


## ------------------------------------------------
## 11. Run dense summary BayesMAP
## ------------------------------------------------

set.seed(123)

time_dense <- system.time({
  fit_dense_chr <- bayesmap_summary_R(
    bhat = bhat,
    LD = LD_full,
    N = N,
    L = L,
    A = A,
    B = B,
    baseline = baseline,
    pathways = pathways,
    phenotype_variance = phenotype_variance,
    niter = 10000,
    burnin = 2000,
    thin = 5,
    store_beta = FALSE,
    store_delta = FALSE,
    verbose = FALSE
  )
})[["elapsed"]]


## ------------------------------------------------
## 12. Run blockwise summary BayesMAP
## ------------------------------------------------

set.seed(123)

time_blocks <- system.time({
  fit_blocks_chr <- bayesmap_summary_blocks_R(
    bhat = bhat,
    LD_blocks = LD_blocks,
    block_indices = block_indices,
    N = N,
    L = L,
    A = A,
    B = B,
    baseline = baseline,
    pathways = pathways,
    phenotype_variance = phenotype_variance,
    niter = 5000,
    burnin = 1000,
    thin = 5,
    store_beta = FALSE,
    store_delta = FALSE,
    verbose = FALSE
  )
})[["elapsed"]]


## ------------------------------------------------
## 13. Compare posterior means
## ------------------------------------------------

dense_means <- colMeans(
  fit_dense_chr$par
)

block_means <- colMeans(
  fit_blocks_chr$par
)

parameter_comparison <- cbind(
  Dense = dense_means,
  Blocks = block_means,
  Difference = block_means - dense_means
)

print(
  round(
    parameter_comparison,
    5
  )
)


## ------------------------------------------------
## 14. Compare PIPs
## ------------------------------------------------

pip_correlations <- c(
  SNP = cor(
    fit_dense_chr$pip_snp,
    fit_blocks_chr$pip_snp
  ),
  Gene = cor(
    fit_dense_chr$pip_gene,
    fit_blocks_chr$pip_gene
  ),
  Cell = cor(
    fit_dense_chr$pip_cell,
    fit_blocks_chr$pip_cell
  )
)

print(pip_correlations)


pip_differences <- rbind(
  SNP = c(
    mean_abs_difference = mean(
      abs(
        fit_dense_chr$pip_snp -
          fit_blocks_chr$pip_snp
      )
    ),
    max_abs_difference = max(
      abs(
        fit_dense_chr$pip_snp -
          fit_blocks_chr$pip_snp
      )
    )
  ),
  Gene = c(
    mean_abs_difference = mean(
      abs(
        fit_dense_chr$pip_gene -
          fit_blocks_chr$pip_gene
      )
    ),
    max_abs_difference = max(
      abs(
        fit_dense_chr$pip_gene -
          fit_blocks_chr$pip_gene
      )
    )
  ),
  Cell = c(
    mean_abs_difference = mean(
      abs(
        fit_dense_chr$pip_cell -
          fit_blocks_chr$pip_cell
      )
    ),
    max_abs_difference = max(
      abs(
        fit_dense_chr$pip_cell -
          fit_blocks_chr$pip_cell
      )
    )
  )
)

print(pip_differences)


## ------------------------------------------------
## 15. Compare runtime
## ------------------------------------------------

runtime_comparison <- c(
  Dense_seconds = time_dense,
  Block_seconds = time_blocks,
  Speedup = time_dense / time_blocks,
  Runtime_reduction_percent =
    100 * (1 - time_blocks / time_dense)
)

print(runtime_comparison)


## ------------------------------------------------
## 16. Check causal-SNP recovery
## ------------------------------------------------

causal_results <- data.frame(
  SNP = causal,
  Chromosome = chr[causal],
  Position = position[causal],
  True_beta = beta_true[causal],
  Dense_PIP = fit_dense_chr$pip_snp[causal],
  Block_PIP = fit_blocks_chr$pip_snp[causal]
)

causal_results <- causal_results[
  order(
    causal_results$Dense_PIP,
    decreasing = TRUE
  ),
]

print(causal_results)


## ------------------------------------------------
## 17. Basic validity checks
## ------------------------------------------------

stopifnot(
  all(is.finite(fit_dense_chr$par)),
  all(is.finite(fit_blocks_chr$par)),
  all(fit_dense_chr$pip_snp >= 0),
  all(fit_dense_chr$pip_snp <= 1),
  all(fit_blocks_chr$pip_snp >= 0),
  all(fit_blocks_chr$pip_snp <= 1)
)

cat(
  "\nIndependent-chromosome validation completed successfully.\n"
)
