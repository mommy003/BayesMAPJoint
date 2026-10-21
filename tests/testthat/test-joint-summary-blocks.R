make_joint_summary_block_test_data <- function() {
  m <- 8L
  G <- 3L
  C <- 2L

  LD1 <- toeplitz(0.20^(0:3))
  LD2 <- toeplitz(0.15^(0:3))

  L <- matrix(0, nrow = m, ncol = G)

  gene_assignment <- rep(
    seq_len(G),
    length.out = m
  )

  L[cbind(seq_len(m), gene_assignment)] <- 1

  A <- matrix(0, nrow = m, ncol = C)

  cell_assignment <- rep(
    seq_len(C),
    length.out = m
  )

  A[cbind(seq_len(m), cell_assignment)] <- 1

  B <- matrix(
    c(
      1, 0,
      1, 1,
      0, 1
    ),
    nrow = G,
    ncol = C,
    byrow = TRUE
  )

  list(
    bhat = c(
      0.05, -0.03, 0.02, 0.01,
      -0.04, 0.03, -0.02, 0.04
    ),
    LD_blocks = list(LD1, LD2),
    block_indices = list(
      1:4,
      5:8
    ),
    N = 1000,
    L = L,
    A = A,
    B = B
  )
}


test_that(
  "joint summary block wrapper produces valid output",
  {
    dat <- make_joint_summary_block_test_data()

    set.seed(101)

    fit <- bayesmap_joint_summary_blocks(
      bhat = dat$bhat,
      LD_blocks = dat$LD_blocks,
      block_indices = dat$block_indices,
      N = dat$N,
      L = dat$L,
      A = dat$A,
      B = dat$B,
      niter = 30,
      burnin = 10,
      joint_blocks_per_iter = 2L,
      verbose = FALSE
    )

    expect_s3_class(
      fit,
      "BayesMAPJointSummaryBlocks"
    )

    expect_s3_class(
      fit,
      "BayesMAPSummaryBlocks"
    )

    expect_equal(fit$nkeep, 20L)
    expect_equal(fit$n_blocks, 2L)
    expect_equal(
      fit$joint_blocks_per_iter,
      2L
    )

    expect_length(fit$pip_snp, 8L)
    expect_length(fit$pip_gene, 3L)
    expect_length(fit$pip_cell, 2L)

    expect_true(all(is.finite(fit$par)))
    expect_true(all(fit$pip_snp >= 0))
    expect_true(all(fit$pip_snp <= 1))
  }
)


test_that(
  "joint summary block prior wrapper records its priors",
  {
    dat <- make_joint_summary_block_test_data()

    set.seed(102)

    fit <- bayesmap_joint_summary_blocks_prior(
      bhat = dat$bhat,
      LD_blocks = dat$LD_blocks,
      block_indices = dat$block_indices,
      N = dat$N,
      L = dat$L,
      A = dat$A,
      B = dat$B,
      niter = 30,
      burnin = 10,
      joint_blocks_per_iter = 2L,
      priorPiGene = 0.10,
      priorRho = 0.05,
      rhoPriorStrength = 20,
      sigmaMuPiSq = 1,
      verbose = FALSE
    )

    expect_s3_class(
      fit,
      "BayesMAPJointSummaryBlocksPrior"
    )

    expect_s3_class(
      fit,
      "BayesMAPSummaryBlocks"
    )

    expect_equal(
      fit$joint_blocks_per_iter,
      2L
    )

    expect_equal(fit$priorPiGene, 0.10)
    expect_equal(fit$priorRho, 0.05)
    expect_equal(
      fit$rhoPriorStrength,
      20
    )

    expect_true(all(is.finite(fit$par)))
    expect_true(all(fit$pip_gene >= 0))
    expect_true(all(fit$pip_gene <= 1))
    expect_true(all(fit$pip_cell >= 0))
    expect_true(all(fit$pip_cell <= 1))
  }
)
