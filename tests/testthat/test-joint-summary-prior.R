test_that("joint summary prior retains sparse-prior controls", {
  set.seed(456)

  N <- 60L
  m <- 15L
  G <- 5L
  C <- 3L
  X_test <- scale(matrix(rnorm(N * m), nrow = N))
  y_test <- as.numeric(scale(rnorm(N)))
  bhat_test <- as.numeric(crossprod(X_test, y_test) / N)
  LD_test <- crossprod(X_test) / N
  L_test <- matrix(rbinom(m * G, 1, 0.20), m, G)
  A_test <- matrix(rbinom(m * C, 1, 0.20), m, C)
  B_test <- matrix(rbinom(G * C, 1, 0.30), G, C)
  L_test[1, 1] <- A_test[1, 1] <- B_test[1, 1] <- 1

  fit <- bayesmap_joint_summary_prior(
    bhat = bhat_test,
    LD = LD_test,
    N = N,
    L = L_test,
    A = A_test,
    B = B_test,
    phenotype_variance = as.numeric(crossprod(y_test) / N),
    niter = 20,
    burnin = 5,
    joint_blocks_per_iter = 1L,
    priorPiGene = 0.01,
    priorRho = 0.05,
    rhoPriorStrength = 20,
    store_beta = FALSE,
    store_delta = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "BayesMAPJointSummaryPrior")
  expect_equal(fit$prior$priorPiGene, 0.01)
  expect_equal(fit$prior$priorRho, 0.05)
  expect_equal(fit$prior$rhoPriorStrength, 20)
  expect_true(all(fit$pip_snp >= 0 & fit$pip_snp <= 1))
  expect_true(all(fit$pip_gene >= 0 & fit$pip_gene <= 1))
  expect_true(all(fit$pip_cell >= 0 & fit$pip_cell <= 1))
})
