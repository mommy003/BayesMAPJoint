test_that("summary joint kernel preserves dimensions and residual", {
  set.seed(321)

  N <- 60L
  m <- 12L
  G <- 4L
  C <- 3L
  X_test <- scale(matrix(rnorm(N * m), nrow = N))
  y_test <- as.numeric(scale(rnorm(N)))
  bhat_test <- as.numeric(crossprod(X_test, y_test) / N)
  LD_test <- crossprod(X_test) / N

  beta_test <- numeric(m)
  delta_test <- rbinom(m, 1L, 0.20)
  Delta_test <- rbinom(G, 1L, 0.25)
  gamma_test <- rbinom(C, 1L, 0.30)
  L_test <- matrix(0, m, G)
  A_test <- matrix(0, m, C)
  B_test <- matrix(0, G, C)
  L_test[cbind(seq_len(m), rep(seq_len(G), length.out = m))] <- 1
  A_test[cbind(seq_len(m), rep(seq_len(C), length.out = m))] <- 1
  B_test[cbind(seq_len(G), rep(seq_len(C), length.out = G))] <- 1

  out <- updateJointIndicatorsSummary_cpp(
    bhat = bhat_test,
    LD = LD_test,
    N = N,
    beta = beta_test,
    delta = delta_test,
    Delta = Delta_test,
    gamma = gamma_test,
    L = L_test,
    A = A_test,
    B = B_test,
    snp_index = c(1L, 2L),
    gene_index = c(1L, 2L),
    cell_index = c(1L, 2L),
    vare = 0.8,
    sigmaBetaSq = 0.05,
    alpha1 = 0.3,
    alpha2 = 0.2,
    alpha3 = 0.4,
    mu_pi = qnorm(0.05),
    mu_Pi = qnorm(0.10),
    baseline = numeric(m),
    pathways = numeric(G),
    rho = 0.10
  )

  expect_length(out$beta, m)
  expect_length(out$delta, m)
  expect_length(out$Delta, G)
  expect_length(out$gamma, C)
  expect_true(all(out$delta %in% 0:1))
  expect_true(all(out$Delta %in% 0:1))
  expect_true(all(out$gamma %in% 0:1))
  expect_equal(out$nnz, sum(out$delta))
  expect_equal(dim(out$sampled_states), c(2L, 3L))
  expect_equal(
    out$rcorr,
    bhat_test - as.vector(LD_test %*% (out$beta * out$delta)),
    tolerance = 1e-10
  )
})


test_that("joint summary wrapper returns valid posterior output", {
  set.seed(123)

  N <- 80L
  m <- 20L
  G <- 6L
  C <- 3L
  X_test <- scale(matrix(rnorm(N * m), nrow = N))
  y_test <- as.numeric(scale(rnorm(N)))
  bhat_test <- as.numeric(crossprod(X_test, y_test) / N)
  LD_test <- crossprod(X_test) / N
  L_test <- matrix(rbinom(m * G, 1, 0.20), m, G)
  A_test <- matrix(rbinom(m * C, 1, 0.20), m, C)
  B_test <- matrix(rbinom(G * C, 1, 0.30), G, C)

  # Guarantee at least one valid connected triplet.
  L_test[1, 1] <- A_test[1, 1] <- B_test[1, 1] <- 1

  fit <- bayesmap_joint_summary(
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
    store_beta = FALSE,
    store_delta = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "BayesMAPJointSummary")
  expect_s3_class(fit, "BayesMAPSummary")
  expect_equal(fit$joint_blocks_per_iter, 1L)
  expect_true(all(fit$pip_snp >= 0 & fit$pip_snp <= 1))
  expect_true(all(fit$pip_gene >= 0 & fit$pip_gene <= 1))
  expect_true(all(fit$pip_cell >= 0 & fit$pip_cell <= 1))
})
