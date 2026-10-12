test_that("joint indicator kernel preserves dimensions and residual", {
  set.seed(321)

  n <- 40L
  m <- 12L
  G <- 4L
  C <- 3L

  X_joint <- scale(matrix(rnorm(n * m), nrow = n))
  y_joint <- as.numeric(scale(rnorm(n)))
  mu_joint <- mean(y_joint)
  beta_joint <- numeric(m)
  delta_joint <- rbinom(m, 1L, 0.20)
  Delta_joint <- rbinom(G, 1L, 0.25)
  gamma_joint <- rbinom(C, 1L, 0.30)

  L_joint <- matrix(0, m, G)
  A_joint <- matrix(0, m, C)
  B_joint <- matrix(0, G, C)
  L_joint[cbind(seq_len(m), rep(seq_len(G), length.out = m))] <- 1
  A_joint[cbind(seq_len(m), rep(seq_len(C), length.out = m))] <- 1
  B_joint[cbind(seq_len(G), rep(seq_len(C), length.out = G))] <- 1

  ycorr_joint <- y_joint - mu_joint -
    as.vector(X_joint %*% (beta_joint * delta_joint))

  out <- updateJointIndicators_cpp(
    X = X_joint,
    ycorr = ycorr_joint,
    beta = beta_joint,
    delta = delta_joint,
    Delta = Delta_joint,
    gamma = gamma_joint,
    L = L_joint,
    A = A_joint,
    B = B_joint,
    snp_index = c(1L, 2L),
    gene_index = c(1L, 2L),
    cell_index = c(1L, 2L),
    xpx = colSums(X_joint^2),
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

  direct_residual <- y_joint - mu_joint -
    as.vector(X_joint %*% (out$beta * out$delta))

  expect_equal(out$ycorr, direct_residual, tolerance = 1e-10)
})


test_that("bayesmap_joint returns the experimental class", {
  set.seed(123)

  fit <- bayesmap_joint(
    X = X,
    y = y,
    L = L,
    A = A,
    B = B,
    niter = 10,
    joint_blocks_per_iter = 1L,
    verbose = FALSE
  )

  expect_s3_class(fit, "BayesMAPJoint")
  expect_equal(fit$joint_blocks_per_iter, 1L)
  expect_true(all(fit$pip_snp >= 0 & fit$pip_snp <= 1))
  expect_true(all(fit$pip_gene >= 0 & fit$pip_gene <= 1))
  expect_true(all(fit$pip_cell >= 0 & fit$pip_cell <= 1))
})


test_that("bayesmap_joint_prior retains sparse-prior controls", {
  set.seed(456)

  fit <- bayesmap_joint_prior(
    X = X,
    y = y,
    L = L,
    A = A,
    B = B,
    niter = 10,
    joint_blocks_per_iter = 1L,
    priorPiGene = 0.01,
    priorRho = 0.05,
    rhoPriorStrength = 20,
    verbose = FALSE
  )

  expect_s3_class(fit, "BayesMAP")
  expect_equal(fit$joint_blocks_per_iter, 1L)
  expect_equal(fit$prior$priorPiGene, 0.01)
  expect_equal(fit$prior$priorRho, 0.05)
  expect_equal(fit$prior$rhoPriorStrength, 20)
})
