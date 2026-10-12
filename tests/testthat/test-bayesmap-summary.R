test_that("bayesmap_summary returns valid posterior output", {
  set.seed(123)

  N <- 80L
  m <- 20L
  G <- 6L
  C <- 3L

  X <- matrix(
    rnorm(N * m),
    nrow = N,
    ncol = m
  )

  X <- scale(
    X,
    center = TRUE,
    scale = TRUE
  )

  beta_true <- numeric(m)
  beta_true[c(2, 7, 14)] <- c(0.20, -0.15, 0.18)

  y <- as.numeric(
    X %*% beta_true +
      rnorm(N)
  )

  y <- as.numeric(
    scale(
      y,
      center = TRUE,
      scale = TRUE
    )
  )

  bhat <- as.numeric(
    crossprod(X, y) / N
  )

  LD <- crossprod(X) / N
  LD <- (LD + t(LD)) / 2

  phenotype_variance <- as.numeric(
    crossprod(y) / N
  )

  L <- matrix(
    rbinom(m * G, 1, 0.15),
    nrow = m,
    ncol = G
  )

  A <- matrix(
    rbinom(m * C, 1, 0.20),
    nrow = m,
    ncol = C
  )

  B <- matrix(
    rbinom(G * C, 1, 0.25),
    nrow = G,
    ncol = C
  )

  set.seed(999)

  fit <- bayesmap_summary(
    bhat = bhat,
    LD = LD,
    N = N,
    L = L,
    A = A,
    B = B,
    phenotype_variance = phenotype_variance,
    niter = 40,
    burnin = 10,
    thin = 2,
    store_beta = FALSE,
    store_delta = FALSE,
    verbose = FALSE
  )

  expect_s3_class(
    fit,
    "BayesMAPSummary"
  )

  expect_equal(
    fit$nkeep,
    15L
  )

  expect_equal(
    dim(fit$par),
    c(15L, 11L)
  )

  expect_null(
    fit$beta
  )

  expect_null(
    fit$delta
  )

  expect_equal(
    dim(fit$Delta),
    c(15L, G)
  )

  expect_equal(
    dim(fit$gamma),
    c(15L, C)
  )

  expect_true(
    all(is.finite(fit$par))
  )

  expect_true(
    all(fit$pip_snp >= 0 & fit$pip_snp <= 1)
  )

  expect_true(
    all(fit$pip_gene >= 0 & fit$pip_gene <= 1)
  )

  expect_true(
    all(fit$pip_cell >= 0 & fit$pip_cell <= 1)
  )

  expect_true(
    all(fit$par[, "Vare"] > 0)
  )

  expect_true(
    all(fit$par[, "Varg"] >= 0)
  )

  expect_true(
    all(fit$par[, "h2"] >= 0 &
          fit$par[, "h2"] <= 1)
  )
})
