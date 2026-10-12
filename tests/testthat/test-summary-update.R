test_that("summary C++ update matches R reference", {
  set.seed(123)

  N <- 100L
  m <- 20L

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

  y <- rnorm(N)
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

  beta_start <- numeric(m)

  delta_start <- rbinom(
    m,
    size = 1L,
    prob = 0.05
  )

  pi_j_start <- rep(
    0.05,
    m
  )

  phenotype_variance <- as.numeric(
    crossprod(y) / N
  )

  vare_start <- phenotype_variance * 0.5

  sigmaBetaSq_start <-
    phenotype_variance * 0.5 / m

  set.seed(999)

  out_R <- updateBetaDeltaSummary_R(
    bhat = bhat,
    LD = LD,
    N = N,
    beta = beta_start,
    delta = delta_start,
    pi_j = pi_j_start,
    vare = vare_start,
    sigmaBetaSq = sigmaBetaSq_start
  )

  set.seed(999)

  out_cpp <- updateBetaDeltaSummary_cpp(
    bhat = bhat,
    LD = LD,
    N = N,
    beta = beta_start,
    delta = delta_start,
    pi_j = pi_j_start,
    vare = vare_start,
    sigmaBetaSq = sigmaBetaSq_start
  )

  expect_identical(
    out_cpp$delta,
    out_R$delta
  )

  expect_equal(
    out_cpp$beta,
    out_R$beta,
    tolerance = 1e-12
  )

  expect_equal(
    out_cpp$rcorr,
    out_R$rcorr,
    tolerance = 1e-12
  )

  expect_equal(
    out_cpp$nnz,
    out_R$nnz
  )

  rcorr_direct <- bhat -
    as.numeric(
      LD %*% (
        out_cpp$beta *
          out_cpp$delta
      )
    )

  expect_equal(
    out_cpp$rcorr,
    rcorr_direct,
    tolerance = 1e-12
  )
})
