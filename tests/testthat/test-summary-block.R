test_that("one LD block reproduces the full LD update", {
  set.seed(123)

  N <- 80L
  m <- 20L

  X <- matrix(
    rnorm(N * m),
    nrow = N,
    ncol = m
  )

  X <- scale(X)

  y <- as.numeric(scale(rnorm(N)))

  bhat <- as.numeric(
    crossprod(X, y) / N
  )

  LD <- crossprod(X) / N
  LD <- (LD + t(LD)) / 2

  beta_initial <- numeric(m)
  delta_initial <- rbinom(m, 1L, 0.05)
  pi_j <- rep(0.05, m)

  phenotype_variance <- as.numeric(
    crossprod(y) / N
  )

  vare <- phenotype_variance * 0.5
  sigmaBetaSq <- phenotype_variance * 0.5 / m

  beta_full <- beta_initial + 0
  delta_full <- as.integer(delta_initial + 0L)

  beta_block <- beta_initial + 0
  delta_block <- as.integer(delta_initial + 0L)

  set.seed(999)

  out_full <- updateBetaDeltaSummary_cpp(
    bhat = bhat,
    LD = LD,
    N = N,
    beta = beta_full,
    delta = delta_full,
    pi_j = pi_j,
    vare = vare,
    sigmaBetaSq = sigmaBetaSq
  )

  set.seed(999)

  out_block <- updateBetaDeltaSummaryBlocks(
    bhat = bhat,
    LD_blocks = list(LD),
    block_indices = list(seq_len(m)),
    N = N,
    beta = beta_block,
    delta = delta_block,
    pi_j = pi_j,
    vare = vare,
    sigmaBetaSq = sigmaBetaSq
  )

  expect_identical(
    out_block$delta,
    out_full$delta
  )

  expect_equal(
    out_block$beta,
    out_full$beta,
    tolerance = 1e-12
  )

  expect_equal(
    out_block$rcorr,
    out_full$rcorr,
    tolerance = 1e-12
  )

  expect_equal(
    out_block$nnz,
    out_full$nnz
  )
})



test_that("multiple LD blocks reproduce a block-diagonal full LD update", {
  set.seed(123)

  N <- 100L
  m1 <- 10L
  m2 <- 12L
  m <- m1 + m2

  ## Generate two independent genotype blocks
  X1 <- matrix(
    rnorm(N * m1),
    nrow = N,
    ncol = m1
  )

  X2 <- matrix(
    rnorm(N * m2),
    nrow = N,
    ncol = m2
  )

  X1 <- scale(X1)
  X2 <- scale(X2)

  y <- as.numeric(scale(rnorm(N)))

  bhat1 <- as.numeric(crossprod(X1, y) / N)
  bhat2 <- as.numeric(crossprod(X2, y) / N)

  bhat <- c(bhat1, bhat2)

  LD1 <- crossprod(X1) / N
  LD1 <- (LD1 + t(LD1)) / 2

  LD2 <- crossprod(X2) / N
  LD2 <- (LD2 + t(LD2)) / 2

  ## Construct exact block-diagonal full LD matrix
  LD_full <- matrix(
    0,
    nrow = m,
    ncol = m
  )

  index1 <- seq_len(m1)
  index2 <- m1 + seq_len(m2)

  LD_full[index1, index1] <- LD1
  LD_full[index2, index2] <- LD2

  beta_initial <- numeric(m)

  delta_initial <- rbinom(
    m,
    size = 1L,
    prob = 0.05
  )

  pi_j <- rep(0.05, m)

  phenotype_variance <- as.numeric(
    crossprod(y) / N
  )

  vare <- phenotype_variance * 0.5
  sigmaBetaSq <- phenotype_variance * 0.5 / m

  ## Independent copies because Rcpp may modify vectors in place
  beta_full <- beta_initial + 0
  delta_full <- as.integer(delta_initial + 0L)

  beta_blocks <- beta_initial + 0
  delta_blocks <- as.integer(delta_initial + 0L)

  set.seed(999)

  out_full <- updateBetaDeltaSummary_cpp(
    bhat = bhat,
    LD = LD_full,
    N = N,
    beta = beta_full,
    delta = delta_full,
    pi_j = pi_j,
    vare = vare,
    sigmaBetaSq = sigmaBetaSq
  )

  set.seed(999)

  out_blocks <- updateBetaDeltaSummaryBlocks(
    bhat = bhat,
    LD_blocks = list(LD1, LD2),
    block_indices = list(index1, index2),
    N = N,
    beta = beta_blocks,
    delta = delta_blocks,
    pi_j = pi_j,
    vare = vare,
    sigmaBetaSq = sigmaBetaSq
  )

  expect_identical(
    out_blocks$delta,
    out_full$delta
  )

  expect_equal(
    out_blocks$beta,
    out_full$beta,
    tolerance = 1e-12
  )

  expect_equal(
    out_blocks$rcorr,
    out_full$rcorr,
    tolerance = 1e-12
  )

  expect_equal(
    out_blocks$nnz,
    out_full$nnz
  )

  rcorr_direct <- bhat -
    as.numeric(
      LD_full %*% (
        out_blocks$beta *
          out_blocks$delta
      )
    )

  expect_equal(
    out_blocks$rcorr,
    rcorr_direct,
    tolerance = 1e-12
  )
})
