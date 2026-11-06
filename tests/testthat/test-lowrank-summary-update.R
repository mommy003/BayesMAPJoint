test_that(
  "full-rank kernel reproduces the exact summary update",
  {
    set.seed(20261106)
    
    n_reference <- 300L
    m <- 12L
    
    X <- matrix(
      rnorm(n_reference * m),
      nrow = n_reference,
      ncol = m
    )
    
    LD <- cor(X)
    LD <- (LD + t(LD)) / 2
    
    LD_eigen <- eigen(
      LD,
      symmetric = TRUE
    )
    
    expect_true(
      all(LD_eigen$values > 0)
    )
    
    lowrank_block <-
      make_lowrank_ld_block(
        U = LD_eigen$vectors,
        lambda = LD_eigen$values,
        block_indices = seq_len(m),
        target_diagonal = diag(LD),
        variance_explained = 1
      )
    
    bhat <- rnorm(
      m,
      sd = 0.02
    )
    
    beta_start <- numeric(m)
    delta_start <- integer(m)
    pi_j <- rep(0.1, m)
    
    set.seed(12345)
    
    exact_update <-
      updateBetaDeltaSummary_cpp(
        bhat = bhat,
        LD = LD,
        N = 5000,
        beta = beta_start + 0,
        delta = delta_start + 0L,
        pi_j = pi_j,
        vare = 0.8,
        sigmaBetaSq = 0.001
      )
    
    set.seed(12345)
    
    lowrank_update <-
      updateBetaDeltaSummaryLowRank_cpp(
        bhat = bhat,
        U = lowrank_block$U,
        lambda = lowrank_block$lambda,
        diagonal_correction =
          lowrank_block$diagonal_correction,
        N = 5000,
        beta = beta_start + 0,
        delta = delta_start + 0L,
        pi_j = pi_j,
        vare = 0.8,
        sigmaBetaSq = 0.001
      )
    
    expect_identical(
      lowrank_update$delta,
      exact_update$delta
    )
    
    expect_identical(
      lowrank_update$nnz,
      exact_update$nnz
    )
    
    expect_equal(
      lowrank_update$beta,
      exact_update$beta,
      tolerance = 1e-12
    )
    
    expect_equal(
      lowrank_update$rcorr,
      exact_update$rcorr,
      tolerance = 1e-12
    )
  }
)


test_that(
  "low-rank summary kernel validates dimensions",
  {
    expect_error(
      updateBetaDeltaSummaryLowRank_cpp(
        bhat = numeric(4),
        U = matrix(0, nrow = 3, ncol = 2),
        lambda = c(1, 0.5),
        diagonal_correction = numeric(4),
        N = 1000,
        beta = numeric(4),
        delta = integer(4),
        pi_j = rep(0.1, 4),
        vare = 1,
        sigmaBetaSq = 0.01
      ),
      "U must have"
    )
  }
)



test_that(
  "full-rank block operations reproduce exact LD blocks",
  {
    set.seed(20261106)
    
    block_indices <- list(
      1:6,
      7:13
    )
    
    make_correlation <- function(m) {
      X <- matrix(
        rnorm(300L * m),
        nrow = 300L,
        ncol = m
      )
      
      LD <- cor(X)
      
      (LD + t(LD)) / 2
    }
    
    exact_blocks <- list(
      make_correlation(6L),
      make_correlation(7L)
    )
    
    lowrank_blocks <- lapply(
      seq_along(exact_blocks),
      function(block_number) {
        LD <- exact_blocks[[block_number]]
        
        decomposition <- eigen(
          LD,
          symmetric = TRUE
        )
        
        make_lowrank_ld_block(
          U = decomposition$vectors,
          lambda = decomposition$values,
          block_indices =
            block_indices[[block_number]],
          target_diagonal = diag(LD),
          variance_explained = 1
        )
      }
    )
    
    expect_true(
      validate_summary_lowrank_blocks(
        lowrank_blocks,
        m = 13L
      )
    )
    
    test_vector <- rnorm(13L)
    
    exact_product <-
      summary_ld_block_matvec(
        vector = test_vector,
        LD_blocks = exact_blocks,
        block_indices = block_indices
      )
    
    lowrank_product <-
      summary_lowrank_block_matvec(
        vector = test_vector,
        lowrank_blocks = lowrank_blocks
      )
    
    expect_equal(
      lowrank_product,
      exact_product,
      tolerance = 1e-12
    )
    
    bhat <- rnorm(
      13L,
      sd = 0.02
    )
    
    beta_start <- numeric(13L)
    delta_start <- integer(13L)
    pi_j <- rep(0.1, 13L)
    
    set.seed(54321)
    
    exact_update <-
      updateBetaDeltaSummaryBlocks(
        bhat = bhat,
        LD_blocks = exact_blocks,
        block_indices = block_indices,
        N = 5000,
        beta = beta_start + 0,
        delta = delta_start + 0L,
        pi_j = pi_j,
        vare = 0.8,
        sigmaBetaSq = 0.001
      )
    
    set.seed(54321)
    
    lowrank_update <-
      updateBetaDeltaSummaryLowRankBlocks(
        bhat = bhat,
        lowrank_blocks = lowrank_blocks,
        N = 5000,
        beta = beta_start + 0,
        delta = delta_start + 0L,
        pi_j = pi_j,
        vare = 0.8,
        sigmaBetaSq = 0.001
      )
    
    expect_identical(
      lowrank_update$delta,
      exact_update$delta
    )
    
    expect_identical(
      lowrank_update$nnz,
      exact_update$nnz
    )
    
    expect_equal(
      lowrank_update$beta,
      exact_update$beta,
      tolerance = 1e-12
    )
    
    expect_equal(
      lowrank_update$rcorr,
      exact_update$rcorr,
      tolerance = 1e-12
    )
  }
)