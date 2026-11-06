test_that(
  "full-rank LD reproduces the exact block MCMC chain",
  {
    set.seed(20261106)
    
    m <- 8L
    G <- 2L
    C <- 2L
    
    X <- matrix(
      rnorm(300L * m),
      nrow = 300L,
      ncol = m
    )
    
    LD <- cor(X)
    LD <- (LD + t(LD)) / 2
    
    decomposition <- eigen(
      LD,
      symmetric = TRUE
    )
    
    lowrank_block <-
      make_lowrank_ld_block(
        U = decomposition$vectors,
        lambda = decomposition$values,
        block_indices = seq_len(m),
        target_diagonal = diag(LD),
        variance_explained = 1
      )
    
    L <- matrix(
      c(
        1, 0,
        1, 0,
        1, 0,
        1, 0,
        0, 1,
        0, 1,
        0, 1,
        0, 1
      ),
      nrow = m,
      ncol = G,
      byrow = TRUE
    )
    
    A <- matrix(
      rep(
        c(1, 0, 0, 1),
        length.out = m * C
      ),
      nrow = m,
      ncol = C,
      byrow = TRUE
    )
    
    B <- diag(2)
    
    bhat <- rnorm(
      m,
      sd = 0.02
    )
    
    common_arguments <- list(
      bhat = bhat,
      N = 2000,
      L = L,
      A = A,
      B = B,
      phenotype_variance = 1,
      estimate_vare = FALSE,
      fixed_vare = 0.8,
      niter = 20,
      burnin = 5,
      thin = 1,
      store_beta = TRUE,
      store_delta = TRUE,
      startH2 = 0.2,
      joint_blocks_per_iter = 0L,
      verbose = FALSE
    )
    
    set.seed(777)
    
    exact_fit <- do.call(
      bayesmap_summary_blocks,
      c(
        common_arguments,
        list(
          LD_blocks = list(LD),
          block_indices =
            list(seq_len(m)),
          lowrank_blocks = NULL
        )
      )
    )
    
    set.seed(777)
    
    lowrank_fit <- do.call(
      bayesmap_summary_blocks,
      c(
        common_arguments,
        list(
          LD_blocks = NULL,
          block_indices = NULL,
          lowrank_blocks =
            list(lowrank_block)
        )
      )
    )
    
    expect_identical(
      exact_fit$ld_representation,
      "exact_blocks"
    )
    
    expect_identical(
      lowrank_fit$ld_representation,
      "lowrank"
    )
    
    expect_identical(
      lowrank_fit$delta,
      exact_fit$delta
    )
    
    expect_identical(
      lowrank_fit$Delta,
      exact_fit$Delta
    )
    
    expect_identical(
      lowrank_fit$gamma,
      exact_fit$gamma
    )
    
    expect_equal(
      lowrank_fit$beta,
      exact_fit$beta,
      tolerance = 1e-12
    )
    
    expect_equal(
      lowrank_fit$par,
      exact_fit$par,
      tolerance = 1e-12
    )
    
    expect_equal(
      lowrank_fit$pip_snp,
      exact_fit$pip_snp,
      tolerance = 1e-12
    )
  }
)


test_that(
  "low-rank mode rejects joint updates until implemented",
  {
    U <- diag(3)
    
    block <- make_lowrank_ld_block(
      U = U,
      lambda = rep(1, 3),
      block_indices = 1:3,
      target_diagonal = 1,
      variance_explained = 1
    )
    
    expect_error(
      bayesmap_summary_blocks(
        bhat = numeric(3),
        LD_blocks = NULL,
        block_indices = NULL,
        N = 1000,
        L = matrix(1, 3, 1),
        A = matrix(1, 3, 1),
        B = matrix(1, 1, 1),
        niter = 2,
        burnin = 0,
        joint_blocks_per_iter = 1L,
        lowrank_blocks = list(block),
        verbose = FALSE
      ),
      "not yet implemented"
    )
  }
)