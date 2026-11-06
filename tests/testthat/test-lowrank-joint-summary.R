test_that(
  "full-rank joint kernel reproduces exact LD joint kernel",
  {
    set.seed(20261106)
    
    m <- 5L
    G <- 2L
    C <- 2L
    
    LD <- diag(m)
    
    block <- make_lowrank_ld_block(
      U = diag(m),
      lambda = rep(1, m),
      block_indices = seq_len(m),
      target_diagonal = 1,
      variance_explained = 1
    )
    
    L <- matrix(
      c(
        1, 0,
        1, 0,
        1, 0,
        0, 1,
        0, 1
      ),
      nrow = m,
      ncol = G,
      byrow = TRUE
    )
    
    A <- matrix(
      c(
        1, 0,
        1, 0,
        0, 1,
        0, 1,
        1, 0
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
    
    make_common_arguments <- function() {
      list(
        N = 5000,
        beta = c(0.01, 0, -0.015, 0, 0.02),
        delta = c(1L, 0L, 1L, 0L, 1L),
        Delta = c(1L, 0L),
        gamma = c(1L, 0L),
        L = L,
        A = A,
        B = B,
        block_indices = seq_len(m),
        snp_index = 2L,
        gene_index = 1L,
        cell_index = 1L,
        vare = 0.8,
        sigmaBetaSq = 0.001,
        alpha1 = 0.2,
        alpha2 = -0.1,
        alpha3 = 0.15,
        mu_pi = qnorm(0.1),
        mu_Pi = qnorm(0.1),
        baseline = numeric(m),
        pathways = numeric(G),
        rho = 0.1
      )
    }
    
    set.seed(24680)
    
    exact_update <- do.call(
      updateJointIndicatorsSummaryBlocks_cpp,
      c(
        list(
          bhat_block = bhat,
          LD_block = LD
        ),
        make_common_arguments()
      )
    )
    
    set.seed(24680)
    
    lowrank_update <- do.call(
      updateJointIndicatorsSummaryLowRankBlocks_cpp,
      c(
        list(
          bhat_block = bhat,
          U = block$U,
          lambda = block$lambda,
          diagonal_correction =
            block$diagonal_correction
        ),
        make_common_arguments()
      )
    )
    
    expect_identical(
      lowrank_update$delta,
      exact_update$delta
    )
    
    expect_identical(
      lowrank_update$Delta,
      exact_update$Delta
    )
    
    expect_identical(
      lowrank_update$gamma,
      exact_update$gamma
    )
    
    expect_identical(
      lowrank_update$sampled_states,
      exact_update$sampled_states
    )
    
    expect_equal(
      lowrank_update$beta,
      exact_update$beta,
      tolerance = 1e-12
    )
    
    expect_equal(
      lowrank_update$rcorr_block,
      exact_update$rcorr_block,
      tolerance = 1e-12
    )
    
    expect_equal(
      lowrank_update$selected_probability,
      exact_update$selected_probability,
      tolerance = 1e-12
    )
  }
)