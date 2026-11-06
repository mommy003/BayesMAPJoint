test_that(
  "diagonal-corrected low-rank LD multiplication is correct",
  {
    set.seed(20261106)
    
    random_matrix <- matrix(
      rnorm(40),
      nrow = 10,
      ncol = 4
    )
    
    U <- qr.Q(
      qr(random_matrix)
    )
    
    lambda <- c(
      0.8,
      0.5,
      0.3,
      0.1
    )
    
    block <- make_lowrank_ld_block(
      U = U,
      lambda = lambda,
      block_indices = seq_len(10),
      target_diagonal = 1,
      variance_explained = 0.90
    )
    
    expect_s3_class(
      block,
      "BayesMAPLowRankLDBlock"
    )
    
    expect_true(
      validate_lowrank_ld_block(block)
    )
    
    reconstructed_LD <-
      reconstruct_lowrank_ld_block(
        block
      )
    
    expect_equal(
      diag(reconstructed_LD),
      rep(1, 10),
      tolerance = 1e-12
    )
    
    expect_equal(
      reconstructed_LD,
      t(reconstructed_LD),
      tolerance = 1e-12
    )
    
    test_vector <- seq_len(10) / 10
    
    expect_equal(
      lowrank_ld_matvec(
        block,
        test_vector
      ),
      as.numeric(
        reconstructed_LD %*%
          test_vector
      ),
      tolerance = 1e-12
    )
    
    expect_equal(
      block$full_dimension,
      10L
    )
    
    expect_equal(
      block$retained_rank,
      4L
    )
    
    expect_equal(
      block$rank_fraction,
      0.4
    )
  }
)


test_that(
  "low-rank LD constructor rejects malformed inputs",
  {
    U <- diag(4)[, 1:2, drop = FALSE]
    
    expect_error(
      make_lowrank_ld_block(
        U = U,
        lambda = 1,
        block_indices = seq_len(4)
      ),
      "length\\(lambda\\)"
    )
    
    expect_error(
      make_lowrank_ld_block(
        U = U,
        lambda = c(1, -0.5),
        block_indices = seq_len(4)
      ),
      "positive and finite"
    )
    
    expect_error(
      make_lowrank_ld_block(
        U = U,
        lambda = c(0.5, 0.25),
        block_indices = c(1, 2, 2, 4)
      ),
      "unique"
    )
    
    nonorthogonal_U <- matrix(
      c(
        1, 0, 0, 0,
        1, 0, 0, 0
      ),
      nrow = 4,
      ncol = 2
    )
    
    expect_error(
      make_lowrank_ld_block(
        U = nonorthogonal_U,
        lambda = c(0.5, 0.25),
        block_indices = seq_len(4)
      ),
      "orthonormal"
    )
  }
)