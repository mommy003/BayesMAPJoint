test_that("make_ld_blocks assigns SNPs correctly", {
  chr <- c("1", "1", "1", "2", "2", "2")
  position <- c(100, 200, 300, 100, 200, 300)

  block_table <- data.frame(
    chr = c("1", "2"),
    start = c(1, 1),
    end = c(1000, 1000)
  )

  out <- make_ld_blocks(
    chr = chr,
    position = position,
    block_table = block_table
  )

  expect_equal(
    out$block_indices,
    list(
      1:3,
      4:6
    )
  )

  expect_length(
    out$unassigned,
    0L
  )

  expect_equal(
    nrow(out$block_table),
    2L
  )
})


test_that("make_ld_blocks reports unassigned SNPs", {
  chr <- c("1", "1", "1", "2")
  position <- c(100, 200, 300, 100)

  block_table <- data.frame(
    chr = "1",
    start = 1,
    end = 250
  )

  out <- make_ld_blocks(
    chr = chr,
    position = position,
    block_table = block_table
  )

  expect_equal(
    out$block_indices,
    list(1:2)
  )

  expect_equal(
    out$unassigned,
    c(3L, 4L)
  )
})


test_that("make_ld_blocks rejects overlapping blocks", {
  chr <- rep("1", 4)
  position <- c(100, 200, 300, 400)

  block_table <- data.frame(
    chr = c("1", "1"),
    start = c(1, 150),
    end = c(250, 350)
  )

  expect_error(
    make_ld_blocks(
      chr = chr,
      position = position,
      block_table = block_table
    ),
    "Some SNPs were assigned to more than one block"
  )
})


test_that("extract_ld_blocks returns correct dimensions", {
  LD <- matrix(
    seq_len(36),
    nrow = 6,
    ncol = 6
  )

  block_indices <- list(
    1:2,
    3:6
  )

  blocks <- extract_ld_blocks(
    LD = LD,
    block_indices = block_indices
  )

  expect_equal(
    length(blocks),
    2L
  )

  expect_equal(
    dim(blocks[[1]]),
    c(2L, 2L)
  )

  expect_equal(
    dim(blocks[[2]]),
    c(4L, 4L)
  )

  expect_equal(
    blocks[[1]],
    LD[1:2, 1:2, drop = FALSE]
  )

  expect_equal(
    blocks[[2]],
    LD[3:6, 3:6, drop = FALSE]
  )
})


test_that("assess_ld_blocks is exact for block-diagonal LD", {
  LD1 <- matrix(
    c(
      1.0, 0.2,
      0.2, 1.0
    ),
    nrow = 2,
    byrow = TRUE
  )

  LD2 <- matrix(
    c(
      1.0, 0.1,
      0.1, 1.0
    ),
    nrow = 2,
    byrow = TRUE
  )

  LD <- matrix(
    0,
    nrow = 4,
    ncol = 4
  )

  LD[1:2, 1:2] <- LD1
  LD[3:4, 3:4] <- LD2

  diagnostics <- assess_ld_blocks(
    LD = LD,
    block_indices = list(
      1:2,
      3:4
    )
  )

  expect_equal(
    unname(diagnostics["mean_abs_cross_r"]),
    0
  )

  expect_equal(
    unname(diagnostics["max_abs_cross_r"]),
    0
  )

  expect_equal(
    unname(diagnostics["relative_frobenius_error"]),
    0
  )
})


test_that("check_ld_block_quality warns for poor blocks", {
  LD <- matrix(
    c(
      1.0, 0.2, 0.3, 0.2,
      0.2, 1.0, 0.25, 0.15,
      0.3, 0.25, 1.0, 0.2,
      0.2, 0.15, 0.2, 1.0
    ),
    nrow = 4,
    byrow = TRUE
  )

  block_indices <- list(
    1:2,
    3:4
  )

  expect_warning(
    expect_warning(
      check_ld_block_quality(
        LD = LD,
        block_indices = block_indices,
        frobenius_warning = 0.01,
        max_cross_r_warning = 0.10
      ),
      "Strong cross-block LD remains"
    ),
    "block approximation discards substantial LD"
  )
})


test_that("check_ld_block_quality is silent for good blocks", {
  LD <- diag(4)

  expect_silent(
    check_ld_block_quality(
      LD = LD,
      block_indices = list(
        1:2,
        3:4
      )
    )
  )
})
