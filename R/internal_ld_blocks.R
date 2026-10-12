# Validate a list of LD blocks and their SNP indices
validate_summary_ld_blocks <- function(
    LD_blocks,
    block_indices,
    m
) {
  if (!is.list(LD_blocks) || length(LD_blocks) < 1L) {
    stop("LD_blocks must be a non-empty list of numeric matrices.")
  }

  if (!is.list(block_indices)) {
    stop("block_indices must be a list.")
  }

  if (length(LD_blocks) != length(block_indices)) {
    stop("LD_blocks and block_indices must have the same length.")
  }

  all_indices <- unlist(
    block_indices,
    use.names = FALSE
  )

  if (length(all_indices) != m) {
    stop("The LD blocks must contain every SNP exactly once.")
  }

  if (anyDuplicated(all_indices)) {
    stop("A SNP cannot occur in more than one LD block.")
  }

  if (!setequal(all_indices, seq_len(m))) {
    stop("block_indices must collectively equal seq_len(length(bhat)).")
  }

  for (block in seq_along(LD_blocks)) {
    index <- as.integer(block_indices[[block]])
    LD_block <- LD_blocks[[block]]

    if (length(index) < 1L) {
      stop("LD blocks cannot be empty.")
    }

    if (any(index < 1L | index > m)) {
      stop("A block index is outside the valid SNP range.")
    }

    if (!is.matrix(LD_block) || !is.numeric(LD_block)) {
      stop("Each LD block must be an ordinary numeric matrix.")
    }

    block_size <- length(index)

    if (any(dim(LD_block) != c(block_size, block_size))) {
      stop(
        "Each LD block dimension must match its block-index length."
      )
    }

    if (any(!is.finite(LD_block))) {
      stop("LD blocks cannot contain non-finite values.")
    }

    if (any(diag(LD_block) <= 0)) {
      stop("All LD-block diagonal elements must be positive.")
    }

    if (max(abs(LD_block - t(LD_block))) > 1e-8) {
      warning(
        "LD block ", block,
        " is not exactly symmetric. ",
        "It will be replaced by (LD + t(LD)) / 2."
      )

      LD_blocks[[block]] <-
        (LD_block + t(LD_block)) / 2
    }
  }

  list(
    LD_blocks = LD_blocks,
    block_indices = lapply(
      block_indices,
      as.integer
    )
  )
}


# Calculate the block-diagonal product LD %*% vector
summary_ld_block_matvec <- function(
    vector,
    LD_blocks,
    block_indices
) {
  vector <- as.numeric(vector)
  result <- numeric(length(vector))

  for (block in seq_along(LD_blocks)) {
    index <- block_indices[[block]]

    result[index] <- as.numeric(
      LD_blocks[[block]] %*% vector[index]
    )
  }

  result
}
