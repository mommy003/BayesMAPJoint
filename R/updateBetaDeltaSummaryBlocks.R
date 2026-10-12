updateBetaDeltaSummaryBlocks <- function(
    bhat,
    LD_blocks,
    block_indices,
    N,
    beta,
    delta,
    pi_j,
    vare,
    sigmaBetaSq
) {
  bhat <- as.numeric(bhat)
  beta <- as.numeric(beta) + 0
  delta <- as.integer(delta) + 0L
  pi_j <- as.numeric(pi_j)

  m <- length(bhat)

  if (!is.list(LD_blocks) || length(LD_blocks) < 1L) {
    stop("LD_blocks must be a non-empty list of numeric matrices.")
  }

  if (!is.list(block_indices)) {
    stop("block_indices must be a list.")
  }

  if (length(LD_blocks) != length(block_indices)) {
    stop("LD_blocks and block_indices must have the same length.")
  }

  if (length(beta) != m) {
    stop("beta must have length length(bhat).")
  }

  if (length(delta) != m) {
    stop("delta must have length length(bhat).")
  }

  if (length(pi_j) != m) {
    stop("pi_j must have length length(bhat).")
  }

  if (length(N) != 1L || !is.finite(N) || N <= 0) {
    stop("N must be a positive finite scalar.")
  }

  if (!is.finite(vare) || vare <= 0) {
    stop("vare must be positive and finite.")
  }

  if (!is.finite(sigmaBetaSq) || sigmaBetaSq <= 0) {
    stop("sigmaBetaSq must be positive and finite.")
  }

  all_indices <- unlist(
    block_indices,
    use.names = FALSE
  )

  if (length(all_indices) != m) {
    stop(
      "The block indices must contain every SNP exactly once."
    )
  }

  if (!setequal(all_indices, seq_len(m))) {
    stop(
      "block_indices must collectively equal seq_len(length(bhat))."
    )
  }

  if (anyDuplicated(all_indices)) {
    stop("A SNP cannot occur in more than one LD block.")
  }

  rcorr <- numeric(m)
  nnz <- 0L

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

    block_update <- updateBetaDeltaSummary_cpp(
      bhat = bhat[index],
      LD = LD_block,
      N = N,
      beta = beta[index],
      delta = delta[index],
      pi_j = pi_j[index],
      vare = vare,
      sigmaBetaSq = sigmaBetaSq
    )

    beta[index] <- block_update$beta
    delta[index] <- block_update$delta
    rcorr[index] <- block_update$rcorr

    nnz <- nnz + block_update$nnz
  }

  list(
    beta = beta,
    delta = delta,
    rcorr = rcorr,
    nnz = nnz
  )
}
