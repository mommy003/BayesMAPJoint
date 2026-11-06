#' Construct a diagonal-corrected low-rank LD block
#'
#' Constructs the representation
#'
#' \deqn{
#' R \approx U \Lambda U^\mathsf{T} + D,
#' }
#'
#' where `D` is diagonal and restores the requested LD diagonal.
#'
#' @param U Numeric matrix containing retained eigenvectors in columns.
#' @param lambda Positive retained eigenvalues.
#' @param block_indices Genome-wide one-based SNP indices for the block.
#' @param target_diagonal Desired diagonal of the reconstructed LD block.
#'   For an LD correlation matrix this is normally a vector of ones.
#' @param variance_explained Optional cumulative LD variance explained by
#'   the retained eigenvalues.
#' @param tolerance Numerical tolerance used for validation.
#'
#' @return An object of class `BayesMAPLowRankLDBlock`.
#'
#' @export
make_lowrank_ld_block <- function(
    U,
    lambda,
    block_indices,
    target_diagonal = 1,
    variance_explained = NA_real_,
    tolerance = 1e-8
) {
  if (!is.matrix(U) || !is.numeric(U)) {
    stop("U must be an ordinary numeric matrix.")
  }
  
  if (nrow(U) < 1L || ncol(U) < 1L) {
    stop("U must have at least one row and one column.")
  }
  
  if (any(!is.finite(U))) {
    stop("U cannot contain non-finite values.")
  }
  
  lambda <- as.numeric(lambda)
  
  if (length(lambda) != ncol(U)) {
    stop("length(lambda) must equal ncol(U).")
  }
  
  if (any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("All retained eigenvalues must be positive and finite.")
  }
  
  block_indices <- as.integer(block_indices)
  
  if (length(block_indices) != nrow(U)) {
    stop("length(block_indices) must equal nrow(U).")
  }
  
  if (
    anyNA(block_indices) ||
    any(block_indices < 1L) ||
    anyDuplicated(block_indices)
  ) {
    stop(
      paste0(
        "block_indices must contain unique, positive, ",
        "one-based integer indices."
      )
    )
  }
  
  if (
    length(tolerance) != 1L ||
    !is.finite(tolerance) ||
    tolerance <= 0
  ) {
    stop("tolerance must be a positive finite scalar.")
  }
  
  if (length(target_diagonal) == 1L) {
    target_diagonal <- rep(
      as.numeric(target_diagonal),
      nrow(U)
    )
  } else {
    target_diagonal <- as.numeric(target_diagonal)
  }
  
  if (length(target_diagonal) != nrow(U)) {
    stop(
      "target_diagonal must have length one or nrow(U)."
    )
  }
  
  if (
    any(!is.finite(target_diagonal)) ||
    any(target_diagonal <= 0)
  ) {
    stop(
      "target_diagonal must contain positive finite values."
    )
  }
  
  if (
    length(variance_explained) != 1L ||
    (!is.na(variance_explained) &&
     (!is.finite(variance_explained) ||
      variance_explained <= 0 ||
      variance_explained > 1))
  ) {
    stop(
      "variance_explained must be NA or in the interval (0, 1]."
    )
  }
  
  orthogonality_error <- max(
    abs(
      crossprod(U) -
        diag(ncol(U))
    )
  )
  
  if (
    !is.finite(orthogonality_error) ||
    orthogonality_error > 1e-6
  ) {
    stop(
      paste0(
        "The columns of U must be orthonormal. ",
        "Maximum orthogonality error = ",
        signif(orthogonality_error, 5),
        "."
      )
    )
  }
  
  lowrank_diagonal <- as.numeric(
    (U^2) %*% lambda
  )
  
  diagonal_correction <-
    target_diagonal -
    lowrank_diagonal
  
  substantially_negative <-
    diagonal_correction < -tolerance
  
  if (any(substantially_negative)) {
    stop(
      paste0(
        "The low-rank approximation exceeds the requested ",
        "target diagonal for at least one SNP. ",
        "Minimum diagonal correction = ",
        signif(min(diagonal_correction), 6),
        "."
      )
    )
  }
  
  # Remove only negligible floating-point negative values.
  diagonal_correction[
    diagonal_correction < 0
  ] <- 0
  
  result <- list(
    U = U,
    lambda = lambda,
    diagonal_correction =
      diagonal_correction,
    block_indices = block_indices,
    target_diagonal =
      target_diagonal,
    variance_explained =
      as.numeric(variance_explained),
    full_dimension = nrow(U),
    retained_rank = ncol(U),
    rank_fraction =
      ncol(U) / nrow(U),
    orthogonality_error =
      orthogonality_error
  )
  
  class(result) <- c(
    "BayesMAPLowRankLDBlock",
    "list"
  )
  
  result
}


# Validate one diagonal-corrected low-rank LD block
validate_lowrank_ld_block <- function(
    block,
    tolerance = 1e-8
) {
  if (
    !is.list(block) ||
    !inherits(
      block,
      "BayesMAPLowRankLDBlock"
    )
  ) {
    stop(
      "block must be a BayesMAPLowRankLDBlock object."
    )
  }
  
  required_components <- c(
    "U",
    "lambda",
    "diagonal_correction",
    "block_indices",
    "target_diagonal",
    "full_dimension",
    "retained_rank"
  )
  
  missing_components <- setdiff(
    required_components,
    names(block)
  )
  
  if (length(missing_components) > 0L) {
    stop(
      "Missing low-rank block components: ",
      paste(missing_components, collapse = ", "),
      "."
    )
  }
  
  n <- nrow(block$U)
  k <- ncol(block$U)
  
  if (
    length(block$lambda) != k ||
    length(block$diagonal_correction) != n ||
    length(block$block_indices) != n ||
    length(block$target_diagonal) != n
  ) {
    stop(
      "The low-rank block components have incompatible dimensions."
    )
  }
  
  reconstructed_diagonal <-
    as.numeric(
      (block$U^2) %*% block$lambda
    ) +
    block$diagonal_correction
  
  diagonal_error <- max(
    abs(
      reconstructed_diagonal -
        block$target_diagonal
    )
  )
  
  if (
    !is.finite(diagonal_error) ||
    diagonal_error > tolerance
  ) {
    stop(
      paste0(
        "The reconstructed LD diagonal is invalid. ",
        "Maximum diagonal error = ",
        signif(diagonal_error, 6),
        "."
      )
    )
  }
  
  invisible(TRUE)
}


# Calculate (U Lambda U' + D) %*% vector without constructing dense LD
lowrank_ld_matvec <- function(
    block,
    vector
) {
  validate_lowrank_ld_block(block)
  
  vector <- as.numeric(vector)
  
  if (length(vector) != nrow(block$U)) {
    stop("vector must have length nrow(block$U).")
  }
  
  projected_vector <- as.numeric(
    crossprod(
      block$U,
      vector
    )
  )
  
  as.numeric(
    block$U %*% (
      block$lambda *
        projected_vector
    ) +
      block$diagonal_correction *
      vector
  )
}


# Reconstruct dense LD for testing only
reconstruct_lowrank_ld_block <- function(
    block
) {
  validate_lowrank_ld_block(block)
  
  weighted_U <- sweep(
    block$U,
    MARGIN = 2L,
    STATS = sqrt(block$lambda),
    FUN = "*"
  )
  
  reconstructed_LD <-
    tcrossprod(weighted_U)
  
  diag(reconstructed_LD) <-
    diag(reconstructed_LD) +
    block$diagonal_correction
  
  reconstructed_LD
}



# Validate a list of diagonal-corrected low-rank LD blocks
validate_summary_lowrank_blocks <- function(
    lowrank_blocks,
    m
) {
  if (
    !is.list(lowrank_blocks) ||
    length(lowrank_blocks) < 1L
  ) {
    stop(
      "lowrank_blocks must be a non-empty list."
    )
  }
  
  if (
    length(m) != 1L ||
    !is.finite(m) ||
    m < 1L ||
    m != as.integer(m)
  ) {
    stop("m must be a positive integer.")
  }
  
  m <- as.integer(m)
  
  for (block_number in seq_along(
    lowrank_blocks
  )) {
    tryCatch(
      validate_lowrank_ld_block(
        lowrank_blocks[[block_number]]
      ),
      error = function(error) {
        stop(
          "Invalid low-rank LD block ",
          block_number,
          ": ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
  }
  
  all_indices <- unlist(
    lapply(
      lowrank_blocks,
      function(block) {
        block$block_indices
      }
    ),
    use.names = FALSE
  )
  
  if (length(all_indices) != m) {
    stop(
      "The low-rank LD blocks must contain every SNP exactly once."
    )
  }
  
  if (anyDuplicated(all_indices)) {
    stop(
      "A SNP cannot occur in more than one low-rank LD block."
    )
  }
  
  if (
    !setequal(
      all_indices,
      seq_len(m)
    )
  ) {
    stop(
      paste0(
        "Low-rank block indices must collectively ",
        "equal seq_len(m)."
      )
    )
  }
  
  invisible(TRUE)
}


# Update beta and delta over a list of low-rank LD blocks
updateBetaDeltaSummaryLowRankBlocks <- function(
    bhat,
    lowrank_blocks,
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
  
  if (length(beta) != m) {
    stop("beta must have length length(bhat).")
  }
  
  if (length(delta) != m) {
    stop("delta must have length length(bhat).")
  }
  
  if (length(pi_j) != m) {
    stop("pi_j must have length length(bhat).")
  }
  
  if (
    !is.list(lowrank_blocks) ||
    length(lowrank_blocks) < 1L
  ) {
    stop(
      "lowrank_blocks must be a non-empty list."
    )
  }
  
  rcorr <- numeric(m)
  nnz <- 0L
  
  for (block_number in seq_along(
    lowrank_blocks
  )) {
    block <-
      lowrank_blocks[[block_number]]
    
    index <- block$block_indices
    
    block_update <-
      updateBetaDeltaSummaryLowRank_cpp(
        bhat = bhat[index],
        U = block$U,
        lambda = block$lambda,
        diagonal_correction =
          block$diagonal_correction,
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


# Calculate block-diagonal low-rank LD times a vector
#
# The low-rank blocks must already have been validated by
# validate_summary_lowrank_blocks().
summary_lowrank_block_matvec <- function(
    vector,
    lowrank_blocks
) {
  vector <- as.numeric(vector)
  result <- numeric(length(vector))
  
  for (block_number in seq_along(
    lowrank_blocks
  )) {
    block <-
      lowrank_blocks[[block_number]]
    
    index <- block$block_indices
    block_vector <- vector[index]
    
    projected_vector <- as.numeric(
      crossprod(
        block$U,
        block_vector
      )
    )
    
    result[index] <- as.numeric(
      block$U %*% (
        block$lambda *
          projected_vector
      ) +
        block$diagonal_correction *
        block_vector
    )
  }
  
  result
}