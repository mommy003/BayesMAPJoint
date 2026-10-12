make_ld_blocks <- function(
    chr,
    position,
    block_table
) {
  chr <- as.character(chr)
  position <- as.numeric(position)

  m <- length(chr)

  if (length(position) != m) {
    stop("chr and position must have the same length.")
  }

  required_columns <- c(
    "chr",
    "start",
    "end"
  )

  if (!all(required_columns %in% names(block_table))) {
    stop(
      "block_table must contain columns: chr, start, end."
    )
  }

  block_indices <- vector(
    "list",
    nrow(block_table)
  )

  for (b in seq_len(nrow(block_table))) {
    block_indices[[b]] <- which(
      chr == as.character(block_table$chr[b]) &
        position >= block_table$start[b] &
        position <= block_table$end[b]
    )
  }

  keep <- lengths(block_indices) > 0L
  block_indices <- block_indices[keep]
  block_table <- block_table[keep, , drop = FALSE]

  all_indices <- unlist(
    block_indices,
    use.names = FALSE
  )

  if (anyDuplicated(all_indices)) {
    stop(
      "Some SNPs were assigned to more than one block."
    )
  }

  unassigned <- setdiff(
    seq_len(m),
    all_indices
  )

  list(
    block_indices = block_indices,
    block_table = block_table,
    unassigned = unassigned
  )
}



####Step 2: create LD submatrices
extract_ld_blocks <- function(
    LD,
    block_indices
) {
  if (!is.matrix(LD) || nrow(LD) != ncol(LD)) {
    stop("LD must be a square matrix.")
  }

  lapply(
    block_indices,
    function(index) {
      LD[index, index, drop = FALSE]
    }
  )
}


#####
assess_ld_blocks <- function(
    LD,
    block_indices
) {
  if (!is.matrix(LD) ||
      !is.numeric(LD) ||
      nrow(LD) != ncol(LD)) {
    stop("LD must be a square numeric matrix.")
  }

  m <- nrow(LD)

  if (!is.list(block_indices) ||
      length(block_indices) < 1L) {
    stop("block_indices must be a non-empty list.")
  }

  all_indices <- unlist(
    block_indices,
    use.names = FALSE
  )

  if (length(all_indices) != m ||
      !setequal(all_indices, seq_len(m))) {
    stop("block_indices must contain every SNP exactly once.")
  }

  if (anyDuplicated(all_indices)) {
    stop("A SNP cannot occur in more than one block.")
  }

  block_id <- integer(m)

  for (b in seq_along(block_indices)) {
    index <- as.integer(block_indices[[b]])

    if (any(index < 1L | index > m)) {
      stop("A block index is outside the valid SNP range.")
    }

    block_id[index] <- b
  }

  cross_mask <- outer(
    block_id,
    block_id,
    FUN = "!="
  )

  cross_values <- LD[cross_mask]

  LD_approx <- matrix(
    0,
    nrow = m,
    ncol = m
  )

  for (b in seq_along(block_indices)) {
    index <- block_indices[[b]]

    LD_approx[index, index] <-
      LD[index, index, drop = FALSE]
  }

  if (length(cross_values) == 0L) {
    mean_abs_cross_r <- 0
    max_abs_cross_r <- 0
    mean_cross_r2 <- 0
    prop_abs_r_gt_005 <- 0
    prop_abs_r_gt_010 <- 0
  } else {
    mean_abs_cross_r <- mean(abs(cross_values))
    max_abs_cross_r <- max(abs(cross_values))
    mean_cross_r2 <- mean(cross_values^2)

    prop_abs_r_gt_005 <-
      mean(abs(cross_values) > 0.05)

    prop_abs_r_gt_010 <-
      mean(abs(cross_values) > 0.10)
  }

  relative_frobenius_error <-
    sqrt(sum((LD - LD_approx)^2)) /
    sqrt(sum(LD^2))

  c(
    mean_abs_cross_r = mean_abs_cross_r,
    max_abs_cross_r = max_abs_cross_r,
    mean_cross_r2 = mean_cross_r2,
    prop_abs_r_gt_005 = prop_abs_r_gt_005,
    prop_abs_r_gt_010 = prop_abs_r_gt_010,
    relative_frobenius_error =
      relative_frobenius_error
  )
}





###Step 3: add a diagnostic function
check_ld_block_quality <- function(
    LD,
    block_indices,
    frobenius_warning = 0.10,
    max_cross_r_warning = 0.10
) {
  diagnostics <- assess_ld_blocks(
    LD = LD,
    block_indices = block_indices
  )

  if (
    diagnostics["relative_frobenius_error"] >
    frobenius_warning
  ) {
    warning(
      "The block approximation discards substantial LD: ",
      "relative Frobenius error = ",
      round(
        diagnostics["relative_frobenius_error"],
        3
      ),
      "."
    )
  }

  if (
    diagnostics["max_abs_cross_r"] >
    max_cross_r_warning
  ) {
    warning(
      "Strong cross-block LD remains: max |r| = ",
      round(
        diagnostics["max_abs_cross_r"],
        3
      ),
      "."
    )
  }

  diagnostics
}






