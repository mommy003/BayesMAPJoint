#' Fit BayesMAP using low-rank LD summary statistics
#'
#' Fits the BayesMAP summary-statistics model using block-specific,
#' diagonal-corrected low-rank LD representations:
#'
#' \deqn{
#' R_b \approx U_b \Lambda_b U_b^\mathsf{T} + D_b.
#' }
#'
#' Each low-rank block should be constructed using
#' [make_lowrank_ld_block()].
#'
#' @inheritParams bayesmap_summary_blocks
#'
#' @param lowrank_blocks A non-empty list of
#'   `BayesMAPLowRankLDBlock` objects covering every SNP exactly once.
#' @param ... Additional arguments passed to [bayesmap_summary_blocks()],
#'   including `joint_blocks_per_iter`.
#'
#' @return A fitted object with class `BayesMAPSummaryLowRank`,
#'   `BayesMAPSummaryBlocks`, `BayesMAPSummary`, and `BayesMAP`.
#'
#' @seealso
#' [make_lowrank_ld_block()],
#' [bayesmap_summary_blocks()]
#'
#' @export
bayesmap_summary_lowrank <- function(
    bhat,
    lowrank_blocks,
    N,
    L,
    A,
    B,
    ...
) {
  fit <- bayesmap_summary_blocks(
    bhat = bhat,
    LD_blocks = NULL,
    block_indices = NULL,
    N = N,
    L = L,
    A = A,
    B = B,
    lowrank_blocks = lowrank_blocks,
    ...
  )
  
  class(fit) <- c(
    "BayesMAPSummaryLowRank",
    class(fit)
  )
  
  fit
}