#' Fit joint-indicator BayesMAP using block-wise summary statistics
#'
#' A convenience wrapper around [bayesmap_summary_blocks()] that activates
#' local joint updates for connected SNP, gene and cell indicator triplets.
#' The block-specific LD matrices are used for SNP-effect likelihoods, while
#' the complete annotation hierarchy is retained for gene and cell updates.
#'
#' @param ... Arguments passed to [bayesmap_summary_blocks()].
#' @param joint_blocks_per_iter Number of connected SNP-gene-cell triplets
#'   jointly updated per MCMC iteration. Small values such as 1--5 are
#'   recommended initially.
#'
#' @return A `BayesMAPSummaryBlocks` fit with class
#'   `BayesMAPJointSummaryBlocks` added.
#'
#' @seealso
#' [bayesmap_summary_blocks()],
#' [bayesmap_joint_summary()]
#'
#' @export
bayesmap_joint_summary_blocks <- function(
    ...,
    joint_blocks_per_iter = 5L
) {
  fit <- bayesmap_summary_blocks(
    ...,
    joint_blocks_per_iter =
      joint_blocks_per_iter
  )

  class(fit) <- c(
    "BayesMAPJointSummaryBlocks",
    class(fit)
  )

  fit
}
