#' Fit joint-indicator BayesMAP using GWAS summary statistics
#'
#' A convenience wrapper around [bayesmap_summary()] that activates local
#' blocked Gibbs updates for connected SNP, gene and cell indicator triplets.
#' The original summary sampler remains available by calling
#' `bayesmap_summary()` with its default `joint_blocks_per_iter = 0`.
#'
#' @param ... Arguments passed to [bayesmap_summary()].
#' @param joint_blocks_per_iter Number of connected triplets jointly updated
#'   per MCMC iteration. Small values such as 1--5 are recommended initially.
#'
#' @return A `BayesMAPSummary` fit with class `BayesMAPJointSummary` added.
#' @export
bayesmap_joint_summary <- function(..., joint_blocks_per_iter = 5L) {
  fit <- bayesmap_summary(
    ...,
    joint_blocks_per_iter = joint_blocks_per_iter
  )

  class(fit) <- c("BayesMAPJointSummary", class(fit))
  fit
}
