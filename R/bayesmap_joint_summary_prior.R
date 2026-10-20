#' Fit sparse-prior joint BayesMAP using GWAS summary statistics
#'
#' Runs the summary-statistics BayesMAP sampler with local joint updates for
#' connected SNP, gene and cell indicators and the same sparse higher-level
#' priors used by [bayesmap_joint_prior()].
#'
#' @param ... Arguments passed to [bayesmap_summary()].
#' @param joint_blocks_per_iter Number of connected triplets jointly updated
#'   per MCMC iteration.
#' @param priorPiGene Prior mean gene inclusion probability.
#' @param priorRho Prior mean cell inclusion probability.
#' @param rhoPriorStrength Concentration of the beta prior on cell inclusion.
#' @param sigmaMuPiSq Prior variance of the gene-level probit intercept.
#'
#' @details The default starting probabilities are 0.01 for SNPs, 0.10 for
#'   genes and 0.10 for cells, with `startH2 = 0.20`. `mu_pi` and `mu_Pi`
#'   are initialized at the corresponding probit values. These are starting
#'   values only; the informative priors are applied to `mu_Pi` and `rho`.
#'
#' @return A `BayesMAPSummary` fit with class
#'   `BayesMAPJointSummaryPrior` added.
#' @export
bayesmap_joint_summary_prior <- function(
    ...,
    joint_blocks_per_iter = 5L,
    priorPiGene = 0.10,
    priorRho = 0.05,
    rhoPriorStrength = 20,
    sigmaMuPiSq = 1
) {
  arguments <- utils::modifyList(
    list(
      startPiSnp = 0.01,
      startPiGene = 0.10,
      startRho = 0.10,
      startH2 = 0.20,
      mu_pi = NULL,
      mu_Pi = NULL
    ),
    list(...)
  )

  arguments <- c(
    arguments,
    list(
    joint_blocks_per_iter = joint_blocks_per_iter,
    priorPiGene = priorPiGene,
    priorRho = priorRho,
    rhoPriorStrength = rhoPriorStrength,
    sigmaMuPiSq = sigmaMuPiSq
    )
  )

  fit <- do.call(bayesmap_summary, arguments)

  class(fit) <- c("BayesMAPJointSummaryPrior", class(fit))
  fit
}
