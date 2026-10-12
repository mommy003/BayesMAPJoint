#' Summarise a fitted BayesMAP model
#'
#' @param object A fitted object returned by [bayesmap()] or
#'   [bayesmap_summary()].
#' @param ... Additional arguments, currently unused.
#'
#' @return A list containing posterior summaries, returned invisibly.
#' @export
summary.BayesMAP <- function(object, ...) {

  if (!inherits(object, "BayesMAP")) {
    stop("`object` must inherit from class 'BayesMAP'.", call. = FALSE)
  }

  par_means <- if (!is.null(object$par)) {
    colMeans(object$par, na.rm = TRUE)
  } else {
    NULL
  }

  cat("\nBayesMAP model summary\n")
  cat("======================\n\n")

  if (!is.null(object$dimensions)) {
    if ("n" %in% names(object$dimensions)) {
      cat("Individuals        :", object$dimensions["n"], "\n")
    }

    cat("SNPs               :", object$dimensions["m"], "\n")
    cat("Genes              :", object$dimensions["G"], "\n")
    cat("Cells/pathways     :", object$dimensions["C"], "\n")
  }

  if (!is.null(object$N)) {
    cat("GWAS sample size   :", object$N, "\n")
  }

  cat("MCMC iterations    :", object$niter, "\n")
  cat("Burn-in            :", object$burnin, "\n")
  cat("Thinning interval  :", object$thin, "\n")
  cat("Retained samples   :", object$nkeep, "\n")

  if (!is.null(par_means)) {
    cat("\nPosterior means\n")
    cat("----------------------\n")

    parameters_to_show <- intersect(
      c(
        "SigmaBetaSq",
        "Vare",
        "Varg",
        "h2",
        "Alpha1",
        "Alpha2",
        "Alpha3",
        "Rho",
        "Nnz",
        "mu_pi",
        "mu_Pi"
      ),
      names(par_means)
    )

    print(
      round(
        par_means[parameters_to_show],
        digits = 4
      )
    )
  }

  cat("\nMean posterior inclusion probabilities\n")
  cat("--------------------------------------\n")

  if (!is.null(object$pip_snp)) {
    cat("SNP PIP            :", round(mean(object$pip_snp), 4), "\n")
  }

  if (!is.null(object$pip_gene)) {
    cat("Gene PIP           :", round(mean(object$pip_gene), 4), "\n")
  }

  if (!is.null(object$pip_cell)) {
    cat("Cell/pathway PIP   :", round(mean(object$pip_cell), 4), "\n")
  }

  out <- list(
    posterior_means = par_means,
    pip_snp = object$pip_snp,
    pip_gene = object$pip_gene,
    pip_cell = object$pip_cell
  )

  invisible(out)
}
