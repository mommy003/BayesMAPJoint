#' Bayesian Multi-layer Association Mapping using GWAS Summary Statistics
#'
#' Fits the BayesMAP hierarchical Bayesian model using GWAS summary statistics
#' and a dense linkage disequilibrium (LD) correlation matrix. The model jointly
#' identifies trait-associated SNPs, genes, and cell/pathway annotations through
#' a hierarchical Bayesian variable selection framework.
#'
#' Posterior inference is performed using Markov chain Monte Carlo (MCMC),
#' while computationally intensive updates are accelerated using Rcpp.
#'
#' @param bhat Vector of marginal SNP effect estimates.
#' @param LD Dense SNP correlation matrix.
#' @param N GWAS sample size.
#' @param L SNP-to-gene annotation matrix.
#' @param A SNP-to-cell/pathway annotation matrix.
#' @param B Gene-to-cell/pathway annotation matrix.
#' @param baseline Optional SNP-level prior covariates.
#' @param pathways Optional gene-level prior covariates.
#' @param phenotype_variance Phenotypic variance (default = 1).
#' @param niter Number of MCMC iterations.
#' @param burnin Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param store_beta Store posterior SNP effects.
#' @param store_delta Store posterior SNP indicators.
#' @param startPiSnp Initial SNP inclusion probability.
#' @param startPiGene Initial gene inclusion probability.
#' @param startRho Initial cell/pathway inclusion probability.
#' @param startH2 Initial SNP heritability.
#' @param mu_pi Initial SNP intercept.
#' @param mu_Pi Initial gene intercept.
#' @param sigmaAlphaSq Prior variance for enrichment parameters.
#' @param nub Prior degrees of freedom for SNP effects.
#' @param nue Prior degrees of freedom for residual variance.
#' @param verbose Logical indicating whether progress is printed.
#'
#' @return
#' An object of class `"BayesMAPSummary"` and `"BayesMAP"` containing
#'
#' * posterior samples of model parameters
#' * posterior SNP effects
#' * SNP posterior inclusion probabilities
#' * gene posterior inclusion probabilities
#' * cell/pathway posterior inclusion probabilities
#' * estimated SNP heritability
#'
#' @export

bayesmap_summary <- function(
    bhat,
    LD,
    N,
    L,
    A,
    B,
    baseline = NULL,
    pathways = NULL,
    phenotype_variance = 1,
    niter = 2000,
    burnin = 0,
    thin = 1,
    store_beta = TRUE,
    store_delta = TRUE,
    startPiSnp = 0.05,
    startPiGene = 0.10,
    startRho = 0.10,
    startH2 = 0.5,
    mu_pi = NULL,
    mu_Pi = NULL,
    sigmaAlphaSq = 1,
    nub = 4,
    nue = 4,
    verbose = TRUE
) {
  ## ================================================================
  ## 1. Dimensions and validation
  ## ================================================================

  bhat <- as.numeric(bhat)

  m <- length(bhat)
  G <- ncol(L)
  C <- ncol(A)

  if (!is.matrix(LD)) {
    stop("LD must currently be an ordinary numeric matrix.")
  }

  if (!is.numeric(LD)) {
    stop("LD must be numeric.")
  }

  if (any(dim(LD) != c(m, m))) {
    stop("LD must be an m by m matrix, where m = length(bhat).")
  }

  if (nrow(L) != m) {
    stop("L must be m by G, where m = length(bhat).")
  }

  if (nrow(A) != m) {
    stop("A must be m by C, where m = length(bhat).")
  }

  if (nrow(B) != G || ncol(B) != C) {
    stop("B must be G by C.")
  }

  if (length(N) != 1L || !is.finite(N) || N <= 1) {
    stop("N must be a finite scalar greater than 1.")
  }

  if (!is.finite(phenotype_variance) ||
      phenotype_variance <= 0) {
    stop("phenotype_variance must be positive and finite.")
  }

  if (any(!is.finite(bhat))) {
    stop("bhat contains non-finite values.")
  }

  if (any(!is.finite(LD))) {
    stop("LD contains non-finite values.")
  }

  if (max(abs(LD - t(LD))) > 1e-8) {
    warning(
      "LD is not exactly symmetric. ",
      "It will be replaced by (LD + t(LD)) / 2."
    )

    LD <- (LD + t(LD)) / 2
  }

  if (any(diag(LD) <= 0)) {
    stop("All diagonal elements of LD must be positive.")
  }

  if (is.null(baseline)) {
    baseline <- numeric(m)
  } else {
    baseline <- as.numeric(baseline)
  }

  if (is.null(pathways)) {
    pathways <- numeric(G)
  } else {
    pathways <- as.numeric(pathways)
  }

  if (length(baseline) != m) {
    stop("baseline must have length m.")
  }

  if (length(pathways) != G) {
    stop("pathways must have length G.")
  }

  if (length(niter) != 1L ||
      !is.finite(niter) ||
      niter < 1 ||
      niter != as.integer(niter)) {
    stop("niter must be a positive integer.")
  }

  niter <- as.integer(niter)

  if (length(burnin) != 1L ||
      !is.finite(burnin) ||
      burnin < 0 ||
      burnin >= niter ||
      burnin != as.integer(burnin)) {
    stop("burnin must be an integer from 0 to niter - 1.")
  }

  burnin <- as.integer(burnin)

  if (length(thin) != 1L ||
      !is.finite(thin) ||
      thin < 1 ||
      thin != as.integer(thin)) {
    stop("thin must be a positive integer.")
  }

  thin <- as.integer(thin)

  if (!is.logical(store_beta) ||
      length(store_beta) != 1L ||
      is.na(store_beta)) {
    stop("store_beta must be TRUE or FALSE.")
  }

  if (!is.logical(store_delta) ||
      length(store_delta) != 1L ||
      is.na(store_delta)) {
    stop("store_delta must be TRUE or FALSE.")
  }

  if (!is.logical(verbose) ||
      length(verbose) != 1L ||
      is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }

  if (!is.finite(startPiSnp) ||
      startPiSnp <= 0 ||
      startPiSnp >= 1) {
    stop("startPiSnp must be strictly between 0 and 1.")
  }

  if (!is.finite(startPiGene) ||
      startPiGene <= 0 ||
      startPiGene >= 1) {
    stop("startPiGene must be strictly between 0 and 1.")
  }

  if (!is.finite(startRho) ||
      startRho <= 0 ||
      startRho >= 1) {
    stop("startRho must be strictly between 0 and 1.")
  }

  if (!is.finite(startH2) ||
      startH2 <= 0 ||
      startH2 >= 1) {
    stop("startH2 must be strictly between 0 and 1.")
  }

  if (!is.finite(sigmaAlphaSq) ||
      sigmaAlphaSq <= 0) {
    stop("sigmaAlphaSq must be positive.")
  }

  if (!is.finite(nub) || nub <= 2) {
    stop("nub must be greater than 2.")
  }

  if (!is.finite(nue) || nue <= 2) {
    stop("nue must be greater than 2.")
  }

  ## ================================================================
  ## 2. Initialise probit intercepts
  ## ================================================================

  if (is.null(mu_pi)) {
    mu_pi <- qnorm(startPiSnp)
  }

  if (is.null(mu_Pi)) {
    mu_Pi <- qnorm(startPiGene)
  }

  if (length(mu_pi) != 1L || !is.finite(mu_pi)) {
    stop("mu_pi must be a finite scalar.")
  }

  if (length(mu_Pi) != 1L || !is.finite(mu_Pi)) {
    stop("mu_Pi must be a finite scalar.")
  }

  ## ================================================================
  ## 3. Burn-in and thinning
  ## ================================================================

  keep_iterations <- seq.int(
    from = burnin + 1L,
    to = niter,
    by = thin
  )

  nkeep <- length(keep_iterations)

  keep_lookup <- integer(niter)
  keep_lookup[keep_iterations] <- seq_len(nkeep)

  ## ================================================================
  ## 4. Initial variance components
  ## ================================================================

  varg <- phenotype_variance * startH2
  vare <- phenotype_variance * (1 - startH2)

  sigmaBetaSq <- varg / m

  scaleb <- ((nub - 2) / nub) * sigmaBetaSq
  scalee <- ((nue - 2) / nue) * vare

  small_positive_value <- 1e-10

  ## ================================================================
  ## 5. Initial state
  ## ================================================================

  beta <- numeric(m)

  delta <- rbinom(
    m,
    size = 1L,
    prob = startPiSnp
  )

  Delta <- rbinom(
    G,
    size = 1L,
    prob = startPiGene
  )

  gamma <- rbinom(
    C,
    size = 1L,
    prob = startRho
  )

  alpha1 <- 0
  alpha2 <- 0
  alpha3 <- 0
  rho <- startRho

  z_pi <- numeric(m)
  z_Pi <- numeric(G)

  ## ================================================================
  ## 6. Storage
  ## ================================================================

  beta_mcmc <- if (store_beta) {
    matrix(0, nrow = nkeep, ncol = m)
  } else {
    NULL
  }

  delta_mcmc <- if (store_delta) {
    matrix(0L, nrow = nkeep, ncol = m)
  } else {
    NULL
  }

  Delta_mcmc <- matrix(
    0L,
    nrow = nkeep,
    ncol = G
  )

  gamma_mcmc <- matrix(
    0L,
    nrow = nkeep,
    ncol = C
  )

  keptIter <- matrix(
    0,
    nrow = nkeep,
    ncol = 11L
  )

  colnames(keptIter) <- c(
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
  )

  pip_snp_sum <- numeric(m)
  pip_gene_sum <- numeric(G)
  pip_cell_sum <- numeric(C)

  ## ================================================================
  ## 7. MCMC
  ## ================================================================

  for (iter in seq_len(niter)) {
    ## --------------------------------------------------------------
    ## 7.1 SNP prior inclusion probabilities
    ## --------------------------------------------------------------

    LDelta <- as.numeric(L %*% Delta)
    Agamma <- as.numeric(A %*% gamma)

    eta_pi <-
      mu_pi +
      LDelta * alpha1 +
      Agamma * alpha2 +
      baseline

    pi_j <- pnorm(eta_pi)

    pi_j <- pmin(
      pmax(pi_j, 1e-12),
      1 - 1e-12
    )

    ## --------------------------------------------------------------
    ## 7.2 Summary-data beta and delta update
    ## --------------------------------------------------------------

    snp_update <- updateBetaDeltaSummary_cpp(
      bhat = bhat,
      LD = LD,
      N = N,
      beta = beta,
      delta = delta,
      pi_j = pi_j,
      vare = vare,
      sigmaBetaSq = sigmaBetaSq
    )

    beta <- snp_update$beta
    delta <- snp_update$delta
    nnz <- snp_update$nnz

    ## --------------------------------------------------------------
    ## 7.3 Latent SNP probit variables
    ## --------------------------------------------------------------

    LDelta <- as.numeric(L %*% Delta)
    Agamma <- as.numeric(A %*% gamma)

    eta_pi <-
      mu_pi +
      LDelta * alpha1 +
      Agamma * alpha2 +
      baseline

    z_pi <- updateZ_cpp(
      indicator = delta,
      eta = eta_pi
    )

    ## --------------------------------------------------------------
    ## 7.4 Gene indicators
    ## --------------------------------------------------------------

    Delta <- updateDelta_cpp(
      Delta = Delta,
      delta = delta,
      L = L,
      gamma = gamma,
      A = A,
      B = B,
      alpha1 = alpha1,
      alpha2 = alpha2,
      alpha3 = alpha3,
      mu_pi = mu_pi,
      mu_Pi = mu_Pi,
      baseline = baseline,
      pathways = pathways,
      rho = rho
    )

    ## --------------------------------------------------------------
    ## 7.5 Latent gene probit variables
    ## --------------------------------------------------------------

    Bgamma <- as.numeric(B %*% gamma)

    xi_g <-
      mu_Pi +
      Bgamma * alpha3 +
      pathways

    z_Pi <- updateZ_cpp(
      indicator = Delta,
      eta = xi_g
    )

    ## --------------------------------------------------------------
    ## 7.6 Cell indicators
    ## --------------------------------------------------------------

    gamma <- updateGamma_cpp(
      gamma = gamma,
      delta = delta,
      Delta = Delta,
      A = A,
      B = B,
      L = L,
      alpha1 = alpha1,
      alpha2 = alpha2,
      alpha3 = alpha3,
      mu_pi = mu_pi,
      mu_Pi = mu_Pi,
      baseline = baseline,
      pathways = pathways,
      rho = rho
    )

    ## --------------------------------------------------------------
    ## 7.7 Cell inclusion probability
    ## --------------------------------------------------------------

    gamma_sum <- sum(gamma)

    rho <- rbeta(
      1L,
      shape1 = gamma_sum + 1,
      shape2 = C - gamma_sum + 1
    )

    ## --------------------------------------------------------------
    ## 7.8 Alpha1 and alpha2
    ## --------------------------------------------------------------

    LDelta <- as.numeric(L %*% Delta)
    Agamma <- as.numeric(A %*% gamma)

    ## Alpha1

    zcorr1 <-
      z_pi -
      mu_pi -
      Agamma * alpha2 -
      baseline

    C1 <- sum(LDelta^2) + 1 / sigmaAlphaSq
    r1 <- sum(LDelta * zcorr1)

    mean1 <- r1 / C1
    sd1 <- sqrt(1 / C1)

    if (!is.finite(mean1)) {
      mean1 <- 0
    }

    if (!is.finite(sd1) || sd1 <= 0) {
      sd1 <- 1
    }

    alpha1 <- rtruncnorm_scalar(
      mean = mean1,
      sd = sd1,
      lower = 0
    )

    ## Alpha2

    zcorr2 <-
      z_pi -
      mu_pi -
      LDelta * alpha1 -
      baseline

    C2 <- sum(Agamma^2) + 1 / sigmaAlphaSq
    r2 <- sum(Agamma * zcorr2)

    mean2 <- r2 / C2
    sd2 <- sqrt(1 / C2)

    if (!is.finite(mean2)) {
      mean2 <- 0
    }

    if (!is.finite(sd2) || sd2 <= 0) {
      sd2 <- 1
    }

    alpha2 <- rtruncnorm_scalar(
      mean = mean2,
      sd = sd2,
      lower = 0
    )

    ## --------------------------------------------------------------
    ## 7.9 Alpha3
    ## --------------------------------------------------------------

    Bgamma <- as.numeric(B %*% gamma)

    zcorr3 <-
      z_Pi -
      mu_Pi -
      pathways

    C3 <- sum(Bgamma^2) + 1 / sigmaAlphaSq
    r3 <- sum(Bgamma * zcorr3)

    mean3 <- r3 / C3
    sd3 <- sqrt(1 / C3)

    if (!is.finite(mean3)) {
      mean3 <- 0
    }

    if (!is.finite(sd3) || sd3 <= 0) {
      sd3 <- 1
    }

    alpha3 <- rtruncnorm_scalar(
      mean = mean3,
      sd = sd3,
      lower = 0
    )

    ## --------------------------------------------------------------
    ## 7.10 Probit intercepts
    ## --------------------------------------------------------------

    zcorr0 <-
      z_pi -
      LDelta * alpha1 -
      Agamma * alpha2 -
      baseline

    C0 <- m + 1 / sigmaAlphaSq
    r0 <- sum(zcorr0)

    mu_pi <- rnorm(
      1L,
      mean = r0 / C0,
      sd = sqrt(1 / C0)
    )

    zcorrP <-
      z_Pi -
      Bgamma * alpha3 -
      pathways

    C0P <- G + 1 / sigmaAlphaSq
    r0P <- sum(zcorrP)

    mu_Pi <- rnorm(
      1L,
      mean = r0P / C0P,
      sd = sqrt(1 / C0P)
    )

    ## --------------------------------------------------------------
    ## 7.11 SNP-effect variance
    ## --------------------------------------------------------------

    active_beta_sum_sq <- sum(
      beta[delta == 1L]^2
    )

    sigmaBetaSq <- (
      active_beta_sum_sq +
        nub * scaleb
    ) / rchisq(
      1L,
      df = nnz + nub
    )

    if (!is.finite(sigmaBetaSq) ||
        sigmaBetaSq <= 0) {
      sigmaBetaSq <- small_positive_value
    }

    ## --------------------------------------------------------------
    ## 7.12 Summary-data residual variance
    ##
    ## SSE/N =
    ## phenotype_variance
    ## - 2 beta' bhat
    ## + beta' LD beta
    ## --------------------------------------------------------------

    effect <- beta * delta
    LD_effect <- as.numeric(LD %*% effect)

    varg <- as.numeric(
      crossprod(effect, LD_effect)
    )

    if (!is.finite(varg) || varg < 0) {
      varg <- 0
    }

    residual_variance_component <-
      phenotype_variance -
      2 * sum(effect * bhat) +
      varg

    residual_variance_component <- max(
      residual_variance_component,
      small_positive_value
    )

    summary_sse <-
      N * residual_variance_component

    vare <- (
      summary_sse +
        nue * scalee
    ) / rchisq(
      1L,
      df = N + nue
    )

    if (!is.finite(vare) || vare <= 0) {
      vare <- small_positive_value
    }

    ## --------------------------------------------------------------
    ## 7.13 Heritability
    ## --------------------------------------------------------------

    h2 <- varg / (varg + vare)

    if (!is.finite(h2)) {
      h2 <- 0
    }

    ## --------------------------------------------------------------
    ## 7.14 Storage
    ## --------------------------------------------------------------

    store_index <- keep_lookup[iter]

    if (store_index > 0L) {
      if (store_beta) {
        beta_mcmc[store_index, ] <- beta
      }

      if (store_delta) {
        delta_mcmc[store_index, ] <- delta
      }

      Delta_mcmc[store_index, ] <- Delta
      gamma_mcmc[store_index, ] <- gamma

      keptIter[store_index, ] <- c(
        sigmaBetaSq,
        vare,
        varg,
        h2,
        alpha1,
        alpha2,
        alpha3,
        rho,
        nnz,
        mu_pi,
        mu_Pi
      )

      pip_snp_sum <- pip_snp_sum + delta
      pip_gene_sum <- pip_gene_sum + Delta
      pip_cell_sum <- pip_cell_sum + gamma
    }

    ## --------------------------------------------------------------
    ## 7.15 Progress
    ## --------------------------------------------------------------

    if (verbose && iter %% 100L == 0L) {
      cat(
        sprintf(
          paste0(
            "\n iter %4d, nnz = %4d, h2 = %6.3f, ",
            "alpha = (%.3f, %.3f, %.3f), ",
            "mu_pi = %.3f, mu_Pi = %.3f\n"
          ),
          iter,
          nnz,
          h2,
          alpha1,
          alpha2,
          alpha3,
          mu_pi,
          mu_Pi
        )
      )
    }
  }

  ## ================================================================
  ## 8. Return
  ## ================================================================

  result <- list(
    par = keptIter,
    beta = beta_mcmc,
    delta = delta_mcmc,
    Delta = Delta_mcmc,
    gamma = gamma_mcmc,
    pip_snp = pip_snp_sum / nkeep,
    pip_gene = pip_gene_sum / nkeep,
    pip_cell = pip_cell_sum / nkeep,
    burnin = burnin,
    thin = thin,
    nkeep = nkeep,
    niter = niter,
    N = N,
    phenotype_variance = phenotype_variance,
    dimensions = c(
      m = m,
      G = G,
      C = C
    )
  )

  class(result) <- c(
    "BayesMAPSummary",
    "BayesMAP"
  )

  result
}
