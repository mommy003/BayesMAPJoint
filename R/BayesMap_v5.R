#' Fit the BayesMAP model
#'
#' Fits a Bayesian multi-level association mapping model.
#'
#' @param X Numeric genotype or design matrix.
#' @param y Numeric phenotype vector.
#' @param L SNP-to-gene annotation matrix.
#' @param A SNP-to-cell annotation matrix.
#' @param B Gene-to-cell annotation matrix.
#' @param baseline Optional SNP-level baseline annotation vector.
#' @param pathways Optional gene-level pathway annotation vector.
#' @param niter Number of MCMC iterations.
#' @param niter Total number of MCMC iterations.
#' @param burnin Number of initial MCMC iterations discarded as burn-in.
#' @param thin Positive integer giving the thinning interval for stored samples.
#' @param store_beta Logical; whether to store the full posterior beta chain.
#' @param store_delta Logical; whether to store the full posterior delta chain.
#' @param joint_blocks_per_iter Number of supplementary local
#'   SNP-gene-cell blocks updated jointly per iteration. The default zero
#'   preserves the original sampler. Use a small positive value through
#'   [bayesmap_joint()] for the experimental blocked sampler.
#' @param startPiSnp Initial SNP inclusion probability.
#' @param startPiSnp Initial SNP inclusion probability.
#' @param startPiGene Initial gene inclusion probability.
#' @param startRho Initial cell inclusion probability.
#' @param startH2 Initial heritability.
#' @param mu_pi Initial SNP-level probit intercept.
#' @param mu_Pi Initial gene-level probit intercept.
#' @param sigmaAlphaSq Prior variance for alpha parameters.
#' @param nub Degrees-of-freedom parameter for SNP-effect variance.
#' @param nue Degrees-of-freedom parameter for residual variance.
#' @param verbose Logical; print progress during sampling.
#'
#' @return A list containing posterior samples and posterior inclusion
#'   probabilities.
#' @export
bayesmap <- function(
    X, y,
    L, A, B,
    baseline = NULL,
    pathways = NULL,
    niter = 2000,
    burnin = 0,
    thin = 1,
    store_beta = TRUE,
    store_delta = TRUE,
    joint_blocks_per_iter = 0L,
    startPiSnp = 0.05,
    startPiGene = 0.10,
    startRho = 0.10,
    startH2 = 0.5,
    mu_pi = 0,
    mu_Pi = 0,
    sigmaAlphaSq = 1,
    nub = 4,
    nue = 4,
    verbose = TRUE
) {
  ## ------------------------------------------------------------------
  ## 1. Dimensions and input validation
  ## ------------------------------------------------------------------
  
  n <- nrow(X)
  m <- ncol(X)
  G <- ncol(L)
  C <- ncol(A)
  
  if (length(y) != n) {
    stop("length(y) must equal nrow(X).")
  }
  
  if (nrow(L) != m) {
    stop("L must be m x G, where m = ncol(X).")
  }
  
  if (nrow(A) != m) {
    stop("A must be m x C, where m = ncol(X).")
  }
  
  if (nrow(B) != G) {
    stop("B must be G x C.")
  }
  
  if (ncol(B) != C) {
    stop("B must be G x C, where C = ncol(A).")
  }
  
  if (is.null(baseline)) {
    baseline <- numeric(m)
  }
  
  if (is.null(pathways)) {
    pathways <- numeric(G)
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
    stop("burnin must be an integer between 0 and niter - 1.")
  }
  
  burnin <- as.integer(burnin)
  
  if (length(thin) != 1L ||
      !is.finite(thin) ||
      thin < 1 ||
      thin != as.integer(thin)) {
    stop("thin must be a positive integer.")
  }
  
  thin <- as.integer(thin)
  
  if (!is.logical(store_beta) || length(store_beta) != 1L ||
      is.na(store_beta)) {
    stop("store_beta must be TRUE or FALSE.")
  }
  
  if (!is.logical(store_delta) || length(store_delta) != 1L ||
      is.na(store_delta)) {
    stop("store_delta must be TRUE or FALSE.")
  }

  if (length(joint_blocks_per_iter) != 1L ||
      !is.finite(joint_blocks_per_iter) ||
      joint_blocks_per_iter < 0 ||
      joint_blocks_per_iter != as.integer(joint_blocks_per_iter)) {
    stop("joint_blocks_per_iter must be a non-negative integer.")
  }

  joint_blocks_per_iter <- as.integer(joint_blocks_per_iter)
  
  if (!is.logical(verbose) || length(verbose) != 1L ||
      is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  
  if (!is.finite(startH2) || startH2 <= 0 || startH2 >= 1) {
    stop("startH2 must be strictly between 0 and 1.")
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
  
  if (!is.finite(sigmaAlphaSq) || sigmaAlphaSq <= 0) {
    stop("sigmaAlphaSq must be positive.")
  }
  
  if (!is.finite(nub) || nub <= 2) {
    stop("nub must be greater than 2.")
  }
  
  if (!is.finite(nue) || nue <= 2) {
    stop("nue must be greater than 2.")
  }
  
  ## ------------------------------------------------------------------
  ## 2. Burn-in, thinning and storage indices
  ## ------------------------------------------------------------------
  
  keep_iterations <- seq.int(
    from = burnin + 1L,
    to = niter,
    by = thin
  )
  
  nkeep <- length(keep_iterations)
  
  keep_lookup <- integer(niter)
  keep_lookup[keep_iterations] <- seq_len(nkeep)
  
  ## ------------------------------------------------------------------
  ## 3. Truncated-normal helper
  ## ------------------------------------------------------------------
  
  rtruncnorm_scalar <- function(
    mean = 0,
    sd = 1,
    lower = -Inf,
    upper = Inf
  ) {
    if (!is.finite(sd) || sd <= 0) {
      if (is.finite(lower)) {
        return(lower)
      }
      
      if (is.finite(upper)) {
        return(upper)
      }
      
      return(mean)
    }
    
    a <- pnorm((lower - mean) / sd)
    b <- pnorm((upper - mean) / sd)
    
    if (!is.finite(a) || !is.finite(b) || a >= b) {
      if (is.finite(lower)) {
        return(lower)
      }
      
      if (is.finite(upper)) {
        return(upper)
      }
      
      return(mean)
    }
    
    u <- runif(1L, min = a, max = b)
    
    if (!is.finite(u)) {
      if (is.finite(lower)) {
        return(lower)
      }
      
      if (is.finite(upper)) {
        return(upper)
      }
      
      return(mean)
    }
    
    out <- mean + sd * qnorm(u)
    
    if (!is.finite(out)) {
      if (is.finite(lower)) {
        return(lower)
      }
      
      if (is.finite(upper)) {
        return(upper)
      }
      
      return(mean)
    }
    
    out
  }
  
  ## ------------------------------------------------------------------
  ## 4. Initial values
  ## ------------------------------------------------------------------
  
  vary <- var(y)
  
  if (!is.finite(vary) || vary <= 0) {
    stop("y must have a positive, finite variance.")
  }
  
  varg <- vary * startH2
  vare <- vary * (1 - startH2)
  
  sigmaBetaSq <- varg / m
  
  scaleb <- ((nub - 2) / nub) * sigmaBetaSq
  scalee <- ((nue - 2) / nue) * vare
  
  mu <- mean(y)
  
  beta <- numeric(m)
  delta <- rbinom(m, size = 1L, prob = startPiSnp)
  Delta <- rbinom(G, size = 1L, prob = startPiGene)
  gamma <- rbinom(C, size = 1L, prob = startRho)
  
  alpha1 <- 0
  alpha2 <- 0
  alpha3 <- 0
  
  rho <- startRho
  
  z_pi <- numeric(m)
  z_Pi <- numeric(G)
  
  ## ------------------------------------------------------------------
  ## 5. Posterior-storage objects
  ## ------------------------------------------------------------------
  
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
    ncol = 12L
  )
  
  colnames(keptIter) <- c(
    "Mu",
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
  
  ## ------------------------------------------------------------------
  ## 6. Quantities used repeatedly
  ## ------------------------------------------------------------------
  
  xb <- as.vector(X %*% (beta * delta))
  ycorr <- y - mu - xb
  xpx <- colSums(X * X)
  
  ## ------------------------------------------------------------------
  ## 7. MCMC sampler
  ## ------------------------------------------------------------------
  
  for (iter in seq_len(niter)) {
    ## --------------------------------------------------------------
    ## 7.1 Intercept mu
    ## --------------------------------------------------------------
    
    ycorr <- ycorr + mu
    
    muHat <- mean(ycorr)
    
    mu <- rnorm(
      1L,
      mean = muHat,
      sd = sqrt(vare / n)
    )
    
    ycorr <- ycorr - mu
    
    ## --------------------------------------------------------------
    ## 7.2 SNP effects beta and SNP indicators delta
    ## --------------------------------------------------------------
    
    LDelta <- as.vector(L %*% Delta)
    Agamma <- as.vector(A %*% gamma)
    
    eta_pi <- mu_pi +
      LDelta * alpha1 +
      Agamma * alpha2 +
      baseline
    
    pi_j <- pnorm(eta_pi)
    
    pi_j <- pmin(
      pmax(pi_j, 1e-12),
      1 - 1e-12
    )
    
    tmp <- updateBetaDelta_cpp(
      X = X,
      ycorr = ycorr,
      beta = beta,
      delta = delta,
      pi_j = pi_j,
      xpx = xpx,
      vare = vare,
      sigmaBetaSq = sigmaBetaSq
    )
    
    beta <- tmp$beta
    delta <- tmp$delta
    ycorr <- tmp$ycorr
    nnz <- tmp$nnz
    
    ## --------------------------------------------------------------
    ## 7.3 Latent SNP probit variables z_pi
    ## --------------------------------------------------------------
    
    LDelta <- as.vector(L %*% Delta)
    Agamma <- as.vector(A %*% gamma)
    
    eta_pi <- mu_pi +
      LDelta * alpha1 +
      Agamma * alpha2 +
      baseline
    
    z_pi <- updateZ_cpp(
      indicator = delta,
      eta = eta_pi
    )
    
    ## --------------------------------------------------------------
    ## 7.4 Gene indicators Delta
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
    ## 7.5 Latent gene probit variables z_Pi
    ## --------------------------------------------------------------
    
    Bgamma <- as.vector(B %*% gamma)
    
    xi_g <- mu_Pi +
      Bgamma * alpha3 +
      pathways
    
    z_Pi <- updateZ_cpp(
      indicator = Delta,
      eta = xi_g
    )
    
    ## --------------------------------------------------------------
    ## 7.6 Cell or annotation indicators gamma
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
    ## 7.6a Experimental local joint indicator updates
    ##
    ## Each block jointly samples one connected
    ## (delta[j], Delta[g], gamma[c]) triplet by enumerating its eight
    ## possible states. The original separate updates above are retained,
    ## and the default joint_blocks_per_iter = 0 reproduces the original
    ## sampler exactly.
    ## --------------------------------------------------------------

    if (joint_blocks_per_iter > 0L) {
      joint_blocks <- .sample_joint_indicator_blocks(
        L = L,
        A = A,
        B = B,
        number = joint_blocks_per_iter
      )

      if (nrow(joint_blocks) > 0L) {
        joint_out <- updateJointIndicators_cpp(
          X = X,
          ycorr = ycorr,
          beta = beta,
          delta = delta,
          Delta = Delta,
          gamma = gamma,
          L = L,
          A = A,
          B = B,
          snp_index = joint_blocks[, "snp"],
          gene_index = joint_blocks[, "gene"],
          cell_index = joint_blocks[, "cell"],
          xpx = xpx,
          vare = vare,
          sigmaBetaSq = sigmaBetaSq,
          alpha1 = alpha1,
          alpha2 = alpha2,
          alpha3 = alpha3,
          mu_pi = mu_pi,
          mu_Pi = mu_Pi,
          baseline = baseline,
          pathways = pathways,
          rho = rho
        )

        beta <- joint_out$beta
        delta <- joint_out$delta
        Delta <- joint_out$Delta
        gamma <- joint_out$gamma
        ycorr <- joint_out$ycorr
        nnz <- joint_out$nnz

        ## Refresh both latent probit vectors because the blocked move can
        ## change indicators at every level.
        LDelta <- as.vector(L %*% Delta)
        Agamma <- as.vector(A %*% gamma)
        eta_pi <- mu_pi +
          LDelta * alpha1 +
          Agamma * alpha2 +
          baseline

        z_pi <- updateZ_cpp(
          indicator = delta,
          eta = eta_pi
        )

        Bgamma <- as.vector(B %*% gamma)
        xi_g <- mu_Pi +
          Bgamma * alpha3 +
          pathways

        z_Pi <- updateZ_cpp(
          indicator = Delta,
          eta = xi_g
        )
      }
    }
    
    ## --------------------------------------------------------------
    ## 7.7 Cell inclusion probability rho
    ## --------------------------------------------------------------
    
    gamma_sum <- sum(gamma)
    
    rho <- rbeta(
      1L,
      shape1 = gamma_sum + 1,
      shape2 = C - gamma_sum + 1
    )
    
    ## --------------------------------------------------------------
    ## 7.8 Alpha1 and Alpha2
    ## --------------------------------------------------------------
    
    LDelta <- as.vector(L %*% Delta)
    Agamma <- as.vector(A %*% gamma)
    
    ## Alpha1
    
    zcorr1 <- z_pi -
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
      lower = 0,
      upper = Inf
    )
    
    ## Alpha2
    
    zcorr2 <- z_pi -
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
      lower = 0,
      upper = Inf
    )
    
    ## --------------------------------------------------------------
    ## 7.9 Alpha3
    ## --------------------------------------------------------------
    
    Bgamma <- as.vector(B %*% gamma)
    
    zcorr3 <- z_Pi -
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
      lower = 0,
      upper = Inf
    )
    
    ## --------------------------------------------------------------
    ## 7.10 Probit intercepts mu_pi and mu_Pi
    ## --------------------------------------------------------------
    
    zcorr0 <- z_pi -
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
    
    zcorrP <- z_Pi -
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
    
    active_beta_sum_sq <- sum(beta[delta == 1L]^2)
    
    sigmaBetaSq <- (
      active_beta_sum_sq +
        nub * scaleb
    ) / rchisq(
      1L,
      df = nnz + nub
    )
    
    ## Numerical safeguard
    
    if (!is.finite(sigmaBetaSq) || sigmaBetaSq <= 0) {
      sigmaBetaSq <- .Machine$double.eps
    }
    
    ## --------------------------------------------------------------
    ## 7.12 Residual variance
    ##
    ## ycorr is already:
    ## y - mu - X %*% (beta * delta)
    ## --------------------------------------------------------------
    
    resid <- ycorr
    
    vare <- (
      sum(resid^2) +
        nue * scalee
    ) / rchisq(
      1L,
      df = n + nue
    )
    
    if (!is.finite(vare) || vare <= 0) {
      vare <- .Machine$double.eps
    }
    
    ## --------------------------------------------------------------
    ## 7.13 Genetic variance and heritability
    ##
    ## Avoid another X %*% (...) multiplication.
    ## --------------------------------------------------------------
    
    ghat <- as.vector(y - mu - ycorr)
    
    varg <- var(ghat)
    
    if (!is.finite(varg) || varg < 0) {
      varg <- 0
    }
    
    h2 <- varg / (varg + vare)
    
    ## --------------------------------------------------------------
    ## 7.14 Store selected iterations only
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
        mu,
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
    ## 7.15 Progress report
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
  
  ## ------------------------------------------------------------------
  ## 8. Return results
  ## ------------------------------------------------------------------
  
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
    dimensions = c(
      n = n,
      m = m,
      G = G,
      C = C
    )
  )
  
  class(result) <- "BayesMAP"
  
  result
}


# Draw a small set of connected SNP-gene-cell blocks without constructing
# a potentially enormous list of every nonzero annotation triplet.
.sample_joint_indicator_blocks <- function(L, A, B, number) {
  if (number <= 0L) {
    return(matrix(integer(), nrow = 0L, ncol = 3L,
                  dimnames = list(NULL, c("snp", "gene", "cell"))))
  }

  m <- nrow(L)
  selected <- matrix(
    NA_integer_,
    nrow = number,
    ncol = 3L,
    dimnames = list(NULL, c("snp", "gene", "cell"))
  )

  filled <- 0L
  attempts <- 0L
  maximum_attempts <- max(100L, 25L * number)

  while (filled < number && attempts < maximum_attempts) {
    attempts <- attempts + 1L
    j <- sample.int(m, 1L)

    gene_candidates <- which(L[j, ] != 0)
    if (length(gene_candidates) == 0L) {
      next
    }

    gene_weights <- abs(L[j, gene_candidates])
    g <- gene_candidates[which.max(gene_weights)]

    cell_score <- abs(A[j, ]) + abs(B[g, ])
    if (!any(cell_score > 0)) {
      next
    }

    c_index <- which.max(cell_score)
    filled <- filled + 1L
    selected[filled, ] <- c(j, g, c_index)
  }

  if (filled == 0L) {
    return(selected[FALSE, , drop = FALSE])
  }

  selected[seq_len(filled), , drop = FALSE]
}


#' Fit BayesMAP with supplementary joint indicator blocks
#'
#' Experimental version of [bayesmap()] that adds exact local blocked Gibbs
#' updates for connected SNP-gene-cell indicator triplets. All other model
#' components and returned objects are unchanged.
#'
#' @param ... Arguments passed to [bayesmap()].
#' @param joint_blocks_per_iter Number of local triplets jointly updated per
#'   MCMC iteration. Small values such as 1--5 are recommended initially.
#'
#' @return A `BayesMAP` fit with an additional `joint_blocks_per_iter` field.
#' @export
bayesmap_joint <- function(..., joint_blocks_per_iter = 5L) {
  fit <- bayesmap(
    ...,
    joint_blocks_per_iter = joint_blocks_per_iter
  )

  fit$joint_blocks_per_iter <- as.integer(joint_blocks_per_iter)
  class(fit) <- c("BayesMAPJoint", class(fit))
  fit
}
