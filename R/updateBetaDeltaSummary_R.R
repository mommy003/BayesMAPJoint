updateBetaDeltaSummary_R <- function(
    bhat,
    LD,
    N,
    beta,
    delta,
    pi_j,
    vare,
    sigmaBetaSq
) {
  m <- length(bhat)

  if (!is.matrix(LD) || any(dim(LD) != c(m, m))) {
    stop("LD must be an m by m numeric matrix.")
  }

  if (length(beta) != m) {
    stop("beta must have length m.")
  }

  if (length(delta) != m) {
    stop("delta must have length m.")
  }

  if (length(pi_j) != m) {
    stop("pi_j must have length m.")
  }

  if (length(N) != 1L || !is.finite(N) || N <= 0) {
    stop("N must be a positive scalar.")
  }

  if (!is.finite(vare) || vare <= 0) {
    stop("vare must be positive.")
  }

  if (!is.finite(sigmaBetaSq) || sigmaBetaSq <= 0) {
    stop("sigmaBetaSq must be positive.")
  }

  joint_effect <- beta * delta

  ## LD-adjusted marginal residual:
  ## bhat - LD %*% joint_effect
  rcorr <- bhat - as.numeric(LD %*% joint_effect)

  nnz <- 0L

  for (j in seq_len(m)) {
    old_effect <- beta[j] * delta[j]

    ## Add the current SNP effect back before updating SNP j
    r_j <- rcorr[j] + LD[j, j] * old_effect

    rhs <- N * r_j / vare

    invLhs <- 1 / (
      N * LD[j, j] / vare +
        1 / sigmaBetaSq
    )

    betaHat <- invLhs * rhs

    prior_probability <- min(
      max(pi_j[j], 1e-12),
      1 - 1e-12
    )

    logp0 <- log1p(-prior_probability)

    logp1 <- 0.5 * (
      log(invLhs) -
        log(sigmaBetaSq) +
        betaHat * rhs
    ) + log(prior_probability)

    log_odds <- logp1 - logp0

    if (log_odds >= 0) {
      p_delta1 <- 1 / (1 + exp(-log_odds))
    } else {
      exp_log_odds <- exp(log_odds)
      p_delta1 <- exp_log_odds / (1 + exp_log_odds)
    }

    if (!is.finite(p_delta1)) {
      p_delta1 <- 0.5
    }

    p_delta1 <- min(
      max(p_delta1, 1e-12),
      1 - 1e-12
    )

    new_delta <- rbinom(
      1L,
      size = 1L,
      prob = p_delta1
    )

    if (new_delta == 1L) {
      new_beta <- rnorm(
        1L,
        mean = betaHat,
        sd = sqrt(invLhs)
      )

      nnz <- nnz + 1L
    } else {
      new_beta <- 0
    }

    new_effect <- new_beta * new_delta
    change <- new_effect - old_effect

    beta[j] <- new_beta
    delta[j] <- new_delta

    ## Incrementally update all residual marginal associations
    if (change != 0) {
      rcorr <- rcorr - LD[, j] * change
    }
  }

  list(
    beta = beta,
    delta = delta,
    rcorr = rcorr,
    nnz = nnz
  )
}
