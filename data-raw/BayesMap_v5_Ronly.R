
bayesmap_R <- function(
    X, y,
    L, A, B,
    baseline = NULL,
    pathways = NULL,
    niter = 2000,
    startPiSnp = 0.05,
    startPiGene = 0.10,
    startRho = 0.10,
    startH2 = 0.5,
    mu_pi = 0,      # initial value (CHANGED: now updated inside)
    mu_Pi = 0,      # initial value (CHANGED: now updated inside)
    sigmaAlphaSq = 1,
    nub = 4,
    nue = 4,
    verbose = TRUE
) {
  n <- nrow(X)
  m <- ncol(X)
  G <- ncol(L)
  C <- ncol(A)
  
  if (nrow(L) != m) stop("L must be m x G, where m = ncol(X).")
  if (nrow(A) != m) stop("A must be m x C, where m = ncol(X).")
  if (nrow(B) != G) stop("B must be G x C.")
  if (ncol(B) != C) stop("B must be G x C, where C = ncol(A).")
  
  if (is.null(baseline)) baseline <- rep(0, m)
  if (is.null(pathways)) pathways <- rep(0, G)
  
  if (length(baseline) != m) stop("baseline must have length m.")
  if (length(pathways) != G) stop("pathways must have length G.")
  
  ## changed for the correction of NA produced
  ## stops α from ever becoming NaN via the truncation step
  rtruncnorm_scalar <- function(mean = 0, sd = 1, lower = -Inf, upper = Inf) {
    if (!is.finite(sd) || sd <= 0) return(lower)
    
    a <- pnorm((lower - mean) / sd)
    b <- pnorm((upper - mean) / sd)
    
    if (!is.finite(a) || !is.finite(b) || a == b) return(lower)
    
    u <- runif(1, min = a, max = b)
    if (!is.finite(u)) return(lower)
    
    out <- mean + sd * qnorm(u)
    if (!is.finite(out)) out <- lower
    
    out
  }
  
  vary <- var(y)
  varg <- vary * startH2
  vare <- vary * (1 - startH2)
  
  sigmaBetaSq <- varg / m
  scaleb <- (nub - 2) / nub * sigmaBetaSq
  scalee <- (nue - 2) / nue * vare
  
  mu <- mean(y)
  
  beta  <- rep(0, m)
  delta <- rbinom(m, 1, startPiSnp)
  Delta <- rbinom(G, 1, startPiGene)
  gamma <- rbinom(C, 1, startRho)
  
  alpha1 <- 0
  alpha2 <- 0
  alpha3 <- 0
  rho    <- startRho
  
  z_pi <- rep(0, m)
  z_Pi <- rep(0, G)
  
  beta_mcmc  <- matrix(0, niter, m)
  delta_mcmc <- matrix(0, niter, m)
  Delta_mcmc <- matrix(0, niter, G)
  gamma_mcmc <- matrix(0, niter, C)
  
  # CHANGED 1: store mu_pi and mu_Pi as well
  keptIter <- matrix(0, niter, 12)
  colnames(keptIter) <- c(
    "Mu", "SigmaBetaSq", "Vare", "Varg", "h2",
    "Alpha1", "Alpha2", "Alpha3", "Rho", "Nnz",
    "mu_pi", "mu_Pi"
  )
  
  xb    <- X %*% (beta * delta)
  ycorr <- y - mu - xb
  xpx   <- colSums(X * X)
  
  for (iter in 1:niter) {
    ## 1. mu
    ycorr <- ycorr + mu
    muHat <- mean(ycorr)
    mu    <- rnorm(1, muHat, sqrt(vare / n))
    ycorr <- ycorr - mu
    
    ## 2. beta, delta
    LDelta <- as.vector(L %*% Delta)
    Agamma <- as.vector(A %*% gamma)
    eta_pi <- mu_pi + LDelta * alpha1 + Agamma * alpha2 + baseline
    pi_j   <- pnorm(eta_pi)
    pi_j   <- pmin(pmax(pi_j, 1e-12), 1 - 1e-12)
    
    nnz <- 0
    
    for (j in 1:m) {
      old_beta   <- beta[j]
      old_delta  <- delta[j]
      old_effect <- old_beta * old_delta
      
      if (old_effect != 0) {
        ycorr <- ycorr + X[, j] * old_effect
      }
      
      logp0 <- log(1 - pi_j[j])
      
      rhs    <- sum(X[, j] * ycorr) / vare
      invLhs <- 1 / (xpx[j] / vare + 1 / sigmaBetaSq)
      betaHat <- invLhs * rhs
      
      logp1 <- 0.5 * (log(invLhs) - log(sigmaBetaSq) + betaHat * rhs) + log(pi_j[j])
      
      p_delta1 <- 1 / (1 + exp(logp0 - logp1))
      ## to avoid NA (SNP inclusion)
      # ---- SAFETY PATCH ----
      if (!is.finite(p_delta1)) p_delta1 <- 0.5 # added 
      p_delta1 <- max(min(p_delta1, 1 - 1e-12), 1e-12) #added
      delta[j] <- rbinom(1, 1, p_delta1)
      
      
      if (delta[j] == 1) {
        beta[j] <- rnorm(1, betaHat, sqrt(invLhs))
        ycorr   <- ycorr - X[, j] * beta[j]
        nnz     <- nnz + 1
      } else {
        beta[j] <- 0
      }
    }
    
    xb <- X %*% (beta * delta)
    
    ## 3. z_pi
    LDelta <- as.vector(L %*% Delta)
    Agamma <- as.vector(A %*% gamma)
    eta_pi <- mu_pi + LDelta * alpha1 + Agamma * alpha2 + baseline
    
    for (j in 1:m) {
      if (delta[j] == 1)
        z_pi[j] <- rtruncnorm_scalar(mean = eta_pi[j], sd = 1, lower = 0)
      else
        z_pi[j] <- rtruncnorm_scalar(mean = eta_pi[j], sd = 1, upper = 0)
    }
    
    ## 4. Delta_g
    Bgamma <- as.vector(B %*% gamma)
    xi_g   <- mu_Pi + Bgamma * alpha3 + pathways
    Pi_g   <- pnorm(xi_g)
    Pi_g   <- pmin(pmax(Pi_g, 1e-12), 1 - 1e-12)
    
    for (g in 1:G) {
      lcol <- L[, g]
      
      eta0 <- mu_pi +
        (as.vector(L %*% Delta) - lcol * Delta[g]) * alpha1 +
        as.vector(A %*% gamma) * alpha2 + baseline
      p0 <- pnorm(eta0)
      p0 <- pmin(pmax(p0, 1e-12), 1 - 1e-12)
      
      logpost0 <- sum(delta * log(p0) + (1 - delta) * log(1 - p0)) +
        log(1 - Pi_g[g])
      
      eta1 <- eta0 + lcol * alpha1
      p1   <- pnorm(eta1)
      p1   <- pmin(pmax(p1, 1e-12), 1 - 1e-12)
      
      logpost1 <- sum(delta * log(p1) + (1 - delta) * log(1 - p1)) +
        log(Pi_g[g])
      
      pr1 <- 1 / (1 + exp(logpost0 - logpost1))
      ## To avoid NA (gene inclusion) 
      # ---- SAFETY PATCH ----
      if (!is.finite(pr1)) pr1 <- 0.5 # added
      pr1 <- max(min(pr1, 1 - 1e-12), 1e-12) # added
      Delta[g] <- rbinom(1, 1, pr1)
    }
    
    ## 5. z_Pi
    Bgamma <- as.vector(B %*% gamma)
    xi_g   <- mu_Pi + Bgamma * alpha3 + pathways
    
    for (g in 1:G) {
      if (Delta[g] == 1)
        z_Pi[g] <- rtruncnorm_scalar(mean = xi_g[g], sd = 1, lower = 0)
      else
        z_Pi[g] <- rtruncnorm_scalar(mean = xi_g[g], sd = 1, upper = 0)
    }
    
    ## 6. gamma_c
    for (c in 1:C) {
      acol <- A[, c]
      bcol <- B[, c]
      
      eta0 <- mu_pi + as.vector(L %*% Delta) * alpha1 +
        (as.vector(A %*% gamma) - acol * gamma[c]) * alpha2 + baseline
      p_snp0 <- pnorm(eta0)
      p_snp0 <- pmin(pmax(p_snp0, 1e-12), 1 - 1e-12)
      
      xi0 <- mu_Pi + (as.vector(B %*% gamma) - bcol * gamma[c]) * alpha3 + pathways
      p_gene0 <- pnorm(xi0)
      p_gene0 <- pmin(pmax(p_gene0, 1e-12), 1 - 1e-12)
      
      logpost0 <- sum(delta * log(p_snp0) + (1 - delta) * log(1 - p_snp0)) +
        sum(Delta * log(p_gene0) + (1 - Delta) * log(1 - p_gene0)) +
        log(1 - rho)
      
      eta1 <- eta0 + acol * alpha2
      p_snp1 <- pnorm(eta1)
      p_snp1 <- pmin(pmax(p_snp1, 1e-12), 1 - 1e-12)
      
      xi1 <- xi0 + bcol * alpha3
      p_gene1 <- pnorm(xi1)
      p_gene1 <- pmin(pmax(p_gene1, 1e-12), 1 - 1e-12)
      
      logpost1 <- sum(delta * log(p_snp1) + (1 - delta) * log(1 - p_snp1)) +
        sum(Delta * log(p_gene1) + (1 - Delta) * log(1 - p_gene1)) +
        log(rho)
      
      pr1 <- 1 / (1 + exp(logpost0 - logpost1))
      ## To avoid NA (for cell inclusion)
      # ---- SAFETY PATCH ----
      if (!is.finite(pr1)) pr1 <- 0.5 # added
      pr1 <- max(min(pr1, 1 - 1e-12), 1e-12) # added
      gamma[c] <- rbinom(1, 1, pr1)
    }
    
    ## 7. rho
    rho <- rbeta(1, sum(gamma) + 1, C - sum(gamma) + 1)
    
    ## 8. alpha1, alpha2
    LDelta <- as.vector(L %*% Delta)
    Agamma <- as.vector(A %*% gamma)
    
    # ---- alpha1 ---- 
    zcorr1 <- z_pi - mu_pi - Agamma * alpha2 - baseline
    C1 <- sum(LDelta^2) + 1 / sigmaAlphaSq
    r1 <- sum(LDelta * zcorr1)
    
    mean1 <- r1 / C1
    sd1   <- sqrt(1 / C1)
    
    if (!is.finite(mean1)) mean1 <- 0
    if (!is.finite(sd1) || sd1 <= 0) sd1 <- 1
    #alpha1 <- rnorm(1, r1 / C1, sqrt(1 / C1))
    alpha1 <- rtruncnorm_scalar(mean = mean1, sd = sd1, lower = 0) # Changed
    
    # ---- alpha2 ----
    zcorr2 <- z_pi - mu_pi - LDelta * alpha1 - baseline
    C2 <- sum(Agamma^2) + 1 / sigmaAlphaSq
    r2 <- sum(Agamma * zcorr2)
    mean2 <- r2 / C2
    sd2   <- sqrt(1 / C2)
    
    if (!is.finite(mean2)) mean2 <- 0
    if (!is.finite(sd2) || sd2 <= 0) sd2 <- 1
    #alpha2 <- rnorm(1, r2 / C2, sqrt(1 / C2))
    alpha2 <- rtruncnorm_scalar(mean = mean2, sd = sd2, lower = 0)   # CHANGED
    
    
    ## 9. alpha3
    Bgamma <- as.vector(B %*% gamma)
    zcorr3 <- z_Pi - mu_Pi - pathways
    C3 <- sum(Bgamma^2) + 1 / sigmaAlphaSq
    r3 <- sum(Bgamma * zcorr3)
    
    mean3 <- r3 / C3
    sd3   <- sqrt(1 / C3)
    
    if (!is.finite(mean3)) mean3 <- 0
    if (!is.finite(sd3) || sd3 <= 0) sd3 <- 1
    
    #alpha3 <- rnorm(1, r3 / C3, sqrt(1 / C3))
    alpha3 <- rtruncnorm_scalar(mean = mean3, sd = sd3, lower = 0)   # CHANGED
    
    
    ## 10. CHANGED 2: update mu_pi and mu_Pi (intercepts)
    # mu_pi | z_pi, LDelta, Agamma, baseline
    zcorr0 <- z_pi - LDelta * alpha1 - Agamma * alpha2 - baseline
    C0 <- m + 1 / sigmaAlphaSq
    r0 <- sum(zcorr0)
    mu_pi <- rnorm(1, r0 / C0, sqrt(1 / C0))
    
    # mu_Pi | z_Pi, Bgamma, pathways
    zcorrP <- z_Pi - Bgamma * alpha3 - pathways
    C0P <- G + 1 / sigmaAlphaSq
    r0P <- sum(zcorrP)
    mu_Pi <- rnorm(1, r0P / C0P, sqrt(1 / C0P))
    
    ## 11. sigmaBetaSq
    sigmaBetaSq <- (sum((beta[delta == 1])^2) + nub * scaleb) /
      rchisq(1, df = nnz + nub)
    
    ## 12. vare
    resid <- y - mu - X %*% (beta * delta)
    vare <- (sum(resid^2) + nue * scalee) / rchisq(1, df = n + nue)
    
    ## 13. store
    ghat <- as.vector(X %*% (beta * delta))
    varg <- var(ghat)
    h2   <- varg / (varg + vare)
    
    beta_mcmc[iter, ]  <- beta
    delta_mcmc[iter, ] <- delta
    Delta_mcmc[iter, ] <- Delta
    gamma_mcmc[iter, ] <- gamma
    
    keptIter[iter, ] <- c(
      mu, sigmaBetaSq, vare, varg, h2,
      alpha1, alpha2, alpha3, rho, nnz,
      mu_pi, mu_Pi
    )
    
    if (verbose && !(iter %% 100)) {
      cat(sprintf(
        "\n iter %4d, nnz = %4d, h2 = %6.3f, alpha = (%.3f, %.3f, %.3f), mu_pi = %.3f, mu_Pi = %.3f\n",
        iter, nnz, h2, alpha1, alpha2, alpha3, mu_pi, mu_Pi
      ))
    }
  }
  
  list(
    par      = keptIter,
    beta     = beta_mcmc,
    delta    = delta_mcmc,
    Delta    = Delta_mcmc,
    gamma    = gamma_mcmc,
    pip_snp  = colMeans(delta_mcmc),
    pip_gene = colMeans(Delta_mcmc),
    pip_cell = colMeans(gamma_mcmc)
  )
}
