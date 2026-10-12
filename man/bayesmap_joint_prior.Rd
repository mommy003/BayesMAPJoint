% Experimental documentation; regenerate with roxygen2 before release.
\name{bayesmap_joint_prior}
\alias{bayesmap_joint_prior}
\title{Fit sparse-prior BayesMAP with joint indicator blocks}
\usage{
bayesmap_joint_prior(
  X,
  y,
  L,
  A,
  B,
  baseline = NULL,
  pathways = NULL,
  niter = 2000,
  burnin = 0,
  thin = 1,
  store_beta = TRUE,
  store_delta = TRUE,
  joint_blocks_per_iter = 5L,
  startPiSnp = 0.05,
  startPiGene = 0.1,
  startRho = 0.1,
  startH2 = 0.5,
  mu_pi = 0,
  mu_Pi = NULL,
  priorPiGene = 0.01,
  priorRho = 0.05,
  rhoPriorStrength = 20,
  sigmaMuPiSq = 1,
  sigmaAlphaSq = 1,
  nub = 4,
  nue = 4,
  verbose = TRUE
)
}
\arguments{
\item{X, y, L, A, B}{Individual-level data and SNP-to-gene,
SNP-to-cell and gene-to-cell annotation matrices.}
\item{baseline, pathways}{Optional SNP- and gene-level annotations.}
\item{niter, burnin, thin}{MCMC storage controls.}
\item{store_beta, store_delta}{Whether to retain the complete SNP chains.}
\item{joint_blocks_per_iter}{Number of connected SNP-gene-cell triplets
sampled jointly at each iteration.}
\item{startPiSnp, startPiGene, startRho, startH2}{Initial values.}
\item{mu_pi, mu_Pi}{Initial probit intercepts.}
\item{priorPiGene, priorRho, rhoPriorStrength, sigmaMuPiSq}{Sparse
higher-level prior controls.}
\item{sigmaAlphaSq, nub, nue}{Other prior controls.}
\item{verbose}{Whether to print progress.}
}
\value{
A \code{BayesMAP} posterior fit containing SNP, gene and cell PIPs.
}
\description{
Experimental sparse-prior sampler that supplements the original separate
indicator updates with exact local eight-state Gibbs updates for connected
SNP-gene-cell triplets.
}
