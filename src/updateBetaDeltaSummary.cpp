#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>

using namespace Rcpp;

// Keep probabilities away from exactly 0 and 1.
inline double clamp_summary_probability(double p) {
  const double eps = 1e-12;

  if (!R_finite(p)) {
    return 0.5;
  }

  if (p < eps) {
    return eps;
  }

  if (p > 1.0 - eps) {
    return 1.0 - eps;
  }

  return p;
}

// Numerically stable logistic function.
inline double logistic_summary(double x) {
  double probability;

  if (x >= 0.0) {
    const double e = std::exp(-x);
    probability = 1.0 / (1.0 + e);
  } else {
    const double e = std::exp(x);
    probability = e / (1.0 + e);
  }

  return clamp_summary_probability(probability);
}

// [[Rcpp::export]]
List updateBetaDeltaSummary_cpp(
    NumericVector bhat,
    NumericMatrix LD,
    double N,
    NumericVector beta,
    IntegerVector delta,
    NumericVector pi_j,
    double vare,
    double sigmaBetaSq
) {
  const int m = bhat.size();

  // ---------------------------------------------------------------
  // Input validation
  // ---------------------------------------------------------------

  if (LD.nrow() != m || LD.ncol() != m) {
    stop("LD must be an m by m matrix, where m = length(bhat).");
  }

  if (beta.size() != m) {
    stop("beta must have length m.");
  }

  if (delta.size() != m) {
    stop("delta must have length m.");
  }

  if (pi_j.size() != m) {
    stop("pi_j must have length m.");
  }

  if (!R_finite(N) || N <= 0.0) {
    stop("N must be positive and finite.");
  }

  if (!R_finite(vare) || vare <= 0.0) {
    stop("vare must be positive and finite.");
  }

  if (!R_finite(sigmaBetaSq) || sigmaBetaSq <= 0.0) {
    stop("sigmaBetaSq must be positive and finite.");
  }

  // ---------------------------------------------------------------
  // Raw pointers
  // ---------------------------------------------------------------

  const double* bhat_ptr = REAL(bhat);
  const double* LD_ptr = REAL(LD);

  double* beta_ptr = REAL(beta);
  int* delta_ptr = INTEGER(delta);
  const double* pi_ptr = REAL(pi_j);

  const double inv_vare = 1.0 / vare;
  const double inv_sigma_beta_sq = 1.0 / sigmaBetaSq;
  const double log_sigma_beta_sq = std::log(sigmaBetaSq);

  // ---------------------------------------------------------------
  // rcorr = bhat - LD %*% (beta * delta)
  // ---------------------------------------------------------------

  NumericVector rcorr(m);
  double* rcorr_ptr = REAL(rcorr);

  for (int i = 0; i < m; ++i) {
    rcorr_ptr[i] = bhat_ptr[i];
  }

  for (int j = 0; j < m; ++j) {
    const double effect =
      beta_ptr[j] * static_cast<double>(delta_ptr[j]);

    if (effect == 0.0) {
      continue;
    }

    // Column j is contiguous because R matrices are column-major.
    const double* LD_col =
      LD_ptr + static_cast<std::size_t>(j) * m;

    for (int i = 0; i < m; ++i) {
      rcorr_ptr[i] -= LD_col[i] * effect;
    }
  }

  int nnz = 0;

  // ---------------------------------------------------------------
  // Sequential SNP Gibbs updates
  // ---------------------------------------------------------------

  for (int j = 0; j < m; ++j) {
    const double* LD_col =
      LD_ptr + static_cast<std::size_t>(j) * m;

    const double LD_jj = LD_col[j];

    if (!R_finite(LD_jj) || LD_jj <= 0.0) {
      stop("All diagonal elements of LD must be positive and finite.");
    }

    const double old_beta = beta_ptr[j];
    const int old_delta = delta_ptr[j];

    const double old_effect =
      old_beta * static_cast<double>(old_delta);

    // Add the current SNP effect back for its conditional update.
    const double r_j =
      rcorr_ptr[j] +
      LD_jj * old_effect;

    const double rhs =
      N * r_j * inv_vare;

    const double precision =
      N * LD_jj * inv_vare +
      inv_sigma_beta_sq;

    const double invLhs = 1.0 / precision;
    const double betaHat = invLhs * rhs;

    const double prior_probability =
      clamp_summary_probability(pi_ptr[j]);

    // log posterior odds:
    // log P(delta_j = 1 | ...) - log P(delta_j = 0 | ...)
    const double log_odds =
      0.5 * (
          std::log(invLhs) -
            log_sigma_beta_sq +
            betaHat * rhs
      ) +
        std::log(prior_probability) -
        std::log1p(-prior_probability);

    const double p_delta1 =
      logistic_summary(log_odds);

    const int new_delta = static_cast<int>(
      R::rbinom(1.0, p_delta1)
    );

    double new_beta = 0.0;

    if (new_delta == 1) {
      new_beta = R::rnorm(
        betaHat,
        std::sqrt(invLhs)
      );

      ++nnz;
    }

    const double new_effect =
      new_beta * static_cast<double>(new_delta);

    const double change =
      new_effect - old_effect;

    beta_ptr[j] = new_beta;
    delta_ptr[j] = new_delta;

    // rcorr = bhat - LD %*% current_effect
    if (change != 0.0) {
      for (int i = 0; i < m; ++i) {
        rcorr_ptr[i] -= LD_col[i] * change;
      }
    }
  }

  return List::create(
    Named("beta") = beta,
    Named("delta") = delta,
    Named("rcorr") = rcorr,
    Named("nnz") = nnz
  );
}
