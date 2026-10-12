#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// Keep probabilities away from exactly zero and one.
inline double clamp_probability_beta_delta(double p) {
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

// Numerically stable logistic transformation.
inline double logistic_beta_delta(double x) {
  double probability;

  if (x >= 0.0) {
    const double e = std::exp(-x);
    probability = 1.0 / (1.0 + e);
  } else {
    const double e = std::exp(x);
    probability = e / (1.0 + e);
  }

  return clamp_probability_beta_delta(probability);
}

// [[Rcpp::export]]
List updateBetaDelta_cpp(
    NumericMatrix X,
    NumericVector ycorr,
    NumericVector beta,
    IntegerVector delta,
    NumericVector pi_j,
    NumericVector xpx,
    double vare,
    double sigmaBetaSq
) {
  const int n = X.nrow();
  const int m = X.ncol();

  // ---------------------------------------------------------------
  // Input validation
  // ---------------------------------------------------------------

  if (ycorr.size() != n) {
    stop("ycorr length must equal nrow(X).");
  }

  if (beta.size() != m) {
    stop("beta length must equal ncol(X).");
  }

  if (delta.size() != m) {
    stop("delta length must equal ncol(X).");
  }

  if (pi_j.size() != m) {
    stop("pi_j length must equal ncol(X).");
  }

  if (xpx.size() != m) {
    stop("xpx length must equal ncol(X).");
  }

  if (!R_finite(vare) || vare <= 0.0) {
    stop("vare must be positive and finite.");
  }

  if (!R_finite(sigmaBetaSq) || sigmaBetaSq <= 0.0) {
    stop("sigmaBetaSq must be positive and finite.");
  }

  // Raw pointers reduce repeated Rcpp indexing overhead.
  double* X_ptr = REAL(X);
  double* ycorr_ptr = REAL(ycorr);
  double* beta_ptr = REAL(beta);
  int* delta_ptr = INTEGER(delta);
  double* pi_ptr = REAL(pi_j);
  double* xpx_ptr = REAL(xpx);

  const double inv_vare = 1.0 / vare;
  const double inv_sigma_beta_sq = 1.0 / sigmaBetaSq;
  const double log_sigma_beta_sq = std::log(sigmaBetaSq);

  int nnz = 0;

  // ---------------------------------------------------------------
  // Sequential SNP Gibbs updates
  // ---------------------------------------------------------------

  for (int j = 0; j < m; ++j) {
    // R matrices are column-major, so this points to the beginning
    // of column j. Its n values are contiguous in memory.
    const double* x_col = X_ptr + static_cast<std::size_t>(j) * n;

    const double old_beta = beta_ptr[j];
    const int old_delta = delta_ptr[j];
    const double old_effect =
      old_beta * static_cast<double>(old_delta);

    // Remove the old SNP effect from the residual.
    if (old_effect != 0.0) {
      for (int i = 0; i < n; ++i) {
        ycorr_ptr[i] += x_col[i] * old_effect;
      }
    }

    // Calculate X_j' ycorr without creating X[, j].
    double cross_product = 0.0;

    for (int i = 0; i < n; ++i) {
      cross_product += x_col[i] * ycorr_ptr[i];
    }

    const double rhs = cross_product * inv_vare;

    const double precision =
      xpx_ptr[j] * inv_vare +
      inv_sigma_beta_sq;

    const double invLhs = 1.0 / precision;
    const double betaHat = invLhs * rhs;

    const double prior_probability =
      clamp_probability_beta_delta(pi_ptr[j]);

    // log posterior odds:
    // log P(delta = 1 | ...) - log P(delta = 0 | ...)
    const double log_odds =
      0.5 * (
        std::log(invLhs) -
        log_sigma_beta_sq +
        betaHat * rhs
      ) +
      std::log(prior_probability) -
      std::log1p(-prior_probability);

    const double p_delta1 =
      logistic_beta_delta(log_odds);

    const int new_delta = static_cast<int>(
      R::rbinom(1.0, p_delta1)
    );

    delta_ptr[j] = new_delta;

    if (new_delta == 1) {
      const double new_beta = R::rnorm(
        betaHat,
        std::sqrt(invLhs)
      );

      beta_ptr[j] = new_beta;

      // Add the newly sampled effect to the model by removing it
      // from the residual.
      for (int i = 0; i < n; ++i) {
        ycorr_ptr[i] -= x_col[i] * new_beta;
      }

      ++nnz;
    } else {
      beta_ptr[j] = 0.0;
    }
  }

  return List::create(
    Named("beta") = beta,
    Named("delta") = delta,
    Named("ycorr") = ycorr,
    Named("nnz") = nnz
  );
}