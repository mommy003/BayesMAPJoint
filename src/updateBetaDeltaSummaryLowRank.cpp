#include <Rcpp.h>
#include <cmath>
#include <cstddef>

using namespace Rcpp;

namespace {

inline double clamp_lowrank_probability(double p) {
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


inline double logistic_lowrank(double x) {
  if (x >= 0.0) {
    const double value = std::exp(-x);

    return clamp_lowrank_probability(
      1.0 / (1.0 + value)
    );
  }

  const double value = std::exp(x);

  return clamp_lowrank_probability(
    value / (1.0 + value)
  );
}

} // namespace


//' Update summary-data SNP effects using low-rank LD
//'
//' Uses the diagonal-corrected representation
//' `R = U diag(lambda) U' + diag(diagonal_correction)`.
//'
//' @keywords internal
// [[Rcpp::export]]
List updateBetaDeltaSummaryLowRank_cpp(
    NumericVector bhat,
    NumericMatrix U,
    NumericVector lambda,
    NumericVector diagonal_correction,
    double N,
    NumericVector beta,
    IntegerVector delta,
    NumericVector pi_j,
    double vare,
    double sigmaBetaSq
) {
  const int m = bhat.size();
  const int retained_rank = U.ncol();

  /*
   * Input validation
   */

  if (U.nrow() != m || retained_rank < 1) {
    stop(
      "U must have length(bhat) rows and at least one column."
    );
  }

  if (lambda.size() != retained_rank) {
    stop("length(lambda) must equal ncol(U).");
  }

  if (diagonal_correction.size() != m) {
    stop(
      "diagonal_correction must have length length(bhat)."
    );
  }

  if (
      beta.size() != m ||
      delta.size() != m ||
      pi_j.size() != m
  ) {
    stop(
      "beta, delta and pi_j must have length length(bhat)."
    );
  }

  if (
      !R_finite(N) ||
      N <= 0.0 ||
      !R_finite(vare) ||
      vare <= 0.0 ||
      !R_finite(sigmaBetaSq) ||
      sigmaBetaSq <= 0.0
  ) {
    stop(
      "N, vare and sigmaBetaSq must be positive and finite."
    );
  }

  for (int component = 0;
       component < retained_rank;
       ++component) {
    if (
        !R_finite(lambda[component]) ||
        lambda[component] <= 0.0
    ) {
      stop(
        "All retained eigenvalues must be positive and finite."
      );
    }
  }

  for (int j = 0; j < m; ++j) {
    if (
        !R_finite(bhat[j]) ||
        !R_finite(beta[j]) ||
        !R_finite(pi_j[j]) ||
        !R_finite(diagonal_correction[j]) ||
        diagonal_correction[j] < 0.0
    ) {
      stop(
        "Low-rank inputs contain invalid or non-finite values."
      );
    }

    if (delta[j] != 0 && delta[j] != 1) {
      stop("delta must contain only zero and one.");
    }
  }

  const double* U_pointer = REAL(U);

  /*
   * Current effects and projected effects:
   *
   * effect = beta * delta
   * projected_effect = U' * effect
   */

  NumericVector effect(m);

  for (int j = 0; j < m; ++j) {
    effect[j] =
      beta[j] *
      static_cast<double>(delta[j]);
  }

  NumericVector projected_effect(
    retained_rank
  );

  for (int component = 0;
       component < retained_rank;
       ++component) {
    const double* U_column =
      U_pointer +
      static_cast<std::size_t>(component) * m;

    double projected_value = 0.0;

    for (int j = 0; j < m; ++j) {
      projected_value +=
        U_column[j] * effect[j];
    }

    projected_effect[component] =
      projected_value;
  }

  /*
   * Diagonal of the approximated LD matrix.
   */

  NumericVector LD_diagonal(m);

  for (int j = 0; j < m; ++j) {
    double diagonal_value =
      diagonal_correction[j];

    for (int component = 0;
         component < retained_rank;
         ++component) {
      const double U_value =
        U_pointer[
          j +
          static_cast<std::size_t>(component) * m
        ];

      diagonal_value +=
        lambda[component] *
        U_value *
        U_value;
    }

    if (
        !R_finite(diagonal_value) ||
        diagonal_value <= 0.0
    ) {
      stop(
        "The reconstructed LD diagonal must be positive and finite."
      );
    }

    LD_diagonal[j] = diagonal_value;
  }

  const double inverse_vare =
    1.0 / vare;

  const double inverse_sigma_beta_sq =
    1.0 / sigmaBetaSq;

  const double log_sigma_beta_sq =
    std::log(sigmaBetaSq);

  int nnz = 0;

  /*
   * Sequential Gibbs updates.
   */

  for (int j = 0; j < m; ++j) {
    const double old_effect = effect[j];

    /*
     * Remove SNP j from U' * effect.
     */

    if (old_effect != 0.0) {
      for (int component = 0;
           component < retained_rank;
           ++component) {
        const double U_value =
          U_pointer[
            j +
            static_cast<std::size_t>(component) * m
          ];

        projected_effect[component] -=
          U_value * old_effect;
      }
    }

    /*
     * Conditional residual excluding SNP j:
     *
     * r_j = bhat_j -
     *       U_j Lambda U_{-j}' effect_{-j}
     *
     * The diagonal correction contributes only to SNP j,
     * which has already been removed.
     */

    double predicted_without_j = 0.0;

    for (int component = 0;
         component < retained_rank;
         ++component) {
      const double U_value =
        U_pointer[
          j +
          static_cast<std::size_t>(component) * m
        ];

      predicted_without_j +=
        U_value *
        lambda[component] *
        projected_effect[component];
    }

    const double residual_j =
      bhat[j] -
      predicted_without_j;

    const double rhs =
      N *
      residual_j *
      inverse_vare;

    const double precision =
      N *
      LD_diagonal[j] *
      inverse_vare +
      inverse_sigma_beta_sq;

    const double inverse_precision =
      1.0 / precision;

    const double beta_hat =
      inverse_precision * rhs;

    const double prior_probability =
      clamp_lowrank_probability(
        pi_j[j]
      );

    const double log_odds =
      0.5 * (
        std::log(inverse_precision) -
        log_sigma_beta_sq +
        beta_hat * rhs
      ) +
      std::log(prior_probability) -
      std::log1p(-prior_probability);

    const double inclusion_probability =
      logistic_lowrank(log_odds);

    const int new_delta =
      static_cast<int>(
        R::rbinom(
          1.0,
          inclusion_probability
        )
      );

    double new_effect = 0.0;

    if (new_delta == 1) {
      new_effect = R::rnorm(
        beta_hat,
        std::sqrt(inverse_precision)
      );

      ++nnz;
    }

    beta[j] = new_effect;
    delta[j] = new_delta;
    effect[j] = new_effect;

    /*
     * Add the new SNP effect to U' * effect.
     */

    if (new_effect != 0.0) {
      for (int component = 0;
           component < retained_rank;
           ++component) {
        const double U_value =
          U_pointer[
            j +
            static_cast<std::size_t>(component) * m
          ];

        projected_effect[component] +=
          U_value * new_effect;
      }
    }
  }

  /*
   * Final residual correlation:
   *
   * rcorr = bhat -
   *         U Lambda U' effect -
   *         D effect
   */

  NumericVector rcorr(m);

  for (int j = 0; j < m; ++j) {
    double fitted_value =
      diagonal_correction[j] *
      effect[j];

    for (int component = 0;
         component < retained_rank;
         ++component) {
      const double U_value =
        U_pointer[
          j +
          static_cast<std::size_t>(component) * m
        ];

      fitted_value +=
        U_value *
        lambda[component] *
        projected_effect[component];
    }

    rcorr[j] =
      bhat[j] -
      fitted_value;
  }

  return List::create(
    Named("beta") = beta,
    Named("delta") = delta,
    Named("rcorr") = rcorr,
    Named("nnz") = nnz,
    Named("projected_effect") =
      projected_effect
  );
}