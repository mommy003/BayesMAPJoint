#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

inline double clamp_probability_delta(double p) {
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

inline double logistic_from_log_odds_delta(double log_odds) {
  double probability;

  if (log_odds >= 0.0) {
    const double e = std::exp(-log_odds);
    probability = 1.0 / (1.0 + e);
  } else {
    const double e = std::exp(log_odds);
    probability = e / (1.0 + e);
  }

  return clamp_probability_delta(probability);
}

// [[Rcpp::export]]
IntegerVector updateDelta_cpp(
    IntegerVector Delta,
    IntegerVector delta,
    NumericMatrix L,
    IntegerVector gamma,
    NumericMatrix A,
    NumericMatrix B,
    double alpha1,
    double alpha2,
    double alpha3,
    double mu_pi,
    double mu_Pi,
    NumericVector baseline,
    NumericVector pathways,
    double rho
) {
  const int m = L.nrow();
  const int G = L.ncol();
  const int C = A.ncol();

  // rho is not used directly in the gene update.
  (void)rho;

  if (delta.size() != m) {
    stop("delta length must equal nrow(L).");
  }

  if (Delta.size() != G) {
    stop("Delta length must equal ncol(L).");
  }

  if (A.nrow() != m) {
    stop("A must have the same number of rows as L.");
  }

  if (gamma.size() != C) {
    stop("gamma length must equal ncol(A).");
  }

  if (B.nrow() != G || B.ncol() != C) {
    stop("B must be ncol(L) by ncol(A).");
  }

  if (baseline.size() != m) {
    stop("baseline length must equal nrow(L).");
  }

  if (pathways.size() != G) {
    stop("pathways length must equal ncol(L).");
  }

  // ---------------------------------------------------------------
  // Cache matrix-vector products
  // ---------------------------------------------------------------

  NumericVector LDelta(m);
  NumericVector Agamma(m);
  NumericVector Bgamma(G);

  // LDelta = L %*% Delta
  for (int g = 0; g < G; ++g) {
    if (Delta[g] == 0) {
      continue;
    }

    for (int j = 0; j < m; ++j) {
      const double l_value = L(j, g);

      if (l_value != 0.0) {
        LDelta[j] += l_value;
      }
    }
  }

  // Agamma = A %*% gamma
  for (int c = 0; c < C; ++c) {
    if (gamma[c] == 0) {
      continue;
    }

    for (int j = 0; j < m; ++j) {
      const double a_value = A(j, c);

      if (a_value != 0.0) {
        Agamma[j] += a_value;
      }
    }
  }

  // Bgamma = B %*% gamma
  for (int c = 0; c < C; ++c) {
    if (gamma[c] == 0) {
      continue;
    }

    for (int g = 0; g < G; ++g) {
      const double b_value = B(g, c);

      if (b_value != 0.0) {
        Bgamma[g] += b_value;
      }
    }
  }

  // ---------------------------------------------------------------
  // Sequential Gibbs updates
  // ---------------------------------------------------------------

  for (int g = 0; g < G; ++g) {
    const int old_Delta_g = Delta[g];

    const double xi_g =
      mu_Pi +
      Bgamma[g] * alpha3 +
      pathways[g];

    double Pi_g = R::pnorm(
      xi_g,
      0.0,
      1.0,
      true,
      false
    );

    Pi_g = clamp_probability_delta(Pi_g);

    // Start with prior log odds:
    // log(Pi_g / (1 - Pi_g))
    double log_odds =
      std::log(Pi_g) -
      std::log1p(-Pi_g);

    // Only SNPs linked to gene g can distinguish Delta[g] = 1
    // from Delta[g] = 0.
    for (int j = 0; j < m; ++j) {
      const double l_value = L(j, g);

      if (l_value == 0.0) {
        continue;
      }

      const double L_without_g =
        LDelta[j] -
        l_value * static_cast<double>(old_Delta_g);

      const double eta0 =
        mu_pi +
        L_without_g * alpha1 +
        Agamma[j] * alpha2 +
        baseline[j];

      const double eta1 =
        eta0 +
        l_value * alpha1;

      double p0 = R::pnorm(
        eta0,
        0.0,
        1.0,
        true,
        false
      );

      double p1 = R::pnorm(
        eta1,
        0.0,
        1.0,
        true,
        false
      );

      p0 = clamp_probability_delta(p0);
      p1 = clamp_probability_delta(p1);

      if (delta[j] == 1) {
        log_odds += std::log(p1) - std::log(p0);
      } else {
        log_odds += std::log1p(-p1) - std::log1p(-p0);
      }
    }

    const double pr1 =
      logistic_from_log_odds_delta(log_odds);

    const int new_Delta_g = static_cast<int>(
      R::rbinom(1.0, pr1)
    );

    Delta[g] = new_Delta_g;

    // Preserve sequential Gibbs updating.
    const int change = new_Delta_g - old_Delta_g;

    if (change != 0) {
      for (int j = 0; j < m; ++j) {
        const double l_value = L(j, g);

        if (l_value != 0.0) {
          LDelta[j] +=
            static_cast<double>(change) * l_value;
        }
      }
    }
  }

  return Delta;
}