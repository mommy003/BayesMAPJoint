#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

inline double clamp_probability_gamma(double p) {
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

inline double probability_from_logpost_gamma(
    double logpost0,
    double logpost1
) {
  const double difference = logpost0 - logpost1;

  double probability;

  if (difference >= 0.0) {
    const double e = std::exp(-difference);
    probability = e / (1.0 + e);
  } else {
    const double e = std::exp(difference);
    probability = 1.0 / (1.0 + e);
  }

  return clamp_probability_gamma(probability);
}

// [[Rcpp::export]]
IntegerVector updateGamma_cpp(
    IntegerVector gamma,
    IntegerVector delta,
    IntegerVector Delta,
    NumericMatrix A,
    NumericMatrix B,
    NumericMatrix L,
    double alpha1,
    double alpha2,
    double alpha3,
    double mu_pi,
    double mu_Pi,
    NumericVector baseline,
    NumericVector pathways,
    double rho
) {
  const int m = A.nrow();
  const int C = A.ncol();
  const int G = B.nrow();

  // ---------------------------------------------------------------
  // Input checks
  // ---------------------------------------------------------------

  if (delta.size() != m) {
    stop("delta length must equal nrow(A).");
  }

  if (gamma.size() != C) {
    stop("gamma length must equal ncol(A).");
  }

  if (B.ncol() != C) {
    stop("B must have the same number of columns as A.");
  }

  if (Delta.size() != G) {
    stop("Delta length must equal nrow(B).");
  }

  if (L.nrow() != m || L.ncol() != G) {
    stop("L must be nrow(A) by nrow(B).");
  }

  if (baseline.size() != m) {
    stop("baseline length must equal nrow(A).");
  }

  if (pathways.size() != G) {
    stop("pathways length must equal nrow(B).");
  }

  rho = clamp_probability_gamma(rho);

  // ---------------------------------------------------------------
  // Cached matrix-vector products
  // ---------------------------------------------------------------

  NumericVector LDelta(m);
  NumericVector Agamma(m);
  NumericVector Bgamma(G);

  // LDelta = L %*% Delta
  for (int j = 0; j < m; ++j) {
    double value = 0.0;

    for (int g = 0; g < G; ++g) {
      value += L(j, g) * static_cast<double>(Delta[g]);
    }

    LDelta[j] = value;
  }

  // Agamma = A %*% gamma
  for (int j = 0; j < m; ++j) {
    double value = 0.0;

    for (int c = 0; c < C; ++c) {
      value += A(j, c) * static_cast<double>(gamma[c]);
    }

    Agamma[j] = value;
  }

  // Bgamma = B %*% gamma
  for (int g = 0; g < G; ++g) {
    double value = 0.0;

    for (int c = 0; c < C; ++c) {
      value += B(g, c) * static_cast<double>(gamma[c]);
    }

    Bgamma[g] = value;
  }

  // ---------------------------------------------------------------
  // Sequential Gibbs updates for gamma
  // ---------------------------------------------------------------

  for (int c = 0; c < C; ++c) {
    const int old_gamma_c = gamma[c];

    double logpost0 = 0.0;
    double logpost1 = 0.0;

    // SNP-level component
    for (int j = 0; j < m; ++j) {
      const double Agamma_without_c =
        Agamma[j] -
        A(j, c) * static_cast<double>(old_gamma_c);

      const double eta0 =
        mu_pi +
        LDelta[j] * alpha1 +
        Agamma_without_c * alpha2 +
        baseline[j];

      const double eta1 =
        eta0 +
        A(j, c) * alpha2;

      double p_snp0 = R::pnorm(
        eta0,
        0.0,
        1.0,
        true,
        false
      );

      double p_snp1 = R::pnorm(
        eta1,
        0.0,
        1.0,
        true,
        false
      );

      p_snp0 = clamp_probability_gamma(p_snp0);
      p_snp1 = clamp_probability_gamma(p_snp1);

      if (delta[j] == 1) {
        logpost0 += std::log(p_snp0);
        logpost1 += std::log(p_snp1);
      } else {
        logpost0 += std::log1p(-p_snp0);
        logpost1 += std::log1p(-p_snp1);
      }
    }

    // Gene-level component
    for (int g = 0; g < G; ++g) {
      const double Bgamma_without_c =
        Bgamma[g] -
        B(g, c) * static_cast<double>(old_gamma_c);

      const double xi0 =
        mu_Pi +
        Bgamma_without_c * alpha3 +
        pathways[g];

      const double xi1 =
        xi0 +
        B(g, c) * alpha3;

      double p_gene0 = R::pnorm(
        xi0,
        0.0,
        1.0,
        true,
        false
      );

      double p_gene1 = R::pnorm(
        xi1,
        0.0,
        1.0,
        true,
        false
      );

      p_gene0 = clamp_probability_gamma(p_gene0);
      p_gene1 = clamp_probability_gamma(p_gene1);

      if (Delta[g] == 1) {
        logpost0 += std::log(p_gene0);
        logpost1 += std::log(p_gene1);
      } else {
        logpost0 += std::log1p(-p_gene0);
        logpost1 += std::log1p(-p_gene1);
      }
    }

    logpost0 += std::log1p(-rho);
    logpost1 += std::log(rho);

    const double pr1 = probability_from_logpost_gamma(
      logpost0,
      logpost1
    );

    const int new_gamma_c = static_cast<int>(
      R::rbinom(1.0, pr1)
    );

    gamma[c] = new_gamma_c;

    // Incrementally update cached products
    const int change = new_gamma_c - old_gamma_c;

    if (change != 0) {
      for (int j = 0; j < m; ++j) {
        Agamma[j] +=
          static_cast<double>(change) * A(j, c);
      }

      for (int g = 0; g < G; ++g) {
        Bgamma[g] +=
          static_cast<double>(change) * B(g, c);
      }
    }
  }

  return gamma;
}