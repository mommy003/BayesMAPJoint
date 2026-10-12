#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace Rcpp;

namespace {

inline double clamp_probability_joint(double p) {
  const double eps = 1e-12;

  if (!R_finite(p)) {
    return 0.5;
  }

  return std::max(eps, std::min(1.0 - eps, p));
}

inline double log_bernoulli_joint(int value, double probability) {
  const double p = clamp_probability_joint(probability);
  return value == 1 ? std::log(p) : std::log1p(-p);
}

inline double log_sum_exp_joint(const NumericVector& values) {
  const double maximum = max(values);
  double total = 0.0;

  for (R_xlen_t i = 0; i < values.size(); ++i) {
    total += std::exp(values[i] - maximum);
  }

  return maximum + std::log(total);
}

} // namespace


//' Jointly update local SNP-gene-cell indicator blocks
//'
//' Each requested block contains one SNP, one gene and one cell. The
//' function enumerates all eight configurations of `(delta_j, Delta_g,
//' gamma_c)` and samples from their exact full conditional while holding
//' all other variables fixed. The SNP effect is analytically integrated
//' when the eight state probabilities are calculated and is subsequently
//' sampled conditionally if the selected SNP state is nonzero.
//'
//' Indices supplied from R are one-based.
//'
//' @keywords internal
// [[Rcpp::export]]
List updateJointIndicators_cpp(
    NumericMatrix X,
    NumericVector ycorr,
    NumericVector beta,
    IntegerVector delta,
    IntegerVector Delta,
    IntegerVector gamma,
    NumericMatrix L,
    NumericMatrix A,
    NumericMatrix B,
    IntegerVector snp_index,
    IntegerVector gene_index,
    IntegerVector cell_index,
    NumericVector xpx,
    double vare,
    double sigmaBetaSq,
    double alpha1,
    double alpha2,
    double alpha3,
    double mu_pi,
    double mu_Pi,
    NumericVector baseline,
    NumericVector pathways,
    double rho
) {
  const int n = X.nrow();
  const int m = X.ncol();
  const int G = L.ncol();
  const int C = A.ncol();
  const int number_of_blocks = snp_index.size();

  if (ycorr.size() != n || beta.size() != m || delta.size() != m) {
    stop("X, ycorr, beta and delta have incompatible dimensions.");
  }

  if (L.nrow() != m || A.nrow() != m || A.ncol() != C) {
    stop("L and A must have one row per SNP.");
  }

  if (Delta.size() != G || gamma.size() != C) {
    stop("Delta or gamma has an incompatible length.");
  }

  if (B.nrow() != G || B.ncol() != C) {
    stop("B must be ncol(L) by ncol(A).");
  }

  if (gene_index.size() != number_of_blocks ||
      cell_index.size() != number_of_blocks) {
    stop("snp_index, gene_index and cell_index must have equal lengths.");
  }

  if (xpx.size() != m || baseline.size() != m || pathways.size() != G) {
    stop("xpx, baseline or pathways has an incompatible length.");
  }

  if (!R_finite(vare) || vare <= 0.0 ||
      !R_finite(sigmaBetaSq) || sigmaBetaSq <= 0.0) {
    stop("vare and sigmaBetaSq must be positive and finite.");
  }

  rho = clamp_probability_joint(rho);

  // Cached predictors for the current state.
  NumericVector LDelta(m);
  NumericVector Agamma(m);
  NumericVector Bgamma(G);

  for (int j = 0; j < m; ++j) {
    double l_value = 0.0;
    double a_value = 0.0;

    for (int g = 0; g < G; ++g) {
      l_value += L(j, g) * static_cast<double>(Delta[g]);
    }

    for (int c = 0; c < C; ++c) {
      a_value += A(j, c) * static_cast<double>(gamma[c]);
    }

    LDelta[j] = l_value;
    Agamma[j] = a_value;
  }

  for (int g = 0; g < G; ++g) {
    double value = 0.0;

    for (int c = 0; c < C; ++c) {
      value += B(g, c) * static_cast<double>(gamma[c]);
    }

    Bgamma[g] = value;
  }

  IntegerMatrix sampled_states(number_of_blocks, 3);
  NumericVector selected_probability(number_of_blocks);

  for (int block = 0; block < number_of_blocks; ++block) {
    const int j = snp_index[block] - 1;
    const int g = gene_index[block] - 1;
    const int c = cell_index[block] - 1;

    if (j < 0 || j >= m || g < 0 || g >= G || c < 0 || c >= C) {
      stop("Joint block indices are outside their valid ranges.");
    }

    const int old_delta = delta[j];
    const int old_Delta = Delta[g];
    const int old_gamma = gamma[c];
    const double old_beta = beta[j];
    const double old_effect = old_beta * static_cast<double>(old_delta);
    const double* x_column = REAL(X) + static_cast<std::size_t>(j) * n;

    // Remove the current SNP effect. ycorr now conditions on all other SNPs.
    if (old_effect != 0.0) {
      for (int i = 0; i < n; ++i) {
        ycorr[i] += x_column[i] * old_effect;
      }
    }

    double cross_product = 0.0;
    for (int i = 0; i < n; ++i) {
      cross_product += x_column[i] * ycorr[i];
    }

    const double rhs = cross_product / vare;
    const double precision = xpx[j] / vare + 1.0 / sigmaBetaSq;
    const double inverse_precision = 1.0 / precision;
    const double beta_hat = inverse_precision * rhs;
    const double log_bayes_factor = 0.5 * (
      std::log(inverse_precision) -
      std::log(sigmaBetaSq) +
      beta_hat * rhs
    );

    NumericVector log_weight(8);

    for (int state = 0; state < 8; ++state) {
      const int candidate_delta = state & 1;
      const int candidate_Delta = (state >> 1) & 1;
      const int candidate_gamma = (state >> 2) & 1;
      double log_posterior = candidate_delta == 1 ? log_bayes_factor : 0.0;

      // The changed gene or cell can alter the SNP prior for every SNP
      // annotated to that gene or cell. Include all affected terms.
      for (int k = 0; k < m; ++k) {
        if (k != j && L(k, g) == 0.0 && A(k, c) == 0.0) {
          continue;
        }

        const double candidate_LDelta =
          LDelta[k] +
          L(k, g) * static_cast<double>(candidate_Delta - old_Delta);

        const double candidate_Agamma =
          Agamma[k] +
          A(k, c) * static_cast<double>(candidate_gamma - old_gamma);

        const double eta =
          mu_pi +
          candidate_LDelta * alpha1 +
          candidate_Agamma * alpha2 +
          baseline[k];

        const double probability = R::pnorm(
          eta, 0.0, 1.0, true, false
        );

        const int indicator = k == j ? candidate_delta : delta[k];
        log_posterior += log_bernoulli_joint(indicator, probability);
      }

      // Changing gamma_c can alter the gene prior for every gene connected
      // to that cell. Include the selected gene even if B(g, c) is zero.
      for (int h = 0; h < G; ++h) {
        if (h != g && B(h, c) == 0.0) {
          continue;
        }

        const double candidate_Bgamma =
          Bgamma[h] +
          B(h, c) * static_cast<double>(candidate_gamma - old_gamma);

        const double xi =
          mu_Pi +
          candidate_Bgamma * alpha3 +
          pathways[h];

        const double probability = R::pnorm(
          xi, 0.0, 1.0, true, false
        );

        const int indicator = h == g ? candidate_Delta : Delta[h];
        log_posterior += log_bernoulli_joint(indicator, probability);
      }

      log_posterior += log_bernoulli_joint(candidate_gamma, rho);
      log_weight[state] = log_posterior;
    }

    const double normalizer = log_sum_exp_joint(log_weight);
    NumericVector state_probability(8);

    for (int state = 0; state < 8; ++state) {
      state_probability[state] = std::exp(log_weight[state] - normalizer);
    }

    const int sampled_state = Rcpp::sample(
      8,
      1,
      false,
      state_probability
    )[0] - 1;

    const int new_delta = sampled_state & 1;
    const int new_Delta = (sampled_state >> 1) & 1;
    const int new_gamma = (sampled_state >> 2) & 1;

    delta[j] = new_delta;
    Delta[g] = new_Delta;
    gamma[c] = new_gamma;

    if (new_delta == 1) {
      const double new_beta = R::rnorm(
        beta_hat,
        std::sqrt(inverse_precision)
      );

      beta[j] = new_beta;

      for (int i = 0; i < n; ++i) {
        ycorr[i] -= x_column[i] * new_beta;
      }
    } else {
      beta[j] = 0.0;
    }

    const int Delta_change = new_Delta - old_Delta;
    if (Delta_change != 0) {
      for (int k = 0; k < m; ++k) {
        LDelta[k] += static_cast<double>(Delta_change) * L(k, g);
      }
    }

    const int gamma_change = new_gamma - old_gamma;
    if (gamma_change != 0) {
      for (int k = 0; k < m; ++k) {
        Agamma[k] += static_cast<double>(gamma_change) * A(k, c);
      }

      for (int h = 0; h < G; ++h) {
        Bgamma[h] += static_cast<double>(gamma_change) * B(h, c);
      }
    }

    sampled_states(block, 0) = new_delta;
    sampled_states(block, 1) = new_Delta;
    sampled_states(block, 2) = new_gamma;
    selected_probability[block] = state_probability[sampled_state];
  }

  const int nnz = std::accumulate(delta.begin(), delta.end(), 0);

  colnames(sampled_states) = CharacterVector::create(
    "delta", "Delta", "gamma"
  );

  return List::create(
    Named("beta") = beta,
    Named("delta") = delta,
    Named("Delta") = Delta,
    Named("gamma") = gamma,
    Named("ycorr") = ycorr,
    Named("nnz") = nnz,
    Named("sampled_states") = sampled_states,
    Named("selected_probability") = selected_probability
  );
}
