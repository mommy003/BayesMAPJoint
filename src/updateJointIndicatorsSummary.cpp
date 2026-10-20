#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>

using namespace Rcpp;

namespace {

inline double clamp_probability_joint_summary(double p) {
  const double eps = 1e-12;
  if (!R_finite(p)) return 0.5;
  return std::max(eps, std::min(1.0 - eps, p));
}

inline double log_bernoulli_joint_summary(int value, double probability) {
  const double p = clamp_probability_joint_summary(probability);
  return value == 1 ? std::log(p) : std::log1p(-p);
}

inline double log_sum_exp_joint_summary(const NumericVector& values) {
  const double maximum = max(values);
  double total = 0.0;
  for (R_xlen_t i = 0; i < values.size(); ++i) {
    total += std::exp(values[i] - maximum);
  }
  return maximum + std::log(total);
}

} // namespace


//' Jointly update SNP-gene-cell indicators using summary statistics
//'
//' Enumerates all eight states of a connected SNP, gene and cell indicator
//' triplet. The SNP effect is integrated out using the same summary-data
//' likelihood as `updateBetaDeltaSummary_cpp()` and sampled conditionally
//' after a state is selected. Indices supplied from R are one-based.
//'
//' @keywords internal
// [[Rcpp::export]]
List updateJointIndicatorsSummary_cpp(
    NumericVector bhat,
    NumericMatrix LD,
    double N,
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
  const int m = bhat.size();
  const int G = L.ncol();
  const int C = A.ncol();
  const int number_of_blocks = snp_index.size();

  if (LD.nrow() != m || LD.ncol() != m) {
    stop("LD must be an m by m matrix, where m = length(bhat).");
  }
  if (beta.size() != m || delta.size() != m ||
      L.nrow() != m || A.nrow() != m) {
    stop("bhat, beta, delta, L and A have incompatible dimensions.");
  }
  if (Delta.size() != G || gamma.size() != C ||
      B.nrow() != G || B.ncol() != C) {
    stop("Delta, gamma or B has an incompatible dimension.");
  }
  if (gene_index.size() != number_of_blocks ||
      cell_index.size() != number_of_blocks) {
    stop("snp_index, gene_index and cell_index must have equal lengths.");
  }
  if (baseline.size() != m || pathways.size() != G) {
    stop("baseline or pathways has an incompatible length.");
  }
  if (!R_finite(N) || N <= 0.0 ||
      !R_finite(vare) || vare <= 0.0 ||
      !R_finite(sigmaBetaSq) || sigmaBetaSq <= 0.0) {
    stop("N, vare and sigmaBetaSq must be positive and finite.");
  }

  rho = clamp_probability_joint_summary(rho);

  NumericVector rcorr = clone(bhat);
  const double* LD_ptr = REAL(LD);
  for (int j = 0; j < m; ++j) {
    const double effect = beta[j] * static_cast<double>(delta[j]);
    if (effect == 0.0) continue;
    const double* LD_col = LD_ptr + static_cast<std::size_t>(j) * m;
    for (int i = 0; i < m; ++i) rcorr[i] -= LD_col[i] * effect;
  }

  NumericVector LDelta(m), Agamma(m), Bgamma(G);
  for (int j = 0; j < m; ++j) {
    for (int g = 0; g < G; ++g) LDelta[j] += L(j, g) * Delta[g];
    for (int c = 0; c < C; ++c) Agamma[j] += A(j, c) * gamma[c];
  }
  for (int g = 0; g < G; ++g) {
    for (int c = 0; c < C; ++c) Bgamma[g] += B(g, c) * gamma[c];
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
    const double old_effect = beta[j] * static_cast<double>(old_delta);
    const double* LD_col = LD_ptr + static_cast<std::size_t>(j) * m;
    const double LD_jj = LD_col[j];
    if (!R_finite(LD_jj) || LD_jj <= 0.0) {
      stop("All diagonal elements of LD must be positive and finite.");
    }

    // Remove the selected SNP effect; rcorr now conditions on all other SNPs.
    if (old_effect != 0.0) {
      for (int i = 0; i < m; ++i) rcorr[i] += LD_col[i] * old_effect;
    }

    const double rhs = N * rcorr[j] / vare;
    const double precision = N * LD_jj / vare + 1.0 / sigmaBetaSq;
    const double inverse_precision = 1.0 / precision;
    const double beta_hat = inverse_precision * rhs;
    const double log_bayes_factor = 0.5 * (
      std::log(inverse_precision) - std::log(sigmaBetaSq) + beta_hat * rhs
    );

    NumericVector log_weight(8);
    for (int state = 0; state < 8; ++state) {
      const int candidate_delta = state & 1;
      const int candidate_Delta = (state >> 1) & 1;
      const int candidate_gamma = (state >> 2) & 1;
      double log_posterior = candidate_delta == 1 ? log_bayes_factor : 0.0;

      for (int k = 0; k < m; ++k) {
        if (k != j && L(k, g) == 0.0 && A(k, c) == 0.0) continue;
        const double candidate_LDelta = LDelta[k] +
          L(k, g) * static_cast<double>(candidate_Delta - old_Delta);
        const double candidate_Agamma = Agamma[k] +
          A(k, c) * static_cast<double>(candidate_gamma - old_gamma);
        const double eta = mu_pi + candidate_LDelta * alpha1 +
          candidate_Agamma * alpha2 + baseline[k];
        const double probability = R::pnorm(eta, 0.0, 1.0, true, false);
        const int indicator = k == j ? candidate_delta : delta[k];
        log_posterior += log_bernoulli_joint_summary(indicator, probability);
      }

      for (int h = 0; h < G; ++h) {
        if (h != g && B(h, c) == 0.0) continue;
        const double candidate_Bgamma = Bgamma[h] +
          B(h, c) * static_cast<double>(candidate_gamma - old_gamma);
        const double xi = mu_Pi + candidate_Bgamma * alpha3 + pathways[h];
        const double probability = R::pnorm(xi, 0.0, 1.0, true, false);
        const int indicator = h == g ? candidate_Delta : Delta[h];
        log_posterior += log_bernoulli_joint_summary(indicator, probability);
      }

      log_posterior += log_bernoulli_joint_summary(candidate_gamma, rho);
      log_weight[state] = log_posterior;
    }

    const double normalizer = log_sum_exp_joint_summary(log_weight);
    NumericVector state_probability(8);
    for (int state = 0; state < 8; ++state) {
      state_probability[state] = std::exp(log_weight[state] - normalizer);
    }
    const int sampled_state = Rcpp::sample(8, 1, false, state_probability)[0] - 1;
    const int new_delta = sampled_state & 1;
    const int new_Delta = (sampled_state >> 1) & 1;
    const int new_gamma = (sampled_state >> 2) & 1;

    delta[j] = new_delta;
    Delta[g] = new_Delta;
    gamma[c] = new_gamma;
    if (new_delta == 1) {
      beta[j] = R::rnorm(beta_hat, std::sqrt(inverse_precision));
      for (int i = 0; i < m; ++i) rcorr[i] -= LD_col[i] * beta[j];
    } else {
      beta[j] = 0.0;
    }

    const int Delta_change = new_Delta - old_Delta;
    if (Delta_change != 0) {
      for (int k = 0; k < m; ++k) LDelta[k] += Delta_change * L(k, g);
    }
    const int gamma_change = new_gamma - old_gamma;
    if (gamma_change != 0) {
      for (int k = 0; k < m; ++k) Agamma[k] += gamma_change * A(k, c);
      for (int h = 0; h < G; ++h) Bgamma[h] += gamma_change * B(h, c);
    }

    sampled_states(block, 0) = new_delta;
    sampled_states(block, 1) = new_Delta;
    sampled_states(block, 2) = new_gamma;
    selected_probability[block] = state_probability[sampled_state];
  }

  colnames(sampled_states) = CharacterVector::create("delta", "Delta", "gamma");
  return List::create(
    Named("beta") = beta,
    Named("delta") = delta,
    Named("Delta") = Delta,
    Named("gamma") = gamma,
    Named("rcorr") = rcorr,
    Named("nnz") = std::accumulate(delta.begin(), delta.end(), 0),
    Named("sampled_states") = sampled_states,
    Named("selected_probability") = selected_probability
  );
}
