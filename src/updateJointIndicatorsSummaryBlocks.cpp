#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

using namespace Rcpp;

namespace {

inline double clamp_probability_joint_summary_blocks(double p) {
  const double eps = 1e-12;
  
  if (!R_finite(p)) {
    return 0.5;
  }
  
  return std::max(
    eps,
    std::min(1.0 - eps, p)
  );
}


inline double log_bernoulli_joint_summary_blocks(
    int value,
    double probability
) {
  const double p =
    clamp_probability_joint_summary_blocks(probability);
  
  return value == 1 ?
  std::log(p) :
    std::log1p(-p);
}


inline double log_sum_exp_joint_summary_blocks(
    const NumericVector& values
) {
  const double maximum = max(values);
  
  double total = 0.0;
  
  for (R_xlen_t i = 0; i < values.size(); ++i) {
    total += std::exp(values[i] - maximum);
  }
  
  return maximum + std::log(total);
}

} // namespace


//' Jointly update indicators using one summary-data LD block
//'
//' Uses a block-specific LD matrix for the SNP likelihood while evaluating
//' the SNP, gene and cell hierarchy using the complete annotation matrices.
//' All indices supplied from R are one-based.
//'
//' @name updateJointIndicatorsSummaryBlocks_cpp
//' @keywords internal
// [[Rcpp::export]]
 List updateJointIndicatorsSummaryBlocks_cpp(
     NumericVector bhat_block,
     NumericMatrix LD_block,
     double N,
     NumericVector beta,
     IntegerVector delta,
     IntegerVector Delta,
     IntegerVector gamma,
     NumericMatrix L,
     NumericMatrix A,
     NumericMatrix B,
     IntegerVector block_indices,
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
   const int block_size = bhat_block.size();
   const int m = beta.size();
   const int G = L.ncol();
   const int C = A.ncol();
   const int number_of_joint_updates = snp_index.size();
   
   /*
    * Dimension checks
    */
   
   if (
       LD_block.nrow() != block_size ||
         LD_block.ncol() != block_size
   ) {
     stop(
       "LD_block must be block_size by block_size."
     );
   }
   
   if (block_indices.size() != block_size) {
     stop(
       "block_indices must have length equal to length(bhat_block)."
     );
   }
   
   if (
       delta.size() != m ||
         L.nrow() != m ||
         A.nrow() != m ||
         baseline.size() != m
   ) {
     stop(
       "beta, delta, L, A and baseline have incompatible dimensions."
     );
   }
   
   if (
       Delta.size() != G ||
         pathways.size() != G ||
         gamma.size() != C ||
         B.nrow() != G ||
         B.ncol() != C
   ) {
     stop(
       "Delta, gamma, B or pathways has an incompatible dimension."
     );
   }
   
   if (
       gene_index.size() != number_of_joint_updates ||
         cell_index.size() != number_of_joint_updates
   ) {
     stop(
       "snp_index, gene_index and cell_index must have equal lengths."
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
   
   /*
    * Map global SNP indices to local LD-block positions.
    */
   
   std::vector<int> global_to_local(
       static_cast<std::size_t>(m),
       -1
   );
   
   for (int local = 0; local < block_size; ++local) {
     const int global = block_indices[local] - 1;
     
     if (global < 0 || global >= m) {
       stop(
         "block_indices contains an invalid genome-wide SNP index."
       );
     }
     
     if (global_to_local[global] != -1) {
       stop(
         "block_indices contains duplicated SNP indices."
       );
     }
     
     global_to_local[global] = local;
   }
   
   rho = clamp_probability_joint_summary_blocks(rho);
   
   /*
    * Calculate the residual correlation vector for this LD block:
    *
    * rcorr = bhat_block - LD_block %*% (beta * delta)
    */
   
   NumericVector rcorr = clone(bhat_block);
   const double* LD_pointer = REAL(LD_block);
   
   for (int local_j = 0; local_j < block_size; ++local_j) {
     const int global_j = block_indices[local_j] - 1;
     
     const double effect =
       beta[global_j] *
       static_cast<double>(delta[global_j]);
     
     if (effect == 0.0) {
       continue;
     }
     
     const double* LD_column =
       LD_pointer +
       static_cast<std::size_t>(local_j) * block_size;
     
     for (int local_i = 0; local_i < block_size; ++local_i) {
       rcorr[local_i] -= LD_column[local_i] * effect;
     }
   }
   
   /*
    * Genome-wide hierarchical predictors
    */
   
   NumericVector LDelta(m);
   NumericVector Agamma(m);
   NumericVector Bgamma(G);
   
   for (int j = 0; j < m; ++j) {
     for (int g = 0; g < G; ++g) {
       LDelta[j] += L(j, g) * Delta[g];
     }
     
     for (int c = 0; c < C; ++c) {
       Agamma[j] += A(j, c) * gamma[c];
     }
   }
   
   for (int g = 0; g < G; ++g) {
     for (int c = 0; c < C; ++c) {
       Bgamma[g] += B(g, c) * gamma[c];
     }
   }
   
   IntegerMatrix sampled_states(
       number_of_joint_updates,
       3
   );
   
   NumericVector selected_probability(
       number_of_joint_updates
   );
   
   /*
    * Joint SNP-gene-cell updates
    */
   
   for (
       int joint_update = 0;
       joint_update < number_of_joint_updates;
       ++joint_update
   ) {
     const int global_j =
       snp_index[joint_update] - 1;
     
     const int g =
       gene_index[joint_update] - 1;
     
     const int c =
       cell_index[joint_update] - 1;
     
     if (
         global_j < 0 ||
           global_j >= m ||
           g < 0 ||
           g >= G ||
           c < 0 ||
           c >= C
     ) {
       stop(
         "Joint-update indices are outside their valid ranges."
       );
     }
     
     const int local_j =
       global_to_local[global_j];
     
     if (local_j < 0) {
       stop(
         "The selected SNP is not present in the supplied LD block."
       );
     }
     
     const int old_delta =
       delta[global_j];
     
     const int old_Delta =
       Delta[g];
     
     const int old_gamma =
       gamma[c];
     
     const double old_effect =
       beta[global_j] *
       static_cast<double>(old_delta);
     
     const double* LD_column =
       LD_pointer +
       static_cast<std::size_t>(local_j) * block_size;
     
     const double LD_jj =
       LD_column[local_j];
     
     if (!R_finite(LD_jj) || LD_jj <= 0.0) {
       stop(
         "All diagonal elements of LD_block must be positive and finite."
       );
     }
     
     /*
      * Remove the selected SNP effect.
      */
     
     if (old_effect != 0.0) {
       for (
           int local_i = 0;
           local_i < block_size;
           ++local_i
       ) {
         rcorr[local_i] +=
           LD_column[local_i] * old_effect;
       }
     }
     
     /*
      * Conditional SNP-effect likelihood
      */
     
     const double rhs =
       N * rcorr[local_j] / vare;
     
     const double precision =
       N * LD_jj / vare +
       1.0 / sigmaBetaSq;
     
     const double inverse_precision =
       1.0 / precision;
     
     const double beta_hat =
       inverse_precision * rhs;
     
     const double log_bayes_factor =
       0.5 * (
           std::log(inverse_precision) -
             std::log(sigmaBetaSq) +
             beta_hat * rhs
       );
     
     /*
      * Enumerate the eight possible joint states.
      */
     
     NumericVector log_weight(8);
     
     for (int state = 0; state < 8; ++state) {
       const int candidate_delta =
         state & 1;
       
       const int candidate_Delta =
         (state >> 1) & 1;
       
       const int candidate_gamma =
         (state >> 2) & 1;
       
       double log_posterior =
         candidate_delta == 1 ?
         log_bayes_factor :
         0.0;
       
       /*
        * Genome-wide SNP-indicator prior.
        */
       
       for (int k = 0; k < m; ++k) {
         if (
             k != global_j &&
               L(k, g) == 0.0 &&
               A(k, c) == 0.0
         ) {
           continue;
         }
         
         const double candidate_LDelta =
           LDelta[k] +
           L(k, g) *
           static_cast<double>(
             candidate_Delta - old_Delta
           );
         
         const double candidate_Agamma =
           Agamma[k] +
           A(k, c) *
           static_cast<double>(
             candidate_gamma - old_gamma
           );
         
         const double eta =
           mu_pi +
           candidate_LDelta * alpha1 +
           candidate_Agamma * alpha2 +
           baseline[k];
         
         const double probability =
           R::pnorm(
             eta,
             0.0,
             1.0,
             true,
             false
           );
         
         const int indicator =
           k == global_j ?
           candidate_delta :
           delta[k];
         
         log_posterior +=
           log_bernoulli_joint_summary_blocks(
             indicator,
             probability
           );
       }
       
       /*
        * Genome-wide gene-indicator prior.
        */
       
       for (int h = 0; h < G; ++h) {
         if (
             h != g &&
               B(h, c) == 0.0
         ) {
           continue;
         }
         
         const double candidate_Bgamma =
           Bgamma[h] +
           B(h, c) *
           static_cast<double>(
             candidate_gamma - old_gamma
           );
         
         const double xi =
           mu_Pi +
           candidate_Bgamma * alpha3 +
           pathways[h];
         
         const double probability =
           R::pnorm(
             xi,
             0.0,
             1.0,
             true,
             false
           );
         
         const int indicator =
           h == g ?
           candidate_Delta :
           Delta[h];
         
         log_posterior +=
           log_bernoulli_joint_summary_blocks(
             indicator,
             probability
           );
       }
       
       log_posterior +=
         log_bernoulli_joint_summary_blocks(
           candidate_gamma,
           rho
         );
       
       log_weight[state] = log_posterior;
     }
     
     /*
      * Convert log weights to probabilities.
      */
     
     const double normalizer =
       log_sum_exp_joint_summary_blocks(
         log_weight
       );
     
     NumericVector state_probability(8);
     
     for (int state = 0; state < 8; ++state) {
       state_probability[state] =
         std::exp(
           log_weight[state] - normalizer
         );
     }
     
     const int sampled_state =
       Rcpp::sample(
         8,
         1,
         false,
         state_probability
       )[0] - 1;
     
     const int new_delta =
       sampled_state & 1;
     
     const int new_Delta =
       (sampled_state >> 1) & 1;
     
     const int new_gamma =
       (sampled_state >> 2) & 1;
     
     /*
      * Update the selected state.
      */
     
     delta[global_j] = new_delta;
     Delta[g] = new_Delta;
     gamma[c] = new_gamma;
     
     if (new_delta == 1) {
       beta[global_j] =
         R::rnorm(
           beta_hat,
           std::sqrt(inverse_precision)
         );
       
       for (
           int local_i = 0;
           local_i < block_size;
           ++local_i
       ) {
         rcorr[local_i] -=
           LD_column[local_i] *
           beta[global_j];
       }
     } else {
       beta[global_j] = 0.0;
     }
     
     /*
      * Update cached hierarchical predictors.
      */
     
     const int Delta_change =
       new_Delta - old_Delta;
     
     if (Delta_change != 0) {
       for (int k = 0; k < m; ++k) {
         LDelta[k] +=
           Delta_change * L(k, g);
       }
     }
     
     const int gamma_change =
       new_gamma - old_gamma;
     
     if (gamma_change != 0) {
       for (int k = 0; k < m; ++k) {
         Agamma[k] +=
           gamma_change * A(k, c);
       }
       
       for (int h = 0; h < G; ++h) {
         Bgamma[h] +=
           gamma_change * B(h, c);
       }
     }
     
     sampled_states(joint_update, 0) =
       new_delta;
     
     sampled_states(joint_update, 1) =
       new_Delta;
     
     sampled_states(joint_update, 2) =
       new_gamma;
     
     selected_probability[joint_update] =
       state_probability[sampled_state];
   }
   
   colnames(sampled_states) =
     CharacterVector::create(
       "delta",
       "Delta",
       "gamma"
     );
   
   return List::create(
     Named("beta") = beta,
     Named("delta") = delta,
     Named("Delta") = Delta,
     Named("gamma") = gamma,
     Named("rcorr_block") = rcorr,
     Named("nnz") =
       std::accumulate(
         delta.begin(),
         delta.end(),
         0
       ),
       Named("sampled_states") =
         sampled_states,
         Named("selected_probability") =
           selected_probability
   );
 }