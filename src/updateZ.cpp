#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// Draw one value from a truncated normal distribution using
// inverse-CDF sampling.
inline double rtruncnorm_one(
    const double mean,
    const double sd,
    const double lower,
    const double upper
) {
  if (!R_finite(mean)) {
    return 0.0;
  }
  
  if (!R_finite(sd) || sd <= 0.0) {
    if (R_finite(lower)) {
      return lower;
    }
    
    if (R_finite(upper)) {
      return upper;
    }
    
    return mean;
  }
  
  double lower_prob;
  
  if (lower == R_NegInf) {
    lower_prob = 0.0;
  } else {
    lower_prob = R::pnorm(
      (lower - mean) / sd,
      0.0,
      1.0,
      true,
      false
    );
  }
  
  double upper_prob;
  
  if (upper == R_PosInf) {
    upper_prob = 1.0;
  } else {
    upper_prob = R::pnorm(
      (upper - mean) / sd,
      0.0,
      1.0,
      true,
      false
    );
  }
  
  const double eps = 1e-15;
  
  lower_prob = std::max(
    eps,
    std::min(lower_prob, 1.0 - eps)
  );
  
  upper_prob = std::max(
    eps,
    std::min(upper_prob, 1.0 - eps)
  );
  
  if (!R_finite(lower_prob) ||
      !R_finite(upper_prob) ||
      lower_prob >= upper_prob) {
    if (R_finite(lower)) {
      return lower;
    }
    
    if (R_finite(upper)) {
      return upper;
    }
    
    return mean;
  }
  
  const double u = R::runif(
    lower_prob,
    upper_prob
  );
  
  double out = mean +
    sd * R::qnorm(
        u,
        0.0,
        1.0,
        true,
        false
    );
  
  if (!R_finite(out)) {
    if (R_finite(lower)) {
      out = lower;
    } else if (R_finite(upper)) {
      out = upper;
    } else {
      out = mean;
    }
  }
  
  if (R_finite(lower) && out < lower) {
    out = lower;
  }
  
  if (R_finite(upper) && out > upper) {
    out = upper;
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector updateZ_cpp(
    IntegerVector indicator,
    NumericVector eta
) {
  const int size = indicator.size();
  
  if (eta.size() != size) {
    stop("indicator and eta must have the same length.");
  }
  
  NumericVector z(size);
  
  for (int i = 0; i < size; ++i) {
    if (indicator[i] == 1) {
      z[i] = rtruncnorm_one(
        eta[i],
           1.0,
           0.0,
           R_PosInf
      );
    } else if (indicator[i] == 0) {
      z[i] = rtruncnorm_one(
        eta[i],
           1.0,
           R_NegInf,
           0.0
      );
    } else {
      stop("indicator must contain only 0 and 1.");
    }
  }
  
  return z;
}