# Internal truncated-normal helper
rtruncnorm_scalar <- function(
    mean = 0,
    sd = 1,
    lower = -Inf,
    upper = Inf
) {
  if (!is.finite(sd) || sd <= 0) {
    if (is.finite(lower)) return(lower)
    if (is.finite(upper)) return(upper)
    return(mean)
  }

  lower_probability <- if (is.infinite(lower) && lower < 0) {
    0
  } else {
    pnorm((lower - mean) / sd)
  }

  upper_probability <- if (is.infinite(upper) && upper > 0) {
    1
  } else {
    pnorm((upper - mean) / sd)
  }

  epsilon <- 1e-15

  lower_probability <- min(
    max(lower_probability, epsilon),
    1 - epsilon
  )

  upper_probability <- min(
    max(upper_probability, epsilon),
    1 - epsilon
  )

  if (!is.finite(lower_probability) ||
      !is.finite(upper_probability) ||
      lower_probability >= upper_probability) {
    if (is.finite(lower)) return(lower)
    if (is.finite(upper)) return(upper)
    return(mean)
  }

  u <- runif(
    1L,
    min = lower_probability,
    max = upper_probability
  )

  out <- mean + sd * qnorm(u)

  if (!is.finite(out)) {
    if (is.finite(lower)) {
      out <- lower
    } else if (is.finite(upper)) {
      out <- upper
    } else {
      out <- mean
    }
  }

  out
}
