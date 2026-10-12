test_that("BayesMAP runs", {
  
  set.seed(123)
  
  fit <- bayesmap(
    X = X,
    y = y,
    L = L,
    A = A,
    B = B,
    niter = 10,
    verbose = FALSE
  )
  
  expect_true(is.list(fit))
  
  expect_true(
    all(
      c(
        "par",
        "beta",
        "delta",
        "Delta",
        "gamma",
        "pip_snp",
        "pip_gene",
        "pip_cell"
      ) %in% names(fit)
    )
  )
  
  expect_equal(
    ncol(fit$gamma),
    ncol(A)
  )
  
})