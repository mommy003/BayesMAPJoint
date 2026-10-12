profile_file <- tempfile(fileext = ".out")

Rprof(
  profile_file,
  interval = 0.001,
  memory.profiling = TRUE
)

# Run BayesMAP here

Rprof(NULL)

profile <- summaryRprof(
  profile_file,
  memory = "both"
)

unlink(profile_file)