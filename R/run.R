run_sims <- function(runs = 10) {
  res <- array(NA, dim = c(2, runs, 3))
  dimnames(res) <- list(c("norm.predict", "norm.nob"),
                          as.character(1:runs),
                          c("estimate", "2.5 %","97.5 %"))
  for (run in 1:runs) {
    cat(paste0("Simulation run ", run, " of ", runs, "\n"))
    data <- create_data(run = run)
    data <- make_missing(data)
    res[1, run, ] <- test_impute(data, method = "norm.predict", m = 2)
    res[2, run, ] <- test_impute(data, method = "norm.nob")
  }
  res
}
