plot_results <- function(res, true_coefs, methods = "all", value = "est") {
  if (is.null(true_coefs))
    stop("true coefficient values must be provided.")

  nruns <- length(res) - 3
  nmeth <- length(res[[1]]) - 2
  methnm <- names(res[[1]])[1:(length(res[[1]]) - 2)]
  if (is.null(methods) | any(methods == "all")) {
    methods <- methnm
  }
  else {
    if (!all(methods %in% methnm))
      stop("the specified methods are not available.")
  }
  valuenm <- c("RB", "PB", "CR", "AW", "RMSE")
  if (is.null(value) | any(value == "all")) {
    value <- valuenm
  }
  else {
    if (!all(value %in% valuenm))
      stop("the specified values to graph are not available.")
  }

  out <- NA #ggplot2::ggplot()

  out
}
