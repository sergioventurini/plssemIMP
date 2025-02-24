make_missing <- function(data, prop = 0.5, method = "ampute", mech = "MCAR", missArgs = list(),
  seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  res <- list()
  res$orig <- data

  if (method == "ampute") {
    patterns <- missArgs$patterns
    freq <- missArgs$freq
    weights <- missArgs$weights
    type <- missArgs$type
    data <- mice::ampute(data = data, prop = prop, patterns = patterns, freq = freq, mech = mech,
      weights = weights, type = type)$amp
  }
  else if (method == "raw") {
    if (mech == "MCAR") {
      n <- nrow(data)
      p <- ncol(data)
      for (j in 1:p) {
        rx <- rbinom(n, 1, prop)
        data[rx == 0, j] <- NA
      }
    }
  }
  else {
    stop("the specified method is not implemented.")
  }

  res$amputed <- data

  res
}
