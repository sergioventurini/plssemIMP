make_missing <- function(data, prop = 0.5, method = "ampute", mech = "MCAR", missArgs = list()) {
  res <- list()

  if (method == "ampute") {
    patterns <- missArgs$patterns
    freq <- missArgs$freq
    weights <- missArgs$weights
    type <- missArgs$type
    cont <- missArgs$cont
    odds <- missArgs$odds
    data_mat <- as.matrix(data)
    data_amp <- mice::ampute(data = data_mat, prop = prop, patterns = patterns,
      freq = freq, mech = mech, weights = weights, type = type, cont = cont,
      odds = odds)$amp
  }
  else if (method == "raw") {
    if (is.null(missArgs$R))
      stop("when using 'raw' to generate missing values, the response indicator R must be provided.")
    R <- missArgs$R
    data_amp <- data
    for (i in 1:nrow(data)) {
      for (j in 1:ncol(data)) {
        if (R[i, j])
          data_amp[i, j] <- NA
      }
    }
  }
  else {
    stop("the specified missing generation method is not available.")
  }

  res$orig <- data
  res$amputed <- data_amp

  res
}
