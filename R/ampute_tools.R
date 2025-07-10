generate_patterns <- function(p, k, names = FALSE) {
  if (p < 2) {
    stop("p must be at least 2.")
  }
  if (k < 1) {
    stop("k must be at least 1.")
  }
  if (k > p) {
    stop("k must be at most equal to p.")
  }
  
  # get all combinations of k variables from 1 to p
  combs <- combn(p, k)
  
  # initialize the result matrix
  result <- matrix(1, nrow = ncol(combs), ncol = p)
  
  # fill the matrix: for each row, set 1 where the variable is present
  for (i in 1:ncol(combs)) {
    pair <- combs[, i]
    result[i, pair] <- 0
  }
  
  # optional: name the rows and columns
  if (names) {
    colnames(result) <- paste0("Var", 1:p)
    rownames(result) <- apply(combs, 2, function(x) paste0("(", x[1], ",", x[2], ")"))
  }
  
  return(data.frame(result))
}

generate_freqs <- function(patterns, reps = NULL, weights = NULL) {
  if (!is.null(weights)) {
    if (any(weights <= 0))
      stop("weights must be nonnegative.")

    weights <- weights/sum(weights)

    if (!is.null(reps)) {
      if (length(reps) != length(weights))
        stop("the reps and weights arguments must have the same number of elements.")
      if (sum(reps) != nrow(patterns))
        stop("the sum of reps must be equal to the number of patterns.")
    }
  }

  if (is.null(reps)) {
    result <- rep(1/nrow(patterns), nrow(patterns))
  }
  else {
    result <- numeric()
    for (k in 1:length(reps)) {
      result <- c(result, rep(weights[k]/reps[k], reps[k]))
    }
  }

  return(result)
}
