mode_function <- function(x, na.rm = TRUE) {
  if (na.rm) {
    ux <- unique(x[!is.na(x)])
  }
  else {
    ux <- unique(x)
  }
  ux[which.max(tabulate(match(x, ux)))]
}

euclidean_dist <- function(X, idx) {
  n <- nrow(X)
  D <- sqrt(rowSums((X - matrix(1, n, 1) %*% as.matrix(X[idx, ]))^2, na.rm = TRUE))
  D
}
