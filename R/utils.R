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

matrix2vec <- function(x, names = TRUE) {
  if (!is.matrix(x)) {
    stop("argument x is not a matrix")
  }
  if (!is.numeric(x)) {
    stop("argument x is not a numeric matrix")
  }
  res <- t(t(as.vector(x)))

  if (names) {
    if (is.null(rownames(x))) {
      rownames(x) <- paste0("r", 1:nrow(x))
    }
    if (is.null(colnames(x))) {
      colnames(x) <- paste0("c", 1:nrow(x))
    }
    row_names <- rownames(x)
    col_names <- colnames(x)
    row_inds <- row(x)
    col_inds <- col(x)
    rownames(res) <- paste(row_names[row_inds], col_names[col_inds], sep = "_")
  }

  res
}
