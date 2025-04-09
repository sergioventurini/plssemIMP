round_exact <- function(x, digits = 0) {
  round(x + sign(x) * 1e-10, digits)
}

rdirichlet <- function (n, alpha) {
  n <- as.numeric(n)
  if (!is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow = 1)
  }
  if (prod(dim(alpha)) == 0) {
    stop2("alpha should be non-empty.")
  }
  if (isTRUE(any(alpha <= 0))) {
    stop2("alpha must be positive.")
  }
  if (n == 1) {
    n <- nrow(alpha)
  }
  if (n > nrow(alpha)) {
    alpha <- matrix(alpha, nrow = n, ncol = ncol(alpha), byrow = TRUE)
  }
  x <- matrix(rgamma(ncol(alpha) * n, alpha), ncol = ncol(alpha))
  x/rowSums(x)
}

list_2_array <- function(obj) {
  if (is.data.frame(obj)) {
    obj
  }
  obj_len <- length(names(obj))
  obj_nm_len <- nchar(names(obj)[1])
  dat <- data.frame()
  for (i in 1:obj_len) {
    dat <- rbind(dat, as.integer(unlist(strsplit(names(obj)[i], "_"))))
  }
  dat <- cbind(dat, unlist(obj, use.names = FALSE))
  rownames(dat) <- names(obj)
  colnames(dat) <- c(paste0("index_", 1:ceiling(obj_nm_len/2)), "prmax")
  x_seq_len <- max(dat[, 1]) + 1
  
  dat
}

row_match <- function(x, table, nomatch = NA) {
  if (inherits(table, "matrix"))
    table <- as.data.frame(table)
  if (is.null(dim(x)))
    x <- as.data.frame(matrix(x, nrow = 1))
  cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
  ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
  match(cx, ct, nomatch = nomatch)
}

aggregate_sum_df <- function(x, by) {
  x <- as.data.frame(x)
  y <- as.data.frame(by, stringsAsFactors = FALSE)
  ident <- function(x) {
    y <- as.factor(x)
    l <- length(levels(y))
    s <- as.character(seq_len(l))
    n <- nchar(s)
    levels(y) <- paste0(strrep("0", n[l] - n), s)
    y
  }
  grp <- lapply(y, ident)
  grp <- if (ncol(y)) {
    names(grp) <- NULL
    do.call(paste, c(rev(grp), list(sep = ".")))
  }
  else integer(NROW(x))
  y <- y[match(sort(unique(grp)), grp, 0L), , drop = FALSE]
  z <- lapply(x, function(e) {
    ans <- lapply(X = unname(split(e, grp)), FUN = "sum")
    if (length(len <- unique(lengths(ans))) == 1L) {
      if (len == 1L) {
        cl <- lapply(ans, oldClass)
        cl1 <- cl[[1L]]
        ans <- if (!is.null(cl1) & all(vapply(cl, identical, NA, y = cl1)))
          do.call(c, ans)
        else unlist(ans, recursive = FALSE, use.names = FALSE)
      }
      else if (len > 1L)
        ans <- matrix(unlist(ans, recursive = FALSE, 
          use.names = FALSE), ncol = len, byrow = TRUE, 
          dimnames = if (!is.null(nms <- names(ans[[1L]]))) 
            list(NULL, nms))
    }
    ans
  })
  len <- length(y)
  for (i in seq_along(z)) y[[len + i]] <- z[[i]]
  y
}

is_square_matrix <- function(x) {
  if (!is.matrix(x)) 
    stop("argument x is not a matrix")
  return(nrow(x) == ncol(x))
}

is_symmetric_matrix <- function(x) {
  if (!is.matrix(x)) {
    stop("argument x is not a matrix")
  }
  if (!is.numeric(x)) {
    stop("argument x is not a numeric matrix")
  }
  if (!is_square_matrix(x)) 
    stop("argument x is not a square numeric matrix")
  return(sum(x == t(x)) == (nrow(x)^2))
}

is_positive_definite <- function(x, tol = 1e-08) {
  if (!is_square_matrix(x)) 
    stop("argument x is not a square matrix")
  if (!is_symmetric_matrix(x)) 
    stop("argument x is not a symmetric matrix")
  if (!is.numeric(x)) 
    stop("argument x is not a numeric matrix")
  eigenvalues <- eigen(x, only.values = TRUE)$values
  n <- nrow(x)
  for (i in 1:n) {
    if (abs(eigenvalues[i]) < tol) {
      eigenvalues[i] <- 0
    }
  }
  if (any(eigenvalues <= 0)) {
    return(FALSE)
  }
  return(TRUE)
}

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
