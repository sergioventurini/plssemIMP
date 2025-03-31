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
        ans <- if (!is.null(cl1) && all(vapply(cl, identical, NA, y = cl1)))
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

combine_vcov <- function(vcov_path, vcov_load) {
  P <- nrow(vcov_path)
  Q <- nrow(vcov_load)
  PQ <- P + Q
  vcov_all <- matrix(0, nrow = PQ, ncol = PQ)
  vcov_all[1:P, 1:P] <- vcov_path
  vcov_all[(P + 1):PQ, (P + 1):PQ] <- vcov_load
  if (!is_positive_definite(vcov_all))
    stop("variance-covariance matrix not positive definite.")

  vcov_all
}

rubin_est <- function(est_list) {
  est_mat <- do.call(rbind, est_list)
  
  colMeans(est_mat)
}

rubin_vcov <- function(est_list, vcov_list) {
  est_mat <- do.call(rbind, est_list)
  m <- nrow(est_mat)
  est_bar <- colMeans(est_mat)
  deviations <- est_mat - matrix(est_bar, nrow = m, ncol = length(est_bar), byrow = TRUE)
  B_hat <- t(deviations) %*% deviations / (m - 1)
  W_hat <- Reduce("+", vcov_list)/m
  T_hat <- W_hat + (1 + 1/m)*B_hat

  list(T = T_hat, B = B_hat, W = W_hat)
}

barnard_rubin <- function(fitMI, what = "path") {
  m <- length(fitMI$FitList)
  if (what == "path") {
    qbar <- rubin_est(fitMI$PathList)
    qbar <- qbar[which(qbar != 0)]
    vcov_mi <- rubin_vcov(fitMI$PathList, fitMI$PathVCOVList)
    dfcom <- fitMI$nobs - length(qbar)
  }
  else if (what == "loading") {
    qbar <- rubin_est(fitMI$LoadingList)
    qbar <- qbar[which(qbar != 0)]
    vcov_mi <- rubin_vcov(fitMI$LoadingList, fitMI$LoadingVCOVList)
    dfcom <- fitMI$nobs - 1
  }
  ubar <- diag(vcov_mi$W)
  b <- diag(vcov_mi$B)
  t <- diag(vcov_mi$T)  # same as ubar + (1 + 1/m) * b
  lambda <- (1 + 1 / m) * b / t
  lambda[lambda < 1e-04] <- 1e-04
  dfold <- (m - 1) / lambda^2
  dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
  dfold * dfobs / (dfold + dfobs)
}

poolMI <- function(fitMI, boot_mi, level = 0.95) {
  if (boot_mi == "miboot") {
    path_est <- rubin_est(fitMI$PathList)
    idx <- which(path_est != 0)
    path_est <- path_est[idx]
    path_vcov <- rubin_vcov(fitMI$PathList, fitMI$PathVCOVList)$T
    path_sd <- sqrt(diag(path_vcov))
    df <- barnard_rubin(fitMI, what = "path")
    tval <- stats::qt(1 - (1 - level)/2, df = df)
    path_lwr <- path_est - tval*path_sd
    path_upr <- path_est + tval*path_sd

    load_est <- rubin_est(fitMI$LoadingList)
    idx <- which(load_est != 0)
    load_est <- load_est[idx]
    load_vcov <- rubin_vcov(fitMI$LoadingList, fitMI$LoadingVCOVList)$T
    load_sd <- sqrt(diag(load_vcov))
    df <- barnard_rubin(fitMI, what = "loading")
    tval <- stats::qt(1 - (1 - level)/2, df = df)
    load_lwr <- load_est - tval*load_sd
    load_upr <- load_est + tval*load_sd
  }
  else if (boot_mi == "bootmi") {
    npath <- sum(cSEM::parseModel(fitMI$model)$structural)
    path <- fitMI$BootMatrix[, 1:npath]
    load <- fitMI$BootMatrix[, (npath + 1):ncol(fitMI$BootMatrix)]
    alpha <- (1 - level)

    path_est <- colMeans(path)
    path_vcov <- var(path)
    path_sd <- sqrt(diag(path_vcov))
    path_ci <- apply(path, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
    path_lwr <- path_ci[1, ]
    path_upr <- path_ci[2, ]

    load_est <- colMeans(load)
    load_vcov <- var(load)
    load_sd <- sqrt(diag(load_vcov))
    load_ci <- apply(load, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
    load_lwr <- load_ci[1, ]
    load_upr <- load_ci[2, ]
  }
  else if (boot_mi == "miboot_pooled") {
    npath <- sum(cSEM::parseModel(fitMI$model)$structural)
    path <- fitMI$PooledSample[, 1:npath]
    load <- fitMI$PooledSample[, (npath + 1):ncol(fitMI$PooledSample)]
    alpha <- (1 - level)

    path_est <- colMeans(path)
    path_vcov <- var(path)
    path_sd <- sqrt(diag(path_vcov))
    path_ci <- apply(path, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
    path_lwr <- path_ci[1, ]
    path_upr <- path_ci[2, ]

    load_est <- colMeans(load)
    load_vcov <- var(load)
    load_sd <- sqrt(diag(load_vcov))
    load_ci <- apply(load, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
    load_lwr <- load_ci[1, ]
    load_upr <- load_ci[2, ]
  }
  else if (boot_mi == "bootmi_pooled") {
    npath <- sum(cSEM::parseModel(fitMI$model)$structural)
    path <- fitMI$PooledSample[, 1:npath]
    load <- fitMI$PooledSample[, (npath + 1):ncol(fitMI$PooledSample)]
    alpha <- (1 - level)

    path_est <- colMeans(path)
    path_vcov <- var(path)
    path_sd <- sqrt(diag(path_vcov))
    path_ci <- apply(path, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
    path_lwr <- path_ci[1, ]
    path_upr <- path_ci[2, ]

    load_est <- colMeans(load)
    load_vcov <- var(load)
    load_sd <- sqrt(diag(load_vcov))
    load_ci <- apply(load, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
    load_lwr <- load_ci[1, ]
    load_upr <- load_ci[2, ]
  }

  out <- list(path_est = path_est, path_vcov = path_vcov, path_sd = path_sd,
              path_lower = path_lwr, path_upper = path_upr,
              load_est = load_est, load_vcov = load_vcov, load_sd = load_sd,
              load_lower = load_lwr, load_upper = load_upr)
  out
}

csem_combine <- function(csem_res, level = 0.95) {
  alpha <- 1 - level

  path_est <- colSums(csem_res$Estimates$Path_estimates)
  path_est <- path_est[which(path_est != 0)]
  path_boot <- csem_res$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled
  path_vcov <- cov(path_boot)
  path_sd <- sqrt(diag(path_vcov))
  path_ci <- apply(path_boot, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
  path_lwr <- path_ci[1, ]
  path_upr <- path_ci[2, ]
  
  load_est <- colSums(csem_res$Estimates$Loading_estimates)
  load_boot <- csem_res$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled
  load_vcov <- cov(load_boot)
  load_sd <- sqrt(diag(load_vcov))
  load_ci <- apply(load_boot, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
  load_lwr <- load_ci[1, ]
  load_upr <- load_ci[2, ]

  out <- list(path_est = path_est, path_vcov = path_vcov, path_sd = path_sd,
              path_lower = path_lwr, path_upper = path_upr,
              load_est = load_est, load_vcov = load_vcov, load_sd = load_sd,
              load_lower = load_lwr, load_upper = load_upr)
  out
}

pool_samples <- function(fitMI) {
  pooled_sample <- cbind(fitMI$FitList[[1]]$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled,
                         fitMI$FitList[[1]]$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled)
  for (i in 2:length(fitMI$FitList)) {
    tmp <- cbind(fitMI$FitList[[i]]$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled,
                 fitMI$FitList[[i]]$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled)
    pooled_sample <- rbind(pooled_sample, tmp)
  }

  pooled_sample
}

rawbias <- function(estimates, true) {
  rowMeans(estimates) - true
}

percentbias <- function(estimates, true) {
  100*rawbias(estimates, true)/true
}

coveragerate <- function(ciestimates, true) {
  lower <- ciestimates[, , 1]
  upper <- ciestimates[, , 2]
  rowMeans(lower < true & true < upper)
}

averagewidth <- function(ciestimates) {
  rowMeans(ciestimates[, , 2] - ciestimates[, , 1])
}

rmse <- function(estimates, true) {
  sqrt(rowMeans((estimates - true)^2))
}

qualitymeasures <- function(estimates, true) {
  data.frame(RB = rawbias(estimates, true),
    PB = percentbias(estimates, true),
    CR = coveragerate(estimates, true),
    AW = averagewidth(estimates),
    RMSE = rmse(estimates, true))
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
