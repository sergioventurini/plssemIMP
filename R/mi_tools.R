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

poolMI <- function(fitMI, boot_mi, level = 0.95, weights = NULL) {
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
  else if (boot_mi == "weighted_bootmi") {
    weights <- weights/sum(weights)
    npath <- sum(cSEM::parseModel(fitMI$model)$structural)
    path <- fitMI$BootMatrix[, 1:npath]
    load <- fitMI$BootMatrix[, (npath + 1):ncol(fitMI$BootMatrix)]
    alpha <- (1 - level)

    path_est <- apply(path, 2, weighted.mean, w = weights)
    path_vcov <- cov.wt(x = path, wt = weights, cor = FALSE, center = TRUE, method = "unbiased")$cov
    path_sd <- sqrt(diag(path_vcov))
    path_ci <- apply(path, 2, Hmisc::wtd.quantile, probs = c(alpha/2, 1 - alpha/2), weights = weights, normwt = TRUE)
    path_lwr <- path_ci[1, ]
    path_upr <- path_ci[2, ]

    load_est <- apply(load, 2, weighted.mean, w = weights)
    load_vcov <- cov.wt(x = load, wt = weights, cor = FALSE, center = TRUE, method = "unbiased")$cov
    load_sd <- sqrt(diag(load_vcov))
    load_ci <- apply(load, 2, Hmisc::wtd.quantile, probs = c(alpha/2, 1 - alpha/2), weights = weights, normwt = TRUE)
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

  path_est <- model_coef_vec(csem_res, estimates = TRUE, what = "path")
  path_boot <- csem_res$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled
  path_vcov <- cov(path_boot)
  path_sd <- sqrt(diag(path_vcov))
  path_ci <- apply(path_boot, 2, quantile, probs = c(alpha/2, 1 - alpha/2))
  path_lwr <- path_ci[1, ]
  path_upr <- path_ci[2, ]
  
  load_est <- model_coef_vec(csem_res, estimates = TRUE, what = "load")
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
