extract_results <- function(res, approach, type = "path", what = "est", simplify = TRUE) {
  if (!(type == "path" | type == "load")) {
    stop("the requested results are not available.")
  }
  if (!(what == "est" | what == "vcov" | what == "sd" | what == "ci")) {
    stop("the requested results are not available.")
  }
  if (what == "vcov" && simplify) {
    simplify <- FALSE
  }

  if (what == "ci") {
    get_str <- paste0(type, "_", c("lower", "upper"))
    npar <- length(res[[1]][[approach]][[get_str[1]]])
    out <- array(NA, dim = c(npar, length(res) - 3, 2))
    for (i in 1:length(get_str)) {
      out[, , i] <- sapply(res[1:(length(res) - 3)],
        function(x, method, str) x[[method]][[str[i]]],
        method = approach, str = get_str,
        simplify = simplify, USE.NAMES = TRUE)
    }
    dimnames(out)[[2]] <- paste0("run_", 1:(length(res) - 3))
    dimnames(out)[[3]] <- c("lower", "upper")
  }
  else {
    get_str <- paste0(type, "_", what)
    out <- sapply(res[1:(length(res) - 3)],
      function(x, method, str) x[[method]][[str]],
      method = approach, str = get_str,
      simplify = simplify, USE.NAMES = TRUE)
    # out <- as.data.frame(out)
  }

  out
}

aggregate_results <- function(res, true_coefs = NULL, methods = "all", qual_meas = "all") {
  if (is.null(true_coefs) | missing(true_coefs) | any(is.na(true_coefs))) {
    true_coefs <- NA
    warning("the quality measures are not available.")
  }

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
  qualnm <- c("RB", "PB", "CR", "AW", "RMSE")
  if (is.null(qual_meas) | any(qual_meas == "all")) {
    qual_meas <- qualnm
  }
  else {
    if (!all(qual_meas %in% qualnm))
      stop("the specified quality measures are not available.")
  }

  out <- data.frame(true = true_coefs)
  for (meth in methods) {
    res_path_est <- extract_results(res, approach = meth, type = "path", what = "est")
    res_path_sd <- extract_results(res, approach = meth, type = "path", what = "sd")
    res_path_ci <- extract_results(res, approach = meth, type = "path", what = "ci")
    res_load_est <- extract_results(res, approach = meth, type = "load", what = "est")
    res_load_sd <- extract_results(res, approach = meth, type = "load", what = "sd")
    res_load_ci <- extract_results(res, approach = meth, type = "load", what = "ci")
    res_est <- rbind(res_path_est, res_load_est)
    res_sd <- rbind(res_path_sd, res_load_sd)
    res_ci <- array(NA, dim = c(dim(res_est), 2))
    res_ci[, , 1] <- rbind(res_path_ci[, , 1], res_load_ci[, , 1])
    res_ci[, , 2] <- rbind(res_path_ci[, , 2], res_load_ci[, , 2])
    # if (length(true_coefs) != length(res_est))
    #   stop("the number of true coefficient values does not match with that of the results.")
    out_cn <- c(colnames(out), paste0(meth, c("_est", "_sd", "_lwr", "_upr")))
    out <- cbind(out, rowMeans(res_est), rowMeans(res_sd), rowMeans(res_ci[, , 1]), rowMeans(res_ci[, , 2]))
    colnames(out) <- out_cn

    if (!is.null(true_coefs) && !missing(true_coefs) && !any(is.na(true_coefs))) {
      out_cn <- c(colnames(out), paste0(qual_meas, "_", meth))
      for (qual in qual_meas) {
        qual_res <- switch(qual,
          RB = rawbias(res_est, true_coefs),
          PB = percentbias(res_est, true_coefs),
          CR = coveragerate(res_ci, true_coefs),
          AW = averagewidth(res_ci),
          RMSE = rmse(res_est, true_coefs))
        out <- cbind(out, qual_res)
      }
      colnames(out) <- out_cn
    }
  }
  attr(out, "methods") <- methods
  if (!is.null(true_coefs) && !missing(true_coefs) && !any(is.na(true_coefs))) {
    attr(out, "qual_meas") <- qual_meas
  }

  out
}

extract_data <- function(res, full = TRUE) {
  out <- list()
  for (i in 1:(length(res) - 3)) {
    if (full) {
      out[[i]] <- res[[i]][["dat_orig"]]
    }
    else {
      out[[i]] <- res[[i]][["dat_miss"]]
    }
  }

  out
}
