mean_impute <- function(X) {
  Ximp <- lapply(X, function(x) {
    if (is.numeric(x)) {
      x[is.na(x)] <- mean(x[!is.na(x)])
    }
    else {
      x[is.na(x)] <- mode_function(x[!is.na(x)])
    }
    x
  })
  Ximp <- as.data.frame(Ximp)

  Ximp
}

meanimp <- function(model, data, csemArgs = list(), verbose = FALSE, level = 0.95, ...) {
  CALL <- match.call()
  dots <- list(...)
  imputedData <- NULL

  if (missing(data)) {
    stop("a data set is needed to run the mean_impute() function.")
  }

  imputedData <-  data
  if (any(is.na(data))) {
    imputedData <- mean_impute(data)
  }

  csemListCall <- list(cSEM::csem, .model = model, .data = imputedData)
  csemListCall <- c(csemListCall, csemArgs)
  fit <- list()
  # fit$FitList <- list(suppressWarnings(eval(as.call(csemListCall))))
  fit$FitList <- list(suppressWarnings(do.call(cSEM::csem, c(list(.model = model, .data = imputedData), csemArgs))))
  fit$PathList <- lapply(fit$FitList, model_coef_vec, estimates = TRUE, what = "path")
  fit$LoadingList <- lapply(fit$FitList, model_coef_vec, estimates = TRUE, what = "load")
  fit$PathVCOVList <- lapply(fit$FitList,
    function(x) {
      cov(x$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled)
    })
  fit$LoadingVCOVList <- lapply(fit$FitList,
    function(x) {
      cov(x$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled)
    })
  fit$DataList <- imputedData
  fit$model <- model
  fit$call <- CALL
  fit$csemListCall <- csemListCall
  fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
  if (!all(unlist(fit$convList))) 
    warning("the model did not converge for any imputed data sets.")
  fit$dots <- dots
  fit$nobs <- nrow(data)
  fit$pooled <- csem_combine(fit$FitList[[1]], level = level)
  fit$pooled
}

knn_impute <- function(X, k = 5, ties = TRUE, FUN = median) {
  n <- nrow(X)
  Q <- ncol(X)
 
  Ximp <- X
  complete_idx <- which(rowSums(is.na(X)) == 0)
  N <- length(complete_idx)          # number of complete observations
  incomplete_idx <- setdiff(1:n, complete_idx)
  R <- n - N                         # number of incomplete observations
 
  for (i in 1:R) {
    vars_nonmiss <- which(!is.na(X[incomplete_idx[i], ]))
    tmp_idx <- which(rowSums(is.na(X[, vars_nonmiss, drop = FALSE])) == 0)
    tmp_idx <- cbind(tmp_idx, 1:length(tmp_idx))
    D_idx <- tmp_idx[which(incomplete_idx[i] == tmp_idx[, 1]), 2]
    disti <- euclidean_dist(X[tmp_idx[, 1], vars_nonmiss, drop = FALSE], D_idx)
    disti[D_idx] <- 1e+100
    disti <- cbind(disti, tmp_idx)
    disti <- disti[order(disti[, 1], disti[, 2]), ]
    disti_uniq <- unique(disti[, 1])
    if (ties && (length(disti_uniq) != nrow(disti))) {
      disti_num <- length(disti_uniq)
      for (j in 1:disti_num) {
        disti_j <- which(disti[, 1] == disti_uniq[j])
        tmp <- disti[disti_j, , drop = FALSE]
        disti[disti_j, ] <- tmp[sample(1:nrow(tmp)), ]  # ties randomly broken
      }
    }    
    if (k > nrow(disti))
      stop("value of k too large.")
    Xsub <- X[disti[1:k, 2], ]
    knew <- k
    while (any(colSums(is.na(Xsub)) == knew)) {
      knew <- knew + 1
      Xsub <- X[disti[1:knew, 2], ]
    }
    for (q in 1:Q) {
      if (is.na(X[incomplete_idx[i], q])) {
        if (is.numeric(X[, q])) {
          Ximp[incomplete_idx[i], q] <- do.call(FUN, list(x = Xsub[, q], na.rm = TRUE))
        }
        else {
          Ximp[incomplete_idx[i], q] <- mode_function(Xsub[, q], na.rm = TRUE)
        }
      }
    }
  }

  Ximp
}

knnimp <- function(model, data, csemArgs = list(), k = 5, method = "euclidean",
  verbose = FALSE, level = 0.95, ...) {
  CALL <- match.call()
  dots <- list(...)
  imputedData <- NULL

  if (missing(data)) {
    stop("a data set is needed to run the knnimp() function.")
  }

  imputedData <-  data
  if (any(is.na(data))) {
    imputedData <- knn_impute(data, k = k)
  }

  csemListCall <- list(cSEM::csem, .model = model, .data = imputedData)
  csemListCall <- c(csemListCall, csemArgs)
  fit <- list()
  # fit$FitList <- list(suppressWarnings(eval(as.call(csemListCall))))
  fit$FitList <- list(suppressWarnings(do.call(cSEM::csem, c(list(.model = model, .data = imputedData), csemArgs))))
  fit$PathList <- lapply(fit$FitList, model_coef_vec, estimates = TRUE, what = "path")
  fit$LoadingList <- lapply(fit$FitList, model_coef_vec, estimates = TRUE, what = "load")
  fit$PathVCOVList <- lapply(fit$FitList,
    function(x) {
      cov(x$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled)
    })
  fit$LoadingVCOVList <- lapply(fit$FitList,
    function(x) {
      cov(x$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled)
    })
  fit$DataList <- imputedData
  fit$model <- model
  fit$call <- CALL
  fit$csemListCall <- csemListCall
  fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
  if (!all(unlist(fit$convList))) 
    warning("the model did not converge for any imputed data sets.")
  fit$dots <- dots
  fit$nobs <- nrow(data)
  fit$k <- k
  fit$pooled <- csem_combine(fit$FitList[[1]], level = level)
  fit$pooled
}
