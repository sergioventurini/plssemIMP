# CHANGE: added boot_mc_cores parameter (default NULL = use all available cores minus one).
# For MIBOOT the bottleneck is NOT MICE (called only once) but the m × R cSEM bootstrap
# calls that run sequentially when .eval_plan = "sequential".  The fix is to parallelise
# across the m imputed datasets: each dataset gets its own cSEM call in a separate fork,
# all running concurrently.  The number of active forks is capped at min(boot_mc_cores, m)
# because there is no benefit in allocating more cores than datasets.
plssemMIBOOT <- function(model, data, ..., m = 5, miArgs = list(),
                         csemArgs = list(), miPackage = "mice",
                         verbose = FALSE, seed = NULL, level = 0.95,
                         boot_mc_cores = NULL) {    # <<< NEW PARAMETER
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }
  imputedData <- NULL

  if (missing(data)) {
    stop("a dataset is needed to run the plssemMIBOOT() function.")
  }

  # CHANGE: resolve effective core count, then cap at m (no gain beyond one core per dataset).
  effective_boot_cores <- if (is.null(boot_mc_cores)) {    # <<< NEW
    max(1L, parallel::detectCores() - 1L)                  # <<< NEW
  } else {                                                  # <<< NEW
    max(1L, as.integer(boot_mc_cores))                     # <<< NEW
  }                                                         # <<< NEW

  if (any(is.na(data))) {
    if (miPackage[1] == "Amelia") {
      requireNamespace("Amelia")
      imputeCall <- c(list(Amelia::amelia, x = data, m = m, p2s = 0), miArgs)
      # imputedData <- unclass(eval(as.call(imputeCall))$imputations)
      imputedData <- unclass(do.call(Amelia::amelia, c(list(x = data, m = m, p2s = 0), miArgs))$imputations)
    }
    else if (miPackage[1] == "mice") {
      requireNamespace("mice")
      pred_matrix <- mice::quickpred(data, mincor = 0.1, minpuc = 0.5)
      imputeCall <- c(list(mice::mice, data = data, m = m,
                           diagnostics = FALSE, printFlag = FALSE), miArgs)
      # miceOut <- eval(as.call(imputeCall))
      miceOut <- do.call(mice::mice, c(list(data = data, m = m,
                         diagnostics = FALSE, printFlag = FALSE), miArgs,
                         predictorMatrix = pred_matrix))
      imputedData <- lapply(as.list(1:m), function(i) mice::complete(data = miceOut,
                            action = i, include = FALSE))
    }
    else stop("currently plssemMIBOOT() only supports imputation by Amelia or mice.")
  }
  else {
    imputedData <- list(data)
  }

  csemListCall <- list(cSEM::csem, .model = model, .data = imputedData)
  csemListCall <- c(csemListCall, csemArgs)
  fit <- list()

  # CHANGE: parallelise across the m imputed datasets instead of calling csem with the full
  # list.  Each individual call forces .eval_plan = "sequential" so that cSEM does not try
  # to fork again inside an already-forked mclapply worker (prevents nested parallelism).
  # A list of m individual cSEMResults objects is structurally identical to iterating over
  # a cSEMResults_multi, so all downstream lapply() calls are unchanged.
  csemArgs_seq <- csemArgs                                       # <<< NEW
  csemArgs_seq$.eval_plan <- "sequential"                        # <<< NEW
  n_cores_for_m <- min(effective_boot_cores, length(imputedData))  # <<< NEW: cap at m
  fit$FitList <- parallel::mclapply(                             # <<< CHANGED
    imputedData,                                                 # <<< CHANGED
    function(d) {                                                # <<< CHANGED
      suppressWarnings(do.call(cSEM::csem,                       # <<< CHANGED
        c(list(.model = model, .data = d), csemArgs_seq)))       # <<< CHANGED
    },                                                           # <<< CHANGED
    mc.cores = n_cores_for_m,                                    # <<< CHANGED
    mc.preschedule = FALSE                                       # <<< CHANGED
  )                                                              # <<< CHANGED

  fit$PathList <- lapply(fit$FitList,
    function(x) {
      path <- x$Estimates$Path_estimates
      path <- matrix2vec(path, names = TRUE)
      rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
      path[path != 0, ]
    })
  fit$LoadingList <- lapply(fit$FitList,
    function(x) {
      colSums(x$Estimates$Loading_estimates)
    })
  fit$PathVCOVList <- lapply(fit$FitList,
    function(x) {
      out <- cov(x$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled)
      rownames(out) <- gsub("eta", "gamma", gsub(".~ eta", "", rownames(out)))
      colnames(out) <- rownames(out)
      out
    })
  fit$LoadingVCOVList <- lapply(fit$FitList,
    function(x) {
      out <- cov(x$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled)
      rownames(out) <- gsub(".* ", "", rownames(out))
      colnames(out) <- rownames(out)
      out
    })
  fit$DataList <- imputedData
  fit$model <- model
  fit$call <- CALL
  fit$csemListCall <- csemListCall
  if (any(is.na(data))) fit$imputeCall <- imputeCall
  fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
  if (!all(unlist(fit$convList)))
    warning("the model did not converge for any imputed data sets.")
  fit$dots <- dots
  fit$nobs <- nrow(data)
  fit$pooled <- poolMI(fit, boot_mi = "miboot", level = level)
  fit$pooled
}

# CHANGE: added boot_mc_cores parameter (default NULL = use all available cores minus one).
# This allows the caller (run_sims) to restrict the number of cores used by the inner
# bootstrap loop, enabling explicit two-level parallelism: runs_parallel outer runs each
# using boot_mc_cores cores for their bootstrap iterations.
plssemBOOTMI <- function(model, data, ..., m = 5, miArgs = list(),
                         csemArgs = list(), miPackage = "mice",
                         verbose = FALSE, seed = NULL, level = 0.95,
                         boot_mc_cores = NULL) {    # <<< NEW PARAMETER
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }

  if (missing(data)) {
    stop("a dataset is needed to run the plssemBOOTMI() function.")
  }

  # CHANGE: resolve effective core count once, using the passed value or falling back
  # to all-cores-minus-one when called standalone (boot_mc_cores = NULL).
  effective_boot_cores <- if (is.null(boot_mc_cores)) {   # <<< NEW
    max(1L, parallel::detectCores() - 1L)                 # <<< NEW
  } else {                                                 # <<< NEW
    max(1L, as.integer(boot_mc_cores))                    # <<< NEW
  }                                                        # <<< NEW

  bootstrap_mi <- function(data, indices, mipkg, miargs, miruns, csemmodel, csemargs, verb) {
    boot_sample <- data[indices, ]

    imputedData <- NULL
    if (any(is.na(data))) {
      if (mipkg[1] == "Amelia") {
        requireNamespace("Amelia")
        imputeCall <- c(list(Amelia::amelia, x = boot_sample, m = miruns, p2s = 0), miargs)
        # imputedData <- unclass(eval(as.call(imputeCall))$imputations)
        imputedData <- unclass(do.call(Amelia::amelia, c(list(x = boot_sample, m = miruns, p2s = 0), miargs))$imputations)
      }
      else if (mipkg[1] == "mice") {
        requireNamespace("mice")
        pred_matrix <- pred_matrix_orig
        imputeCall <- c(list(mice::mice, data = boot_sample, m = miruns,
                             diagnostics = FALSE, printFlag = FALSE), miargs)
        # miceOut <- eval(as.call(imputeCall))
        miceOut <- do.call(mice::mice, c(list(data = boot_sample, m = miruns,
                           diagnostics = FALSE, printFlag = FALSE), miargs,
                           predictorMatrix = pred_matrix))
        imputedData <- lapply(as.list(1:miruns), function(i) mice::complete(data = miceOut,
                              action = i, include = FALSE))
      }
      else stop("currently plssemBOOTMI() only supports imputation by Amelia or mice.")
    }
    else {
      imputedData <- list(data)
    }

    csemListCall <- list(cSEM::csem, .model = csemmodel, .data = imputedData)
    csemListCall <- c(csemListCall, csemargs)
    fit <- list()
    # fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
    fit$FitList <- suppressWarnings(do.call(cSEM::csem, c(list(.model = csemmodel, .data = imputedData), csemargs)))
    fit$PathList <- lapply(fit$FitList,
      function(x) {
        path <- x$Estimates$Path_estimates
        path <- matrix2vec(path, names = TRUE)
        rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
        path[path != 0, ]
      })
    fit$LoadingList <- lapply(fit$FitList,
      function(x) {
        colSums(x$Estimates$Loading_estimates)
      })
    fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
    if (!all(unlist(fit$convList)))
      warning("the model did not converge for any imputed data sets.")
    fit$DataList <- imputedData
    c(rubin_est(fit$PathList), rubin_est(fit$LoadingList))
  }

  boot_fit <- list()

  # generate bootstrap indices
  n <- nrow(data)
  R <- csemArgs$.R
  indices_list <- rbind(
    1:n,
    t(replicate(R, sample.int(n, replace = TRUE)))
  )

  pred_matrix_orig <- if (miPackage[1] == "mice" && any(is.na(data))) {
    mice::quickpred(data, mincor = 0.1, minpuc = 0.5)
  } else NULL
  # CHANGE: use effective_boot_cores instead of hardcoded detectCores() - 1
  bootRes <- parallel::mclapply(seq_len(nrow(indices_list)),      ## WARNING: this only works on unix-based OS
                                function(r) {
                                  idx <- indices_list[r, ]
                                  bootstrap_mi(data = data, indices = idx, mipkg = miPackage,
                                               miargs = miArgs, miruns = m, csemmodel = model,
                                               csemargs = csemArgs, verb = verbose)
                                },
                                mc.cores = effective_boot_cores,    # <<< CHANGED
                                mc.preschedule = FALSE              # <<< CHANGED
  )

  boot_fit$BootOrig <- bootRes[[1]]   # these are the cSEM results on the original dataset after imputation
  boot_fit$BootMatrix <- do.call(rbind, bootRes[2:length(bootRes)])

  boot_fit$model <- model
  boot_fit$call <- CALL
  boot_fit$dots <- dots
  boot_fit$nobs <- nrow(data)
  boot_fit$pooled <- poolMI(boot_fit, boot_mi = "bootmi", level = level)
  boot_fit$pooled
}

# CHANGE: added boot_mc_cores parameter (same rationale as plssemMIBOOT above).
# MIBOOT_PS has the same structure as MIBOOT: MICE once, then cSEM bootstraps each of the
# m datasets independently.  Parallelising over m datasets gives the same ~m× speedup.
plssemMIBOOT_PS <- function(model, data, ..., m = 5, miArgs = list(),
                            csemArgs = list(), miPackage = "mice",
                            verbose = FALSE, seed = NULL, level = 0.95,
                            boot_mc_cores = NULL) {    # <<< NEW PARAMETER
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }
  imputedData <- NULL

  if (missing(data)) {
    stop("a dataset is needed to run the plssemMIBOOT_PS() function.")
  }

  # CHANGE: resolve effective core count, then cap at m.
  effective_boot_cores <- if (is.null(boot_mc_cores)) {    # <<< NEW
    max(1L, parallel::detectCores() - 1L)                  # <<< NEW
  } else {                                                  # <<< NEW
    max(1L, as.integer(boot_mc_cores))                     # <<< NEW
  }                                                         # <<< NEW

  if (any(is.na(data))) {
    if (miPackage[1] == "Amelia") {
      requireNamespace("Amelia")
      imputeCall <- c(list(Amelia::amelia, x = data, m = m, p2s = 0), miArgs)
      # imputedData <- unclass(eval(as.call(imputeCall))$imputations)
      imputedData <- unclass(do.call(Amelia::amelia, c(list(x = data, m = m, p2s = 0), miArgs))$imputations)
    }
    else if (miPackage[1] == "mice") {
      requireNamespace("mice")
      pred_matrix <- mice::quickpred(data, mincor = 0.1, minpuc = 0.5)
      imputeCall <- c(list(mice::mice, data = data, m = m,
                           diagnostics = FALSE, printFlag = FALSE), miArgs)
      # miceOut <- eval(as.call(imputeCall))
      miceOut <- do.call(mice::mice, c(list(data = data, m = m,
                         diagnostics = FALSE, printFlag = FALSE), miArgs,
                         predictorMatrix = pred_matrix))
      imputedData <- lapply(as.list(1:m), function(i) mice::complete(data = miceOut,
                            action = i, include = FALSE))
    }
    else stop("currently plssemMIBOOT_PS() only supports imputation by Amelia or mice.")
  }
  else {
    imputedData <- list(data)
  }

  csemListCall <- list(cSEM::csem, .model = model, .data = imputedData)
  csemListCall <- c(csemListCall, csemArgs)
  fit <- list()

  # CHANGE: same as plssemMIBOOT — parallelise the m independent cSEM calls with mclapply,
  # forcing .eval_plan = "sequential" inside each worker to prevent nested forking.
  csemArgs_seq <- csemArgs                                       # <<< NEW
  csemArgs_seq$.eval_plan <- "sequential"                        # <<< NEW
  n_cores_for_m <- min(effective_boot_cores, length(imputedData))  # <<< NEW: cap at m
  fit$FitList <- parallel::mclapply(                             # <<< CHANGED
    imputedData,                                                 # <<< CHANGED
    function(d) {                                                # <<< CHANGED
      suppressWarnings(do.call(cSEM::csem,                       # <<< CHANGED
        c(list(.model = model, .data = d), csemArgs_seq)))       # <<< CHANGED
    },                                                           # <<< CHANGED
    mc.cores = n_cores_for_m,                                    # <<< CHANGED
    mc.preschedule = FALSE                                       # <<< CHANGED
  )                                                              # <<< CHANGED

  fit$PathList <- lapply(fit$FitList,
    function(x) {
      path <- x$Estimates$Path_estimates
      path <- matrix2vec(path, names = TRUE)
      rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
      path[path != 0, ]
    })
  fit$LoadingList <- lapply(fit$FitList,
    function(x) {
      colSums(x$Estimates$Loading_estimates)
    })
  fit$PooledSample <- pool_samples(fit)
  fit$DataList <- imputedData
  fit$model <- model
  fit$call <- CALL
  fit$csemListCall <- csemListCall
  if (any(is.na(data))) fit$imputeCall <- imputeCall
  fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
  if (!all(unlist(fit$convList)))
    warning("the model did not converge for any imputed data sets.")
  fit$dots <- dots
  fit$nobs <- nrow(data)
  fit$pooled <- poolMI(fit, boot_mi = "miboot_pooled", level = level)
  fit$pooled
}

# CHANGE: added boot_mc_cores parameter (same rationale as plssemBOOTMI above).
plssemBOOTMI_PS <- function(model, data, ..., m = 5, miArgs = list(),
                            csemArgs = list(), miPackage = "mice",
                            verbose = FALSE, seed = NULL, level = 0.95,
                            boot_mc_cores = NULL) {    # <<< NEW PARAMETER
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }

  if (missing(data)) {
    stop("a dataset is needed to run the plssemBOOTMI_PS() function.")
  }

  # CHANGE: resolve effective core count
  effective_boot_cores <- if (is.null(boot_mc_cores)) {   # <<< NEW
    max(1L, parallel::detectCores() - 1L)                 # <<< NEW
  } else {                                                 # <<< NEW
    max(1L, as.integer(boot_mc_cores))                    # <<< NEW
  }                                                        # <<< NEW

  bootstrap_mi_ps <- function(data, indices, mipkg, miargs, miruns, csemmodel, csemargs, verb) {
    boot_sample <- data[indices, ]

    imputedData <- NULL
    if (any(is.na(data))) {
      if (mipkg[1] == "Amelia") {
        requireNamespace("Amelia")
        imputeCall <- c(list(Amelia::amelia, x = boot_sample, m = miruns, p2s = 0), miargs)
        # imputedData <- unclass(eval(as.call(imputeCall))$imputations)
        imputedData <- unclass(do.call(Amelia::amelia, c(list(x = boot_sample, m = miruns, p2s = 0), miargs))$imputations)
      }
      else if (mipkg[1] == "mice") {
        requireNamespace("mice")
        pred_matrix <- mice::quickpred(boot_sample, mincor = 0.1, minpuc = 0.5)
        imputeCall <- c(list(mice::mice, data = boot_sample, m = miruns,
                             diagnostics = FALSE, printFlag = FALSE), miargs)
        # miceOut <- eval(as.call(imputeCall))
        miceOut <- do.call(mice::mice, c(list(data = boot_sample, m = miruns,
                           diagnostics = FALSE, printFlag = FALSE), miargs,
                           predictorMatrix = pred_matrix))
        imputedData <- lapply(as.list(1:miruns), function(i) mice::complete(data = miceOut,
                              action = i, include = FALSE))
      }
      else stop("currently plssemBOOTMI_PS() only supports imputation by Amelia or mice.")
    }
    else {
      imputedData <- list(data)
    }

    csemListCall <- list(cSEM::csem, .model = csemmodel, .data = imputedData)
    csemListCall <- c(csemListCall, csemargs)
    fit <- list()
    # fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
    fit$FitList <- suppressWarnings(do.call(cSEM::csem, c(list(.model = csemmodel, .data = imputedData), csemargs)))
    fit$ParamList <- lapply(fit$FitList,
      function(x) {
        path <- x$Estimates$Path_estimates
        path <- matrix2vec(path, names = TRUE)
        rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
        path <- path[path != 0, ]
        load <- colSums(x$Estimates$Loading_estimates)
        c(path, load)
      })
    fit$PooledMI <- do.call(cbind, fit$ParamList)
    fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
    if (!all(unlist(fit$convList)))
      warning("the model did not converge for any imputed data sets.")
    fit$DataList <- imputedData
    fit$PooledMI
  }

  boot_fit <- list()

  # generate bootstrap indices
  n <- nrow(data)
  R <- csemArgs$.R
  indices_list <- rbind(
    1:n,
    t(replicate(R, sample.int(n, replace = TRUE)))
  )

  # CHANGE: use effective_boot_cores instead of hardcoded detectCores() - 1
  bootRes <- parallel::mclapply(seq_len(nrow(indices_list)),      ## WARNING: this only works on unix-based OS
                                function(r) {
                                  idx <- indices_list[r, ]
                                  bootstrap_mi_ps(data = data, indices = idx, mipkg = miPackage,
                                                  miargs = miArgs, miruns = m, csemmodel = model,
                                                  csemargs = csemArgs, verb = verbose)
                                },
                                mc.cores = effective_boot_cores,   # <<< CHANGED
                                mc.preschedule = FALSE             # <<< CHANGED
  )

  boot_fit$BootOrig <- bootRes[[1]]   # these are the cSEM results on the original dataset after imputation
  boot_fit$BootMatrix <- t(sapply(bootRes[2:length(bootRes)],
                                  function(x) as.vector(x),
                                  simplify = TRUE))

  npars <- ncol(boot_fit$BootMatrix)/m
  bootmi_mat <- matrix(NA, nrow = 0, ncol = npars)
  for (i in 1:m) {
    bootmi_mat <- rbind(bootmi_mat, boot_fit$BootMatrix[, ((i - 1)*npars + 1):(npars*i)])
  }
  boot_fit$PooledSample <- bootmi_mat
  boot_fit$model <- model
  boot_fit$call <- CALL
  boot_fit$dots <- dots
  boot_fit$nobs <- nrow(data)
  boot_fit$pooled <- poolMI(boot_fit, boot_mi = "bootmi_pooled", level = level)
  boot_fit$pooled
}

# CHANGE: added boot_mc_cores parameter (same rationale as plssemBOOTMI above).
plssemWGT_BOOTMI <- function(model, data, ..., m = 5, miArgs = list(),
                             csemArgs = list(), miPackage = "mice",
                             verbose = FALSE, seed = NULL, level = 0.95, wgtType = "rows",
                             boot_mc_cores = NULL) {    # <<< NEW PARAMETER
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }

  if (missing(data)) {
    stop("a dataset is needed to run the plssemWGT_BOOTMI() function.")
  }

  # CHANGE: resolve effective core count
  effective_boot_cores <- if (is.null(boot_mc_cores)) {   # <<< NEW
    max(1L, parallel::detectCores() - 1L)                 # <<< NEW
  } else {                                                 # <<< NEW
    max(1L, as.integer(boot_mc_cores))                    # <<< NEW
  }                                                        # <<< NEW

  w_bootstrap_mi <- function(data, indices, mipkg, miargs, miruns, csemmodel, csemargs, verb, wtype) {
    boot_sample <- data[indices, ]
    if (wtype == "rows") {
      wgt <- sum(apply(boot_sample, 1, function(x) any(is.na(x))))  # num of rows with at least one missing
    }
    else if (wtype == "all") {
      wgt <- sum(is.na(boot_sample))                                # num of missing data cells
    }

    imputedData <- NULL
    if (any(is.na(data))) {
      if (mipkg[1] == "Amelia") {
        requireNamespace("Amelia")
        imputeCall <- c(list(Amelia::amelia, x = boot_sample, m = miruns, p2s = 0), miargs)
        # imputedData <- unclass(eval(as.call(imputeCall))$imputations)
        imputedData <- unclass(do.call(Amelia::amelia, c(list(x = boot_sample, m = miruns, p2s = 0), miargs))$imputations)
      }
      else if (mipkg[1] == "mice") {
        requireNamespace("mice")
        pred_matrix <- pred_matrix_orig
        imputeCall <- c(list(mice::mice, data = boot_sample, m = miruns,
                             diagnostics = FALSE, printFlag = FALSE), miargs)
        # miceOut <- eval(as.call(imputeCall))
        miceOut <- do.call(mice::mice, c(list(data = boot_sample, m = miruns,
                           diagnostics = FALSE, printFlag = FALSE), miArgs,
                           predictorMatrix = pred_matrix))
        imputedData <- lapply(as.list(1:miruns), function(i) mice::complete(data = miceOut,
                              action = i, include = FALSE))
      }
      else stop("currently plssemWGT_BOOTMI() only supports imputation by Amelia or mice.")
    }
    else {
      imputedData <- list(data)
    }

    csemListCall <- list(cSEM::csem, .model = csemmodel, .data = imputedData)
    csemListCall <- c(csemListCall, csemargs)
    fit <- list()
    # fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
    fit$FitList <- suppressWarnings(do.call(cSEM::csem, c(list(.model = csemmodel, .data = imputedData), csemargs)))
    fit$PathList <- lapply(fit$FitList,
      function(x) {
        path <- x$Estimates$Path_estimates
        path <- matrix2vec(path, names = TRUE)
        rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
        path[path != 0, ]
      })
    fit$LoadingList <- lapply(fit$FitList,
      function(x) {
        colSums(x$Estimates$Loading_estimates)
      })
    fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
    if (!all(unlist(fit$convList)))
      warning("the model did not converge for any imputed data sets.")
    fit$DataList <- imputedData
    c(wgt, rubin_est(fit$PathList), rubin_est(fit$LoadingList))
  }
  boot_fit <- list()

  # generate bootstrap indices
  n <- nrow(data)
  R <- csemArgs$.R
  indices_list <- rbind(
    1:n,
    t(replicate(R, sample.int(n, replace = TRUE)))
  )

  pred_matrix_orig <- if (miPackage[1] == "mice" && any(is.na(data))) {
    mice::quickpred(data, mincor = 0.1, minpuc = 0.5)
  } else NULL
  # CHANGE: use effective_boot_cores instead of hardcoded detectCores() - 1
  bootRes <- parallel::mclapply(seq_len(nrow(indices_list)),      ## WARNING: this only works on unix-based OS
                                function(r) {
                                  idx <- indices_list[r, ]
                                  w_bootstrap_mi(data = data, indices = idx, mipkg = miPackage,
                                                 miargs = miArgs, miruns = m, csemmodel = model,
                                                 csemargs = csemArgs, verb = verbose, wtype = wgtType)
                                },
                                mc.cores = effective_boot_cores,    # <<< CHANGED
                                mc.preschedule = FALSE              # <<< CHANGED
  )

  boot_fit$BootOrig <- (bootRes[[1]])[-1]   # these are the cSEM results on the original dataset after imputation
  boot_fit$BootMatrix <- (do.call(rbind, bootRes[2:length(bootRes)]))[, -1, drop = FALSE]
  wgts <- (do.call(rbind, bootRes[2:length(bootRes)]))[, 1]

  wgts[wgts == 0] <- 1
  if (wgtType == "rows") {
    wgts <- nrow(data)/wgts
  }
  else if (wgtType == "all") {
    wgts <- prod(dim(data))/wgts
  }

  boot_fit$model <- model
  boot_fit$call <- CALL
  boot_fit$dots <- dots
  boot_fit$nobs <- nrow(data)
  boot_fit$pooled <- poolMI(boot_fit, boot_mi = "weighted_bootmi", level = level, weights = wgts)
  boot_fit$pooled
}
