plssemMIBOOT <- function(model, data, ..., m = 5, miArgs = list(),
                         csemArgs = list(), miPackage = "mice",
                         verbose = FALSE, seed = NULL, level = 0.95) {
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }
  imputedData <- NULL

  if (missing(data)) {
    stop("a dataset is needed to run the plssemMIBOOT() function.")
  }

  if (any(is.na(data))) {
    if (miPackage[1] == "Amelia") {
      requireNamespace("Amelia")
      if (!"package:Amelia" %in% search()) 
        attachNamespace("Amelia")
      imputeCall <- c(list(Amelia::amelia, x = data, m = m, p2s = 0), miArgs)
      imputedData <- unclass(eval(as.call(imputeCall))$imputations)
    }
    else if (miPackage[1] == "mice") {
      requireNamespace("mice")
      if (!"package:mice" %in% search()) 
        attachNamespace("mice")
      imputeCall <- c(list(mice::mice, data = data, m = m, 
                           diagnostics = FALSE, printFlag = FALSE), miArgs)
      miceOut <- eval(as.call(imputeCall))
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
  fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
  fit$PathList <- lapply(fit$FitList,
    function(x) {
      path <- x$Estimates$Path_estimates
      path <- matrix2vec(path, names = TRUE)
      rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
      path[path != 0, ]
    })
  fit$LoadingList <- lapply(fit$FitList,
    function(x) {
      load <- colSums(x$Estimates$Loading_estimates)
      names(load) <- gsub("y", "lambda", names(load))
      load
    })
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
  if (any(is.na(data))) fit$imputeCall <- imputeCall
  fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
  if (!all(unlist(fit$convList))) 
    warning("the model did not converge for all imputed data sets.")
  fit$dots <- dots
  fit$nobs <- nrow(data)
  fit$pooled <- poolMI(fit, boot_mi = "miboot", level = level)
  fit$pooled
}

plssemBOOTMI <- function(model, data, ..., m = 5, miArgs = list(),
                         csemArgs = list(), miPackage = "mice",
                         bootArgs = list(), verbose = FALSE, seed = NULL,
                         level = 0.95) {
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }

  if (missing(data)) {
    stop("a dataset is needed to run the plssemBOOTMI() function.")
  }

  bootstrap_mi <- function(data, indices, mipkg, miargs, miruns, csemmodel, csemargs, verb) {
    boot_sample <- data[indices, ]

    imputedData <- NULL
    if (any(is.na(data))) {
      if (mipkg[1] == "Amelia") {
        requireNamespace("Amelia")
        if (!"package:Amelia" %in% search()) 
          attachNamespace("Amelia")
        imputeCall <- c(list(Amelia::amelia, x = boot_sample, m = miruns, p2s = 0), miargs)
        imputedData <- unclass(eval(as.call(imputeCall))$imputations)
      }
      else if (mipkg[1] == "mice") {
        requireNamespace("mice")
        if (!"package:mice" %in% search()) 
          attachNamespace("mice")
        imputeCall <- c(list(mice::mice, data = boot_sample, m = miruns, 
                             diagnostics = FALSE, printFlag = FALSE), miargs)
        miceOut <- eval(as.call(imputeCall))
        imputedData <- lapply(as.list(1:m), function(i) mice::complete(data = miceOut,
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
    fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
    fit$PathList <- lapply(fit$FitList,
      function(x) {
        path <- x$Estimates$Path_estimates
        path <- matrix2vec(path, names = TRUE)
        rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
        path[path != 0, ]
      })
    fit$LoadingList <- lapply(fit$FitList,
      function(x) {
        load <- colSums(x$Estimates$Loading_estimates)
        names(load) <- gsub("y", "lambda", names(load))
        load
      })
    fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
    if (!all(unlist(fit$convList))) 
      warning("the model did not converge for all imputed data sets.")
    fit$DataList <- imputedData
    c(rubin_est(fit$PathList), rubin_est(fit$LoadingList))
  }

  if (!is.null(bootArgs) && !missing(bootArgs)) {
    if (bootArgs$parallel != "no" && bootArgs$ncpus > 1)
      verbose <- FALSE
  }
  boot_fit <- list()

  ###
  # new implementation
  # fictitious step only needed to generate the bootstrap indices
  bootListCall <- list(boot::boot, data = data, statistic = function(d, i) 0, R = csemArgs$.R)
  bootListCall <- c(bootListCall, bootArgs)
  bootRes <- eval(as.call(bootListCall))
  indices_list <- rbind(1:nrow(data), boot::boot.array(bootRes, indices = TRUE))
  if (.Platform$OS.type == "unix") {
    bootRes <- parallel::mclapply(seq_len(nrow(indices_list)),  # this can be used only with macOS computers
                                  function(r) {
                                    idx <- indices_list[r, ]
                                    bootstrap_mi(data = data, indices = idx, mipkg = miPackage,
                                                 miargs = miArgs, miruns = m, csemmodel = model,
                                                 csemargs = csemArgs, verb = verbose)
                                  },
                                  mc.cores = argsBOOT$ncpus)
  }
  else {
    bootRes <- lapply(seq_len(nrow(indices_list)),  # under Windows we should use parLapply() but it is too cumbersome
                      function(r) {
                        idx <- indices_list[r, ]
                        bootstrap_mi(data = data, indices = idx, mipkg = miPackage,
                                     miargs = miArgs, miruns = m, csemmodel = model,
                                     csemargs = csemArgs, verb = verbose)
                        })
  }
  boot_fit$BootOrig <- bootRes[[1]]   # these are the cSEM results on the original dataset after imputation
  boot_fit$BootMatrix <- do.call(rbind, bootRes[2:length(bootRes)])
  ###

  ###
  # old implementation (to delete at a certain point)
  # bootListCall <- list(boot::boot, data = data, statistic = bootstrap_mi,
  #   R = csemArgs$.R, mipkg = miPackage, miargs = miArgs, miruns = m,
  #   csemmodel = model, csemargs = csemArgs, verb = verbose)
  # bootListCall <- c(bootListCall, bootArgs)
  # bootRes <- eval(as.call(bootListCall))
  # boot_fit$BootOrig <- bootRes$t0   # these are the cSEM results on the original dataset after imputation
  # boot_fit$BootMatrix <- bootRes$t
  ###

  boot_fit$model <- model
  boot_fit$call <- CALL
  boot_fit$dots <- dots
  boot_fit$nobs <- nrow(data)
  boot_fit$pooled <- poolMI(boot_fit, boot_mi = "bootmi", level = level)
  boot_fit$pooled
}

plssemMIBOOT_PS <- function(model, data, ..., m = 5, miArgs = list(),
                            csemArgs = list(), miPackage = "mice",
                            verbose = FALSE, seed = NULL, level = 0.95) {
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }
  imputedData <- NULL

  if (missing(data)) {
    stop("a dataset is needed to run the plssemMIBOOT_PS() function.")
  }

  if (any(is.na(data))) {
    if (miPackage[1] == "Amelia") {
      requireNamespace("Amelia")
      if (!"package:Amelia" %in% search()) 
        attachNamespace("Amelia")
      imputeCall <- c(list(Amelia::amelia, x = data, m = m, p2s = 0), miArgs)
      imputedData <- unclass(eval(as.call(imputeCall))$imputations)
    }
    else if (miPackage[1] == "mice") {
      requireNamespace("mice")
      if (!"package:mice" %in% search()) 
        attachNamespace("mice")
      imputeCall <- c(list(mice::mice, data = data, m = m, 
                           diagnostics = FALSE, printFlag = FALSE), miArgs)
      miceOut <- eval(as.call(imputeCall))
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
  fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
  fit$PathList <- lapply(fit$FitList,
    function(x) {
      path <- x$Estimates$Path_estimates
      path <- matrix2vec(path, names = TRUE)
      rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
      path[path != 0, ]
    })
  fit$LoadingList <- lapply(fit$FitList,
    function(x) {
      load <- colSums(x$Estimates$Loading_estimates)
      names(load) <- gsub("y", "lambda", names(load))
      load
    })
  fit$PooledSample <- pool_samples(fit)
  fit$DataList <- imputedData
  fit$model <- model
  fit$call <- CALL
  fit$csemListCall <- csemListCall
  if (any(is.na(data))) fit$imputeCall <- imputeCall
  fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
  if (!all(unlist(fit$convList))) 
    warning("the model did not converge for all imputed data sets.")
  fit$dots <- dots
  fit$nobs <- nrow(data)
  fit$pooled <- poolMI(fit, boot_mi = "miboot_pooled", level = level)
  fit$pooled
}

plssemBOOTMI_PS <- function(model, data, ..., m = 5, miArgs = list(),
                            csemArgs = list(), miPackage = "mice",
                            bootArgs = list(), verbose = FALSE,
                            seed = NULL, level = 0.95) {
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }

  if (missing(data)) {
    stop("a dataset is needed to run the plssemBOOTMI_PS() function.")
  }

  bootstrap_mi_ps <- function(data, indices, mipkg, miargs, miruns, csemmodel, csemargs, verb) {
    boot_sample <- data[indices, ]

    imputedData <- NULL
    if (any(is.na(data))) {
      if (mipkg[1] == "Amelia") {
        requireNamespace("Amelia")
        if (!"package:Amelia" %in% search()) 
          attachNamespace("Amelia")
        imputeCall <- c(list(Amelia::amelia, x = boot_sample, m = miruns, p2s = 0), miargs)
        imputedData <- unclass(eval(as.call(imputeCall))$imputations)
      }
      else if (mipkg[1] == "mice") {
        requireNamespace("mice")
        if (!"package:mice" %in% search()) 
          attachNamespace("mice")
        imputeCall <- c(list(mice::mice, data = boot_sample, m = miruns, 
                             diagnostics = FALSE, printFlag = FALSE), miargs)
        miceOut <- eval(as.call(imputeCall))
        imputedData <- lapply(as.list(1:m), function(i) mice::complete(data = miceOut,
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
    fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
    fit$ParamList <- lapply(fit$FitList,
      function(x) {
        path <- x$Estimates$Path_estimates
        path <- matrix2vec(path, names = TRUE)
        rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
        path <- path[path != 0, ]
        load <- colSums(x$Estimates$Loading_estimates)
        names(load) <- gsub("y", "lambda", names(load))
        c(path, load)
      })
    fit$PooledMI <- do.call(cbind, fit$ParamList)
    fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
    if (!all(unlist(fit$convList))) 
      warning("the model did not converge for all imputed data sets.")
    fit$DataList <- imputedData
    fit$PooledMI
  }

  if (!is.null(bootArgs) && !missing(bootArgs)) {
    if (bootArgs$parallel != "no" && bootArgs$ncpus > 1)
      verbose <- FALSE
  }
  boot_fit <- list()

  ###
  # new implementation
  # fictitious step only needed to generate the bootstrap indices
  bootListCall <- list(boot::boot, data = data, statistic = function(d, i) 0, R = csemArgs$.R)
  bootListCall <- c(bootListCall, bootArgs)
  bootRes <- eval(as.call(bootListCall))
  indices_list <- rbind(1:nrow(data), boot::boot.array(bootRes, indices = TRUE))
  if (.Platform$OS.type == "unix") {
    bootRes <- parallel::mclapply(seq_len(nrow(indices_list)),  # this can be used only with macOS computers
                                  function(r) {
                                    idx <- indices_list[r, ]
                                    bootstrap_mi_ps(data = data, indices = idx, mipkg = miPackage,
                                                    miargs = miArgs, miruns = m, csemmodel = model,
                                                    csemargs = csemArgs, verb = verbose)
                                  },
                                  mc.cores = argsBOOT$ncpus)
  }
  else {
    bootRes <- lapply(seq_len(nrow(indices_list)),  # under Windows we should use parLapply() but it is too cumbersome
                      function(r) {
                        idx <- indices_list[r, ]
                        bootstrap_mi_ps(data = data, indices = idx, mipkg = miPackage,
                                        miargs = miArgs, miruns = m, csemmodel = model,
                                        csemargs = csemArgs, verb = verbose)
    })
  }
  boot_fit$BootOrig <- bootRes[[1]]   # these are the cSEM results on the original dataset after imputation
  boot_fit$BootMatrix <- t(sapply(bootRes[2:length(bootRes)],
                                  function(x) as.vector(x),
                                  simplify = TRUE))
  ###

  ###
  # old implementation (to delete at a certain point)
  # bootListCall <- list(boot::boot, data = data, statistic = bootstrap_mi_ps,
  #   R = csemArgs$.R, mipkg = miPackage, miargs = miArgs, miruns = m,
  #   csemmodel = model, csemargs = csemArgs, verb = verbose)
  # bootListCall <- c(bootListCall, bootArgs)
  # bootRes <- eval(as.call(bootListCall))
  # boot_fit$BootOrig <- bootRes$t0   # these are the cSEM results on the original dataset after imputation
  # boot_fit$BootMatrix <- bootRes$t
  ###

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

plssemWGT_BOOTMI <- function(model, data, ..., m = 5, miArgs = list(),
                             csemArgs = list(), miPackage = "mice",
                             bootArgs = list(), verbose = FALSE, seed = NULL,
                             level = 0.95, wgtType = "rows") {
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) && !missing(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }

  if (missing(data)) {
    stop("a dataset is needed to run the plssemWGT_BOOTMI() function.")
  }

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
        if (!"package:Amelia" %in% search()) 
          attachNamespace("Amelia")
        imputeCall <- c(list(Amelia::amelia, x = boot_sample, m = miruns, p2s = 0), miargs)
        imputedData <- unclass(eval(as.call(imputeCall))$imputations)
      }
      else if (mipkg[1] == "mice") {
        requireNamespace("mice")
        if (!"package:mice" %in% search()) 
          attachNamespace("mice")
        imputeCall <- c(list(mice::mice, data = boot_sample, m = miruns, 
                             diagnostics = FALSE, printFlag = FALSE), miargs)
        miceOut <- eval(as.call(imputeCall))
        imputedData <- lapply(as.list(1:m), function(i) mice::complete(data = miceOut,
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
    fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
    fit$PathList <- lapply(fit$FitList,
      function(x) {
        path <- x$Estimates$Path_estimates
        path <- matrix2vec(path, names = TRUE)
        rownames(path) <- gsub("eta", "gamma", gsub("_eta", "", rownames(path), fixed = TRUE), fixed = TRUE)
        path[path != 0, ]
      })
    fit$LoadingList <- lapply(fit$FitList,
      function(x) {
        load <- colSums(x$Estimates$Loading_estimates)
        names(load) <- gsub("y", "lambda", names(load))
        load
      })
    fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
    if (!all(unlist(fit$convList))) 
      warning("the model did not converge for all imputed data sets.")
    fit$DataList <- imputedData
    c(wgt, rubin_est(fit$PathList), rubin_est(fit$LoadingList))
  }
  boot_fit <- list()

  ###
  # new implementation
  # fictitious step only needed to generate the bootstrap indices
  bootListCall <- list(boot::boot, data = data, statistic = function(d, i) 0, R = csemArgs$.R)
  bootListCall <- c(bootListCall, bootArgs)
  bootRes <- eval(as.call(bootListCall))
  indices_list <- rbind(1:nrow(data), boot::boot.array(bootRes, indices = TRUE))
  if (.Platform$OS.type == "unix") {  # parallel computing gives an error sometimes, so we revert to non-parallel
    bootRes <- parallel::mclapply(seq_len(nrow(indices_list)),
                                  function(r) {
                                    idx <- indices_list[r, ]
                                    w_bootstrap_mi(data = data, indices = idx, mipkg = miPackage,
                                                   miargs = miArgs, miruns = m, csemmodel = model,
                                                   csemargs = csemArgs, verb = verbose, wtype = wgtType)
                                  },
                                  mc.cores = argsBOOT$ncpus)
  }
  else {
    bootRes <- lapply(seq_len(nrow(indices_list)),
                      function(r) {
                        idx <- indices_list[r, ]
                        w_bootstrap_mi(data = data, indices = idx, mipkg = miPackage,
                                       miargs = miArgs, miruns = m, csemmodel = model,
                                       csemargs = csemArgs, verb = verbose, wtype = wgtType)
                      })
  }
  boot_fit$BootOrig <- (bootRes[[1]])[-1]   # these are the cSEM results on the original dataset after imputation
  boot_fit$BootMatrix <- (do.call(rbind, bootRes[2:length(bootRes)]))[, -1, drop = FALSE]
  wgts <- (do.call(rbind, bootRes[2:length(bootRes)]))[, 1]
  ###

  ###
  # old implementation (to delete at a certain point)
  # bootListCall <- list(boot::boot, data = data, statistic = w_bootstrap_mi,
  #   R = csemArgs$.R, mipkg = miPackage, miargs = miArgs, miruns = m,
  #   csemmodel = model, csemargs = csemArgs, verb = FALSE, wtype = wgtType)
  # bootListCall <- c(bootListCall, bootArgs)
  # bootRes <- eval(as.call(bootListCall))
  # boot_fit$BootOrig <- bootRes$t0   # these are the cSEM results on the original dataset after imputation
  # boot_fit$BootMatrix <- bootRes$t
  ###

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
