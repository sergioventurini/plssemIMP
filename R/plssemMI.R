plssemMIBOOT <- function(model, data, ..., m = 5, miArgs = list(),
                         csemArgs = list(), miPackage = "mice",
                         verbose = FALSE, seed = NULL, level = 0.95) {
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) & !missing(seed)) {
    set.seed(seed = seed)
  }
  imputedData <- NULL

  if (missing(data)) {
    stop("a dataset is needed to run the plssemMIBOOT() function.")
  }

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
    imputedData <- list()
    for (i in 1:m) {
      imputedData[[i]] <- mice::complete(data = miceOut, 
                                         action = i, include = FALSE)
    }
  }
  else stop("currently plssemMIBOOT() only supports imputation by Amelia or mice.")

  csemListCall <- list(cSEM::csem, .model = model, .data = imputedData)
  csemListCall <- c(csemListCall, csemArgs)
  fit <- list()
  fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
  fit$PathList <- lapply(fit$FitList,
    function(x) {
      path <- x$Estimates$Path_estimates
      path <- path[, colSums(path) != 0]
      colSums(path)
    })
  fit$LoadingList <- lapply(fit$FitList,
    function(x) {
      colSums(x$Estimates$Loading_estimates)
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
  fit$imputeCall <- imputeCall
  fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
  if (!all(unlist(fit$convList))) 
    warning("the model did not converge for any imputed data sets.")
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
  if (!is.null(seed) & !missing(seed)) {
    set.seed(seed = seed)
  }

  if (missing(data)) {
    stop("a dataset is needed to run the plssemBOOTMI() function.")
  }

  boot_i <- 0
  bootstrap_mi <- function(data, indices, mipkg, miargs, miruns, csemmodel, csemargs, verb) {
    boot_i <<- boot_i + 1
    if (verb)
      cat(paste0("  - bootstrap sample ", boot_i, "\n"))
    boot_sample <- data[indices, ]

    imputedData <- NULL
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
      imputedData <- list()
      for (i in 1:m) {
        imputedData[[i]] <- mice::complete(data = miceOut, 
                                           action = i, include = FALSE)
      }
    }
    else stop("currently plssemBOOTMI() only supports imputation by Amelia or mice.")

    csemListCall <- list(cSEM::csem, .model = csemmodel, .data = imputedData)
    csemListCall <- c(csemListCall, csemargs)
    fit <- list()
    fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
    fit$PathList <- lapply(fit$FitList,
      function(x) {
        path <- x$Estimates$Path_estimates
        path <- path[, colSums(path) != 0]
        colSums(path)
      })
    fit$LoadingList <- lapply(fit$FitList,
      function(x) {
        colSums(x$Estimates$Loading_estimates)
      })
    fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
    if (!all(unlist(fit$convList))) 
      warning("the model did not converge for any imputed datasets.")
    fit$DataList <- imputedData
    c(rubin_est(fit$PathList), rubin_est(fit$LoadingList))
  }

  if (!is.null(bootArgs) & !missing(bootArgs)) {
    if (bootArgs$parallel != "no" & bootArgs$ncpus > 1)
      verbose <- FALSE
  }
  boot_fit <- list()
  bootListCall <- list(boot::boot, data = data, statistic = bootstrap_mi,
    R = csemArgs$.R, mipkg = miPackage, miargs = miArgs, miruns = m,
    csemmodel = model, csemargs = csemArgs, verb = verbose)
  bootListCall <- c(bootListCall, bootArgs)
  bootRes <- eval(as.call(bootListCall))
  boot_fit$BootOrig <- bootRes$t0   # these are the cSEM results on the original dataset after imputation
  boot_fit$BootMatrix <- bootRes$t
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
  if (!is.null(seed) & !missing(seed)) {
    set.seed(seed = seed)
  }
  imputedData <- NULL

  if (missing(data)) {
    stop("a dataset is needed to run the plssemMIBOOT() function.")
  }

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
    imputedData <- list()
    for (i in 1:m) {
      imputedData[[i]] <- mice::complete(data = miceOut, 
                                         action = i, include = FALSE)
    }
  }
  else stop("currently plssemMIBOOT_PS() only supports imputation by Amelia or mice.")

  csemListCall <- list(cSEM::csem, .model = model, .data = imputedData)
  csemListCall <- c(csemListCall, csemArgs)
  fit <- list()
  fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
  fit$PooledSample <- pool_samples(fit)
  fit$DataList <- imputedData
  fit$model <- model
  fit$call <- CALL
  fit$csemListCall <- csemListCall
  fit$imputeCall <- imputeCall
  fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
  if (!all(unlist(fit$convList))) 
    warning("the model did not converge for any imputed data sets.")
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
  if (!is.null(seed) & !missing(seed)) {
    set.seed(seed = seed)
  }

  if (missing(data)) {
    stop("a dataset is needed to run the plssemBOOTMI() function.")
  }

  boot_i <- 0
  bootstrap_mi_ps <- function(data, indices, mipkg, miargs, miruns, csemmodel, csemargs, verb) {
    boot_i <<- boot_i + 1
    if (verb)
      cat(paste0("  - bootstrap sample ", boot_i, "\n"))
    boot_sample <- data[indices, ]

    imputedData <- NULL
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
      imputedData <- list()
      for (i in 1:m) {
        imputedData[[i]] <- mice::complete(data = miceOut, 
                                           action = i, include = FALSE)
      }
    }
    else stop("currently plssemBOOTMI_PS() only supports imputation by Amelia or mice.")

    csemListCall <- list(cSEM::csem, .model = csemmodel, .data = imputedData)
    csemListCall <- c(csemListCall, csemargs)
    fit <- list()
    fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
    fit$ParamList <- lapply(fit$FitList,
      function(x) {
        path <- x$Estimates$Path_estimates
        path <- path[, colSums(path) != 0]
        path <- colSums(path)
        load <- colSums(x$Estimates$Loading_estimates)
        c(path, load)
      })
    fit$PooledMI <- do.call(cbind, fit$ParamList)
    fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
    if (!all(unlist(fit$convList))) 
      warning("the model did not converge for any imputed datasets.")
    fit$DataList <- imputedData
    fit$PooledMI
  }

  if (!is.null(bootArgs) & !missing(bootArgs)) {
    if (bootArgs$parallel != "no" & bootArgs$ncpus > 1)
      verbose <- FALSE
  }
  boot_fit <- list()
  bootListCall <- list(boot::boot, data = data, statistic = bootstrap_mi_ps,
    R = csemArgs$.R, mipkg = miPackage, miargs = miArgs, miruns = m,
    csemmodel = model, csemargs = csemArgs, verb = verbose)
  bootListCall <- c(bootListCall, bootArgs)
  bootRes <- eval(as.call(bootListCall))
  boot_fit$BootOrig <- bootRes$t0   # these are the cSEM results on the original dataset after imputation
  boot_fit$BootMatrix <- bootRes$t
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
                             level = 0.95) {
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(seed) & !missing(seed)) {
    set.seed(seed = seed)
  }

  if (missing(data)) {
    stop("a dataset is needed to run the plssemBOOTMI() function.")
  }

  boot_i <- 0
  wgts <- numeric(csemArgs$.R)
  bootstrap_mi <- function(data, indices, mipkg, miargs, miruns, csemmodel, csemargs, verb) {
    boot_i <<- boot_i + 1
    if (verb)
      cat(paste0("  - bootstrap sample ", boot_i, "\n"))
    boot_sample <- data[indices, ]
    wgts[boot_i] <<- sum(is.na(boot_sample))

    imputedData <- NULL
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
      imputedData <- list()
      for (i in 1:m) {
        imputedData[[i]] <- mice::complete(data = miceOut, 
                                           action = i, include = FALSE)
      }
    }
    else stop("currently plssemBOOTMI() only supports imputation by Amelia or mice.")

    csemListCall <- list(cSEM::csem, .model = csemmodel, .data = imputedData)
    csemListCall <- c(csemListCall, csemargs)
    fit <- list()
    fit$FitList <- suppressWarnings(eval(as.call(csemListCall)))
    fit$PathList <- lapply(fit$FitList,
      function(x) {
        path <- x$Estimates$Path_estimates
        path <- path[, colSums(path) != 0]
        colSums(path)
      })
    fit$LoadingList <- lapply(fit$FitList,
      function(x) {
        colSums(x$Estimates$Loading_estimates)
      })
    fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
    if (!all(unlist(fit$convList))) 
      warning("the model did not converge for any imputed datasets.")
    fit$DataList <- imputedData
    c(rubin_est(fit$PathList), rubin_est(fit$LoadingList))
  }

  boot_fit <- list()
  bootListCall <- list(boot::boot, data = data, statistic = bootstrap_mi,
    R = csemArgs$.R, mipkg = miPackage, miargs = miArgs, miruns = m,
    csemmodel = model, csemargs = csemArgs, verb = FALSE)
  bootListCall <- c(bootListCall, bootArgs)
  bootRes <- eval(as.call(bootListCall))
  # wgts <- prod(dim(data))/wgts[-1]
  wgts <- nrow(data)/wgts[-1]
  boot_fit$BootOrig <- bootRes$t0   # these are the cSEM results on the original dataset after imputation
  boot_fit$BootMatrix <- bootRes$t
  boot_fit$model <- model
  boot_fit$call <- CALL
  boot_fit$dots <- dots
  boot_fit$nobs <- nrow(data)
  boot_fit$pooled <- poolMI(boot_fit, boot_mi = "weighted_bootmi", level = level, weights = wgts)
  boot_fit$pooled
}
