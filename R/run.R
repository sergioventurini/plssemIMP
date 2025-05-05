run_sims <- function(
  runs = 10,
  argsCD = list(n = 1e3, method = "model", pkg = "cSEM.DGP"),
  argsMM = list(prop = .5, mech = "MCAR", method = "ampute"),
  argsMI = list(m = 5, methods = c("pmm", "norm"), pkg = "mice"),
  argscSEM = list(),
  argsBOOT = list(),
  boot_mi = "miboot", # accepted values are 'miboot', 'bootmi', 'miboot_pooled', 'bootmi_pooled' and 'weighted_bootmi'
  verbose = FALSE,
  seed = NULL,
  level = 0.95,
  meanimp = TRUE,
  knnimp = TRUE, argsKNN = list(k = 5, method = "euclidean"),
  listwise = TRUE, fulloriginal = TRUE,
  datalist = NULL) {

  CALL <- match.call()

  nCD <- argsCD$n
  methodCD <- argsCD$method
  pkgCD <- argsCD$pkg
  argsCD <- argsCD[setdiff(names(argsCD), c("method", "pkg", "n"))]

  propMM <- argsMM$prop
  mechMM <- argsMM$mech
  methodMM <- argsMM$method
  argsMM <- argsMM[setdiff(names(argsMM), c("prop", "mech", "method"))]

  methodsMI <- argsMI$methods
  mMI <- argsMI$m
  pkgMI <- argsMI$pkg
  modelMI <- argsMI$model
  argsMI <- argsMI[setdiff(names(argsMI), c("methods", "m", "pkg", "model"))]

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  start_seed <- .Random.seed

  if (!is.null(datalist)) {
    if (length(datalist) < runs)
      stop("the number of datasets provided is smaller than the number of runs.")
  }

  meth_robust <- c("gamlssBI", "gamlssGA", "gamlssJSU", "gamlssNO",
                   "gamlssPO", "gamlssTF", "gamlssZIBI", "gamlssZIP")
  if (any(methodsMI %in% meth_robust) &&
      !"package:ImputeRobust" %in% search()) {
    suppressMessages(require("ImputeRobust", quietly = TRUE))
  }

  res_all <- res <- run_seeds <- list()
  for (run in 1:runs) {
    run_seeds[[run]] <- .Random.seed
    if (verbose)
      cat(paste0("Simulation run ", run, " of ", runs, "\n"))

    ## STEP 1: generate complete data
    if (is.null(datalist)) {
      dat <- create_data(n = nCD, method = methodCD, pkg = pkgCD, args = argsCD)
    }
    else {
      dat <- datalist[[run]]
      if (any(is.na(dat)))
        stop("the datasets provided must be complete.")
    }
    dat_orig <- dat
    
    ## STEP 2: generate missing data
    dat_miss <- make_missing(dat, prop = propMM, mech = mechMM, method = methodMM,
      missArgs = argsMM)
    dat <- as.data.frame(dat_miss$amputed)
    
    ## STEP 3: perform multiple imputation & bootstrap
    res <- list()
    for (methMI in methodsMI) {
      miArgs <- c(method = methMI, argsMI)
      if (verbose)
        cat(paste0("  - multiple imputation package/method: ", pkgMI, "/", methMI, "\n"))
      if (boot_mi == "miboot") {
        if (is.null(argscSEM$.resample_method)) {
          argscSEM <- c(argscSEM, .resample_method = "bootstrap")
          warning("the .resample_method option has been set to 'bootstrap'.")
        }
        else if (argscSEM$.resample_method != "bootstrap") {
          argscSEM$.resample_method <- "bootstrap"
          warning("the .resample_method option has been set to 'bootstrap'.")
        }
        res[[methMI]] <- plssemMIBOOT(model = modelMI, data = dat, m = mMI,
                                      miArgs = miArgs, miPackage = pkgMI, csemArgs = argscSEM,
                                      verbose = verbose, seed = NULL, level = level)
      }
      else if (boot_mi == "bootmi") {
        if (!is.null(argscSEM$.resample_method) && argscSEM$.resample_method != "none") {
          argscSEM$.resample_method <- "none"
          warning("the .resample_method option has been set to 'none'.")
        }
        res[[methMI]] <- plssemBOOTMI(model = modelMI, data = dat, m = mMI,
                                      miArgs = miArgs, miPackage = pkgMI,
                                      csemArgs = argscSEM, bootArgs = argsBOOT,
                                      verbose = verbose, seed = NULL, level = level)
      }
      else if (boot_mi == "miboot_pooled") {
        if (is.null(argscSEM$.resample_method)) {
          argscSEM <- c(argscSEM, .resample_method = "bootstrap")
          warning("the .resample_method option has been set to 'bootstrap'.")
        }
        else if (argscSEM$.resample_method != "bootstrap") {
          argscSEM$.resample_method <- "bootstrap"
          warning("the .resample_method option has been set to 'bootstrap'.")
        }
        res[[methMI]] <- plssemMIBOOT_PS(model = modelMI, data = dat, m = mMI,
                                         miArgs = miArgs, miPackage = pkgMI,
                                         csemArgs = argscSEM, verbose = verbose,
                                         seed = NULL, level = level)
      }
      else if (boot_mi == "bootmi_pooled") {
        if (!is.null(argscSEM$.resample_method) && argscSEM$.resample_method != "none") {
          argscSEM$.resample_method <- "none"
          warning("the .resample_method option has been set to 'none'.")
        }
        res[[methMI]] <- plssemBOOTMI_PS(model = modelMI, data = dat, m = mMI,
                                         miArgs = miArgs, miPackage = pkgMI,
                                         csemArgs = argscSEM, bootArgs = argsBOOT,
                                         verbose = verbose, seed = NULL, level = level)
      }
      else if (boot_mi == "weighted_bootmi") {
        if (!is.null(argscSEM$.resample_method) && argscSEM$.resample_method != "none") {
          argscSEM$.resample_method <- "none"
          warning("the .resample_method option has been set to 'none'.")
        }
        if (argsBOOT$parallel != "no") {
          argsBOOT$parallel <- "no"
          warning("the weighted version of 'bootmi' can't use parallel computation and it will take longer.")
        }
        res[[methMI]] <- plssemWGT_BOOTMI(model = modelMI, data = dat, m = mMI,
                                          miArgs = miArgs, miPackage = pkgMI,
                                          csemArgs = argscSEM, bootArgs = argsBOOT,
                                          verbose = verbose, seed = NULL, level = level)
      }
      else {
        stop("the selected bootstrap/multiple imputation approach is not available.")
      }
    }

    ## STEP 4: perform single imputation
    if (boot_mi == "bootmi" | boot_mi == "bootmi_pooled" | boot_mi == "weighted_bootmi")
      argscSEM$.resample_method <- "bootstrap"

    # mean imputation
    if (meanimp) {
      if (verbose)
        cat("  - single imputation method: mean\n")
      res[["mean"]] <- meanimp(model = modelMI, data = dat,
                               csemArgs = argscSEM, verbose = verbose,
                               level = level)
    }

    # k-nearest neighbor imputation
    if (knnimp) {
      if (is.null(argsKNN$method))
        argsKNN$method <- "euclidean"
      if (is.null(argsKNN$k)) {
        stop("specify at least on k value.")
      }
      else if (any(!is.numeric(argsKNN$k))) {
        stop("k must contain numeric values.")
      }
      else if (any(floor(argsKNN$k) < 1)) {
        stop("k must contain integer values greater than 0.")
      }
      argsKNN$k <- floor(argsKNN$k)
      for (j in 1:length(argsKNN$k)) {
        if (verbose)
          cat(paste0("  - single imputation method: k-nearest neighbors (k = ", argsKNN$k[j], ")\n"))
        res[[paste0("knn", argsKNN$k[j])]] <- knnimp(model = modelMI, data = dat,
                                                     csemArgs = argscSEM,
                                                     k = argsKNN$k[j],
                                                     method = argsKNN$method[1],
                                                     verbose = verbose,
                                                     level = level)
      }
    }

    # perform complete-case (i.e. listwise deletion) analysis
    if (listwise) {
      if (verbose)
        cat("  - complete-case (listwise deletion) analysis\n")
      res[["listwise"]] <- fulldata(model = modelMI, data = na.omit(dat),
                                    csemArgs = argscSEM, verbose = verbose,
                                    level = level)
    }

    # perform analysis on original full data (i.e. before generating NAs)
    if (fulloriginal) {
      if (verbose)
        cat("  - analysis on the original full data (i.e. before generating NAs)\n")
      res[["fulloriginal"]] <- fulldata(model = modelMI, data = dat_orig,
                                        csemArgs = argscSEM, verbose = verbose,
                                        level = level)
    }

    methods_names <- methodsMI
    if (meanimp) methods_names <- c(methods_names, "mean")
    if (knnimp) methods_names <- c(methods_names, paste0("knn_", argsKNN$k))
    if (fulloriginal) methods_names <- c(methods_names, "fulloriginal")
    if (listwise) methods_names <- c(methods_names, "listwise")
    names(res) <- methods_names
    res$dat_orig <- dat_orig
    res$dat_miss <- dat
    res_all[[run]] <- res
  }

  names(res_all) <- paste0("run_", 1:runs)
  res_all$start_seed <- start_seed
  res_all$run_seeds <- run_seeds
  res_all$call <- CALL

  res_all
}
