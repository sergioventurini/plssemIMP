run_sims <- function(
  runs = 1e2,
  argsCD = list(n = 1e3, method = "model", pkg = "cSEM.DGP"),
  argsMM = list(prop = .5, mech = "MCAR", method = "ampute"),
  argsMI = list(m = 5, methods = c("pmm", "norm"), pkg = "mice"),
  argscSEM = list(),
  argsBOOT = list(),
  boot_mi = "miboot",  # accepted values are 'miboot', 'bootmi', 'miboot_pooled', 'bootmi_pooled' and 'weighted_bootmi'
  wgtType = "rows",    # accepted values are 'rows' and 'all'
  verbose = FALSE,
  level = 0.95,
  meanimp = TRUE,
  knnimp = TRUE, argsKNN = list(k = 5, method = "euclidean"),
  listwise = TRUE, fulloriginal = TRUE,
  datalist = NULL,
  datamisslist = NULL,
  store_data = FALSE,
  runALL = TRUE,
  log_file = NULL) {

  CALL <- match.call()
  arg_names <- names(formals(sys.function()))
  args_list <- mget(arg_names, envir = environment(), inherits = FALSE)

  nCD <- argsCD$n
  methodCD <- argsCD$method
  pkgCD <- argsCD$pkg
  argsCD <- argsCD[setdiff(names(argsCD), c("method", "pkg", "n"))]

  propMM <- argsMM$prop
  mechMM <- argsMM$mech
  methodMM <- argsMM$method
  argsMM <- argsMM[setdiff(names(argsMM), c("prop", "mech", "method"))]

  if (is.null(argsMI$methods)) {
    methodsMI <- "pmm"
  }
  else {
    methodsMI <- argsMI$methods
  }
  modelMI <- argsMI$model
  argsMI <- argsMI[setdiff(names(argsMI), c("methods", "model"))]

  start_seed <- .Random.seed

  if (!is.null(datalist)) {
    if (length(datalist) < runs)
      stop("the number of data sets provided is smaller than the number of runs.")
  }
  if (!is.null(datamisslist)) {
    if (length(datamisslist) < runs)
      stop("the number of data sets provided is smaller than the number of runs.")
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
    if (verbose) {
      cat(paste0("Simulation run ", run, " of ", runs, "\n"))
    }
    log_msg(paste0("Simulation run ", run, " of ", runs), log_file)

    #################################
    ## STEP 1: generate complete data
    #################################
    if (is.null(datalist)) {
      dat <- create_data(n = nCD, method = methodCD, pkg = pkgCD, args = argsCD)
    }
    else {
      dat <- datalist[[run]]
      if (any(is.na(dat)))
        stop("the data sets provided must be complete.")
    }
    dat_orig <- dat
    
    ################################
    ## STEP 2: generate missing data
    ################################
    if (is.null(datamisslist)) {
      dat_miss <- make_missing(dat, prop = propMM, mech = mechMM, method = methodMM,
                               missArgs = argsMM)
      dat <- as.data.frame(dat_miss$amputed)
    }
    else {
      dat <- datamisslist[[run]]
    }
    
    ##################################################
    ## STEP 3: perform multiple imputation & bootstrap
    ##################################################
    if (runALL) {
      for (methMI in methodsMI) {
        miArgs <- c(method = methMI, argsMI)
        res_mi <- plssemIMP(data = dat,
                            model = modelMI,
                            argsMI = miArgs,
                            argscSEM = argscSEM,
                            argsBOOT = argsBOOT,
                            boot_mi = boot_mi,
                            wgtType = wgtType,
                            verbose = verbose,
                            seed = NULL,
                            level = level,
                            log_file = log_file)
        res[[methMI]] <- res_mi$res
      }

      ####################################
      ## STEP 4: perform single imputation
      ####################################
      if (boot_mi %in% c("bootmi", "bootmi_pooled", "weighted_bootmi"))
        argscSEM$.resample_method <- "bootstrap"

      # mean imputation
      if (meanimp) {
        if (verbose) {
          cat("  - single imputation method: mean\n")
        }
        log_msg("  - single imputation method: mean", log_file)
        res[["mean"]] <- meanimp(model = modelMI, data = dat,
                                 csemArgs = argscSEM, verbose = verbose,
                                 level = level)
      }

      # k-nearest-neighbour imputation
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
          if (verbose) {
            cat(paste0("  - single imputation method: k-nearest neighbors (k = ", argsKNN$k[j], ")\n"))
          }
          log_msg(paste0("  - single imputation method: k-nearest neighbors (k = ", argsKNN$k[j], ")"), log_file)
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
        if (verbose) {
          cat("  - complete-case (listwise deletion) analysis\n")
        }
        log_msg("  - complete-case (listwise deletion) analysis", log_file)
        res[["listwise"]] <- fulldata(model = modelMI, data = na.omit(dat),
                                      csemArgs = argscSEM, verbose = verbose,
                                      level = level)
      }

      # perform analysis on original full data (i.e. before generating NAs)
      if (fulloriginal) {
        if (verbose){
          cat("  - analysis on the original full data (i.e. before generating NAs)\n")
        }
        log_msg("  - analysis on the original full data (i.e. before generating NAs)", log_file)
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
    }

    if (store_data) {
      res$dat_orig <- dat_orig
      res$dat_miss <- dat
    }

    res_all[[run]] <- res
  }

  names(res_all) <- paste0("run_", 1:runs)
  res_all$start_seed <- start_seed
  res_all$run_seeds <- run_seeds
  res_all$call <- CALL
  res_all$args <- args_list

  res_all
}
