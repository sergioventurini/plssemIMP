run_sims <- function(
  runs = 1e2,
  argsCD = list(n = 1e3, method = "model", pkg = "cSEM.DGP"),
  argsMM = list(prop = .5, mech = "MCAR", method = "ampute"),
  argsMI = list(m = 5, methods = c("pmm", "norm"), pkg = "mice"),
  argscSEM = list(),
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
  log_file = NULL,
  runs_parallel = 1L,    # <<< NEW: number of concurrent simulation runs.
                         # Must satisfy runs_parallel * boot_mc_cores <= total cores.
                         # When > 1: verbose is forced FALSE and per-run logging is
                         # disabled to avoid garbled output from concurrent processes.
  boot_mc_cores = NULL)  # <<< NEW: cores per run for the inner bootstrap mclapply.
                         # NULL = use all available cores minus one (standalone behaviour).
                         # Typically set to floor(total_cores / runs_parallel).
{
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

  # CHANGE: resolve effective run-level parallelism
  effective_runs_parallel <- max(1L, as.integer(runs_parallel))  # <<< NEW

  # CHANGE: when runs are parallelised, cSEM's own .eval_plan must be "sequential"
  # to prevent a third level of forked processes competing for the same cores.
  # We override it here rather than requiring the caller to remember this rule.
  if (effective_runs_parallel > 1L &&                               # <<< NEW
      !is.null(argscSEM$.eval_plan) &&                             # <<< NEW
      argscSEM$.eval_plan != "sequential") {                       # <<< NEW
    argscSEM$.eval_plan <- "sequential"                            # <<< NEW
    message("run_sims: runs_parallel > 1, forcing argscSEM$.eval_plan = 'sequential'.")
  }                                                                 # <<< NEW

  # CHANGE: define the body of a single simulation run as a self-contained function.
  # This is the function that mclapply will call for each run index.
  # Returning a named list allows us to recover both the results and the RNG seed
  # used by each worker (for reproducibility diagnostics).
  run_one <- function(run) {                                        # <<< NEW FUNCTION

    # Capture the RNG state at the start of this worker. With mc.set.seed = TRUE,
    # each forked worker inherits a distinct L'Ecuyer-CMRG stream, so these seeds
    # will differ across runs even when runs_parallel > 1.
    run_seed_captured <- .Random.seed

    # Suppress verbose output and per-run logging when running in parallel:
    # multiple processes writing to the same stream simultaneously produces garbled output.
    run_verbose  <- if (effective_runs_parallel > 1L) FALSE else verbose
    run_log_file <- if (effective_runs_parallel > 1L) NULL  else log_file

    if (run_verbose) {
      cat(paste0("Simulation run ", run, " of ", runs, "\n"))
    }
    log_msg(paste0("Simulation run ", run, " of ", runs), run_log_file)

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

    res <- list()

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
                            boot_mi = boot_mi,
                            wgtType = wgtType,
                            verbose = run_verbose,
                            seed = NULL,
                            level = level,
                            log_file = run_log_file,
                            boot_mc_cores = boot_mc_cores)    # <<< NEW: forwarded
        res[[methMI]] <- res_mi$res
      }

      ####################################
      ## STEP 4: perform single imputation
      ####################################
      if (boot_mi %in% c("bootmi", "bootmi_pooled", "weighted_bootmi"))
        argscSEM$.resample_method <- "bootstrap"

      # mean imputation
      if (meanimp) {
        if (run_verbose) {
          cat("  - single imputation method: mean\n")
        }
        log_msg("  - single imputation method: mean", run_log_file)
        res[["mean"]] <- meanimp(model = modelMI, data = dat,
                                 csemArgs = argscSEM, verbose = run_verbose,
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
          if (run_verbose) {
            cat(paste0("  - single imputation method: k-nearest neighbors (k = ", argsKNN$k[j], ")\n"))
          }
          log_msg(paste0("  - single imputation method: k-nearest neighbors (k = ", argsKNN$k[j], ")"), run_log_file)
          res[[paste0("knn", argsKNN$k[j])]] <- knnimp(model = modelMI, data = dat,
                                                       csemArgs = argscSEM,
                                                       k = argsKNN$k[j],
                                                       method = argsKNN$method[1],
                                                       verbose = run_verbose,
                                                       level = level)
        }
      }

      # perform complete-case (i.e. listwise deletion) analysis
      if (listwise) {
        if (run_verbose) {
          cat("  - complete-case (listwise deletion) analysis\n")
        }
        log_msg("  - complete-case (listwise deletion) analysis", run_log_file)
        res[["listwise"]] <- fulldata(model = modelMI, data = na.omit(dat),
                                      csemArgs = argscSEM, verbose = run_verbose,
                                      level = level)
      }

      # perform analysis on original full data (i.e. before generating NAs)
      if (fulloriginal) {
        if (run_verbose){
          cat("  - analysis on the original full data (i.e. before generating NAs)\n")
        }
        log_msg("  - analysis on the original full data (i.e. before generating NAs)", run_log_file)
        res[["fulloriginal"]] <- fulldata(model = modelMI, data = dat_orig,
                                          csemArgs = argscSEM, verbose = run_verbose,
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

    # Return both the results and the captured seed for reproducibility diagnostics
    list(res = res, seed = run_seed_captured)                      # <<< NEW
  }  # end run_one                                                  # <<< NEW

  # CHANGE: dispatch runs either sequentially (runs_parallel = 1) or in parallel.
  # mclapply with mc.set.seed = TRUE propagates distinct L'Ecuyer-CMRG substreams
  # to each forked worker, ensuring reproducible and independent RNG per run.
  # WARNING: mclapply (fork-based) only works on Unix/macOS; on Windows it silently
  # falls back to sequential. Use the runs_parallel = 1 path on Windows.
  if (effective_runs_parallel == 1L) {                             # <<< NEW
    run_results <- lapply(1:runs, run_one)                         # <<< NEW (no fork overhead)
  } else {                                                         # <<< NEW
    run_results <- parallel::mclapply(                             # <<< NEW
      1:runs,                                                      # <<< NEW
      run_one,                                                     # <<< NEW
      mc.cores    = effective_runs_parallel,                       # <<< NEW
      mc.preschedule = FALSE,                                      # <<< NEW
      mc.set.seed = TRUE                                           # <<< NEW
    )                                                              # <<< NEW
  }                                                                # <<< NEW

###
# surface the real child error before the reassembly crash masks it
failed <- sapply(run_results, inherits, "try-error")
if (any(failed)) {
  failed_idx  <- which(failed)
  real_errors <- sapply(run_results[failed], as.character)
  stop(sprintf(
    "mclapply: %d run(s) failed. First failure was run %d.\nReal error:\n%s",
    sum(failed), failed_idx[1], real_errors[1]
  ))
}
###

  # CHANGE: reassemble res_all and run_seeds from the list returned by mclapply
  res_all   <- lapply(run_results, `[[`, "res")                    # <<< NEW
  run_seeds <- lapply(run_results, `[[`, "seed")                   # <<< NEW

  names(res_all) <- paste0("run_", 1:runs)
  res_all$start_seed <- start_seed
  res_all$run_seeds  <- run_seeds
  res_all$call       <- CALL
  res_all$args       <- args_list

  res_all
}
