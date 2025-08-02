plssemIMP <- function(
  data,
  model,
  argsMI = list(m = 5, method = "pmm", pkg = "mice"),
  argscSEM = list(),
  argsBOOT = list(),
  boot_mi = "miboot",  # accepted values are 'miboot', 'bootmi', 'miboot_pooled', 'bootmi_pooled' and 'weighted_bootmi'
  wgtType = "rows",    # accepted values are 'rows' and 'all'
  verbose = FALSE,
  seed = NULL,
  level = 0.95)
{
  CALL <- match.call()
  arg_names <- names(formals(sys.function()))
  args_list <- mget(arg_names, envir = environment(), inherits = FALSE)

  if (is.null(data)) {
    stop("a data set must be provided.")
  }
  if (is.null(model)) {
    stop("a model to estimate must be provided.")
  }
  if (!any(is.na(data)))
    stop("the data set provided is complete and doesn't need imputation.")
  if (wgtType != "rows" && wgtType != "all")
    stop("wgtType can be equal to either 'rows' or 'all'.")
  if (level <= 0 | level >= 1)
    stop("level must be in between 0 and 1.")

  if (is.null(argsMI$method)) {
    methodMI <- "pmm"
  }
  else {
    methodMI <- argsMI$method
  }
  if (is.null(argsMI$m)) {
    mMI <- 5
  }
  else {
    mMI <- argsMI$m
  }
  if (is.null(argsMI$pkg)) {
    pkgMI <- "mice"
  }
  else {
    pkgMI <- argsMI$pkg
  }
  argsMI <- argsMI[setdiff(names(argsMI), c("method", "m", "pkg"))]

  start_seed <- .Random.seed

  meth_robust <- c("gamlssBI", "gamlssGA", "gamlssJSU", "gamlssNO",
                   "gamlssPO", "gamlssTF", "gamlssZIBI", "gamlssZIP")
  if (any(methodMI %in% meth_robust) &&
      !"package:ImputeRobust" %in% search()) {
    suppressMessages(require("ImputeRobust", quietly = TRUE))
  }

  res <- res_all <- list()

  ## perform multiple imputation & bootstrap
  miArgs <- c(method = methodMI, argsMI)
  if (verbose)
    cat(paste0("  - multiple imputation package/method: ", pkgMI, "/", methodMI, "\n"))
  if (boot_mi == "miboot") {
    if (is.null(argscSEM$.resample_method)) {
      argscSEM <- c(argscSEM, .resample_method = "bootstrap")
      if (verbose) warning("the .resample_method option has been set to 'bootstrap'.")
    }
    else if (argscSEM$.resample_method != "bootstrap") {
      argscSEM$.resample_method <- "bootstrap"
      if (verbose) warning("the .resample_method option has been set to 'bootstrap'.")
    }
    res <- plssemMIBOOT(model = model, data = data, m = mMI,
                        miArgs = miArgs, miPackage = pkgMI, csemArgs = argscSEM,
                        verbose = verbose, seed = seed, level = level)
  }
  else if (boot_mi == "bootmi") {
    if (!is.null(argscSEM$.resample_method) && argscSEM$.resample_method != "none") {
      argscSEM$.resample_method <- "none"
      if (verbose) warning("the .resample_method option has been set to 'none'.")
    }
    res <- plssemBOOTMI(model = model, data = data, m = mMI,
                        miArgs = miArgs, miPackage = pkgMI,
                        csemArgs = argscSEM, bootArgs = argsBOOT,
                        verbose = verbose, seed = seed, level = level)
  }
  else if (boot_mi == "miboot_pooled") {
    if (is.null(argscSEM$.resample_method)) {
      argscSEM <- c(argscSEM, .resample_method = "bootstrap")
      if (verbose) warning("the .resample_method option has been set to 'bootstrap'.")
    }
    else if (argscSEM$.resample_method != "bootstrap") {
      argscSEM$.resample_method <- "bootstrap"
      if (verbose) warning("the .resample_method option has been set to 'bootstrap'.")
    }
    res <- plssemMIBOOT_PS(model = model, data = data, m = mMI,
                           miArgs = miArgs, miPackage = pkgMI,
                           csemArgs = argscSEM, verbose = verbose,
                           seed = seed, level = level)
  }
  else if (boot_mi == "bootmi_pooled") {
    if (!is.null(argscSEM$.resample_method) && argscSEM$.resample_method != "none") {
      argscSEM$.resample_method <- "none"
      if (verbose) warning("the .resample_method option has been set to 'none'.")
    }
    res <- plssemBOOTMI_PS(model = model, data = data, m = mMI,
                           miArgs = miArgs, miPackage = pkgMI,
                           csemArgs = argscSEM, bootArgs = argsBOOT,
                           verbose = verbose, seed = seed, level = level)
  }
  else if (boot_mi == "weighted_bootmi") {
    if (!is.null(argscSEM$.resample_method) && argscSEM$.resample_method != "none") {
      argscSEM$.resample_method <- "none"
      if (verbose) warning("the .resample_method option has been set to 'none'.")
    }
    # if (argsBOOT$parallel != "no") {
    #   argsBOOT$parallel <- "no"
    #   if (verbose) warning("the weighted version of 'bootmi' can't use parallel computation and it will take longer.")
    # }
    res <- plssemWGT_BOOTMI(model = model, data = data, m = mMI,
                            miArgs = miArgs, miPackage = pkgMI,
                            csemArgs = argscSEM, bootArgs = argsBOOT,
                            verbose = verbose, seed = seed, level = level,
                            wgtType = wgtType)
  }
  else {
    stop("the specified bootstrap/MI method is not available.")
  }

  res_all$res <- res
  res_all$data <- data
  res_all$start_seed <- start_seed
  res_all$call <- CALL
  res_all$args <- args_list

  res_all
}
