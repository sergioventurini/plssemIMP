plssemMI <- function(model, data, ..., m = 5, miArgs = list(),
                     miPackage = "mice", seed = 12345) {
  CALL <- match.call()
  dots <- list(...)
  if (all(!is.null(dots$test),
    tolower(dots$test) %in% c("boot", "bootstrap", "bollen.stine")) ||
    all(!is.null(dots$se), tolower(dots$se) %in% c("boot", "bootstrap"))) {
    stop("Bootstraping unavailable (and not recommended) in combination with ", 
         "multiple imputations. For robust confidence intervals of indirect", 
         " effects, see the ?semTools::monteCarloMed help page. To bootstrap ", 
         "within each imputation, users can pass a custom function to the ", 
         "FUN= argument (see ?lavaanList) to save bootstrap distributions in ", 
         "the @funList slot, then manually combine afterward.")
  }
  if (!is.null(seed)) {
    seed <- as.integer(seed[1])
  }
  else stop("A seed is required to let you reproduce the simulation results")
  imputedData <- NULL
  if (missing(data)) {
    stop("A seed is needed to reproduce the simulation results")
  }
  else if (is.data.frame(data)) {
    if (miPackage[1] == "Amelia") {
      requireNamespace("Amelia")
      if (!"package:Amelia" %in% search()) 
        attachNamespace("Amelia")
      imputeCall <- c(list(Amelia::amelia, x = data, m = m, p2s = 0), miArgs)
      set.seed(seed)
      imputedData <- unclass(eval(as.call(imputeCall))$imputations)
    }
    else if (miPackage[1] == "mice") {
      requireNamespace("mice")
      if (!"package:mice" %in% search()) 
        attachNamespace("mice")
      imputeCall <- c(list(mice::mice, data = data, m = m, 
                           diagnostics = FALSE, printFlag = FALSE), miArgs)
      set.seed(seed)
      miceOut <- eval(as.call(imputeCall))
      imputedData <- list()
      for (i in 1:m) {
        imputedData[[i]] <- mice::complete(data = miceOut, 
                                           action = i, include = FALSE)
      }
    }
    else stop("Currently runMI only supports imputation by Amelia or mice")
  }
  else if (is.list(data)) {
    if (requireNamespace("mice", quietly = TRUE)) {
      if (mice::is.mids(data)) {
        m <- data$m
        imputedData <- list()
        for (i in 1:m) {
          imputedData[[i]] <- mice::complete(data, action = i, include = FALSE)
        }
        imputeCall <- list()
      }
      else {
        seed <- integer(length = 0)
        imputeCall <- list()
        imputedData <- data
        m <- length(data)
        class(imputedData) <- "list"
      }
    }
    else {
      seed <- integer(length = 0)
      imputeCall <- list()
      imputedData <- data
      m <- length(data)
      class(imputedData) <- "list"
    }
  }
  else if (is(data, "lavaan.mi")) {
    seed <- data@seed
    imputeCall <- data@imputeCall
    imputedData <- data@DataList
    m <- length(imputedData)
  }
  else stop("data is not a valid input type: a partially observed data.frame,", 
            " a list of imputed data.frames, or previous lavaan.mi object")
  lavListCall <- list(cSEM::csem, .model = model, .data = imputedData)
  lavListCall <- c(lavListCall, dots)
  fit <- list()
  fit$FitList <- eval(as.call(lavListCall))
  fit$ParamList <- lapply(fit$FitList,
    function(x) {
      path <- x$Estimates$Path_estimates
      path <- path[, colSums(path) != 0]
      path <- colSums(path)
      load <- colSums(x$Estimates$Loading_estimates)
      c(path, load)
    })
  fit$VCOVList <- lapply(fit$FitList,
    function(x) {
      vcov_path <- cov(x$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled)
      vcov_load <- cov(x$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled)
      combine_vcov(vcov_path, vcov_load)
    })
  fit$DataList <- imputedData
  fit$call <- CALL
  fit$lavListCall <- lavListCall
  fit$imputeCall <- imputeCall
  fit$convList <- lapply(fit$FitList, function(x) x$Information$Weight_info$Convergence_status)
  if (!all(unlist(fit$convList))) 
    warning("the model did not converge for any imputed data sets.")
  fit$dots <- dots
  fit
}
