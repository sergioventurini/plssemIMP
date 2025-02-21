plssemMI <- function(model, data, fun = "csem", ..., m, miArgs = list(),
                     miPackage = "Amelia", seed = NULL) {
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
  if (!is.null(seed))
    seed <- as.integer(seed[1])
  imputedData <- NULL
  if (missing(data)) {
  }
  else if (is.data.frame(data)) {
    if (miPackage[1] == "Amelia") {
      requireNamespace("Amelia")
      if (!"package:Amelia" %in% search()) 
        attachNamespace("Amelia")
      imputeCall <- c(list(Amelia::amelia, x = data, m = m, 
                           p2s = 0), miArgs)
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
          imputedData[[i]] <- mice::complete(data, action = i, 
                                             include = FALSE)
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
  .getOutput. <- function(obj) {
    converged <- lavaan::lavInspect(obj, "converged")
    if (converged) {
      se <- lavaan::parTable(obj)$se
      se.test <- all(!is.na(se)) & all(se >= 0) & any(se != 
                                                        0)
      if (lavaan::lavInspect(obj, "ngroups") == 1L && lavaan::lavInspect(obj, 
                                                                         "nlevels") == 1L) {
        Heywood.lv <- det(lavaan::lavInspect(obj, "cov.lv")) <= 
          0
        Heywood.ov <- det(lavaan::lavInspect(obj, "theta")) <= 
          0
      }
      else {
        Heywood.lv <- !all(sapply(lavaan::lavInspect(obj, 
                                                     "cov.lv"), det) > 0)
        Heywood.ov <- !all(sapply(lavaan::lavInspect(obj, 
                                                     "theta"), det) > 0)
      }
    }
    else {
      se.test <- Heywood.lv <- Heywood.ov <- NA
    }
    list(sampstat = lavaan::lavInspect(obj, "sampstat"), 
         coefMats = lavaan::lavInspect(obj, "est"),
         satPT = data.frame(lavaan::lav_partable_unrestricted(obj),
          stringsAsFactors = FALSE), modindices = try(lavaan::modindices(obj),
          silent = TRUE), cov.lv = lavaan::lavInspect(obj, "cov.lv"),
          converged = converged, SE = se.test, Heywood.lv = Heywood.lv,
          Heywood.ov = Heywood.ov)
  }
  lavListCall <- list(lavaan::lavaanList, model = model, dataList = imputedData, 
                      cmd = fun)
  lavListCall <- c(lavListCall, dots)
  lavListCall$store.slots <- c("partable", "vcov", "test", 
                               "h1", "baseline")
  lavListCall$FUN <- if (is.null(dots$FUN)) 
    .getOutput.
  else function(obj) {
    temp1 <- .getOutput.(obj)
    temp2 <- dots$FUN(obj)
    if (!is.list(temp2)) 
      temp2 <- list(userFUN1 = temp2)
    if (is.null(names(temp2))) 
      names(temp2) <- paste0("userFUN", 1:length(temp2))
    duplicatedNames <- which(sapply(names(temp2), function(x) {
      x %in% c("sampstat", "coefMats", "satPT", "modindices", 
               "converged", "SE", "Heywood.lv", "Heywood.ov", 
               "cov.lv")
    }))
    for (i in duplicatedNames) names(temp2)[i] <- paste0("userFUN", 
                                                         i)
    c(temp1, temp2)
  }
  fit <- eval(as.call(lavListCall))
  fit@SampleStatsList <- lapply(fit@funList, "[[", i = "sampstat")
  fit@DataList <- imputedData
  for (i in 1:m) fit@h1List[[i]] <- c(fit@h1List[[i]], list(PT = fit@funList[[i]]$satPT))
  fit <- as(fit, "lavaan.mi")
  fit@coefList <- lapply(fit@funList, "[[", i = "coefMats")
  fit@miList <- lapply(fit@funList, "[[", i = "modindices")
  fit@phiList <- lapply(fit@funList, "[[", i = "cov.lv")
  fit@seed <- seed
  fit@call <- CALL
  fit@lavListCall <- lavListCall
  fit@imputeCall <- imputeCall
  convList <- lapply(fit@funList, "[", i = c("converged", "SE", 
                                             "Heywood.lv", "Heywood.ov"))
  nonConv <- which(sapply(convList, is.null))
  if (length(nonConv)) 
    for (i in nonConv) {
      convList[[i]] <- list(converged = FALSE, SE = NA, 
                            Heywood.lv = NA, Heywood.ov = NA)
    }
  fit@convergence <- lapply(convList, function(x) do.call(c, 
                                                          x))
  conv <- which(sapply(fit@convergence, "[", i = "converged"))
  if (!length(conv)) 
    warning("The model did not converge for any imputed data sets.")
  funNames <- names(fit@funList[[1]])
  keepIndex <- which(!sapply(funNames, function(x) {
    x %in% c("sampstat", "coefMats", "satPT", "modindices", 
             "converged", "SE", "Heywood.lv", "Heywood.ov", "cov.lv")
  }))
  if (length(keepIndex)) {
    fit@funList <- lapply(fit@funList, "[", i = keepIndex)
    if (length(keepIndex) > 1L) {
      keepNames <- funNames[keepIndex]
      noNames <- which(keepNames == "")
      for (i in seq_along(noNames)) keepNames[noNames[i]] <- paste0("userFUN", 
                                                                    i)
      fit@funList <- lapply(fit@funList, "names<-", value = keepNames)
    }
  }
  else fit@funList <- list()
  NewStartVals <- try(getMethod("coef", "lavaan.mi")(fit, type = "user", 
                                                     labels = FALSE), silent = TRUE)
  if (!inherits(NewStartVals, "try-error")) 
    fit@ParTable$start <- NewStartVals
  fit
}
