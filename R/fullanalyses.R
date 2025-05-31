fulldata <- function(model, data, csemArgs = list(), verbose = FALSE, level = 0.95, ...) {
  CALL <- match.call()
  dots <- list(...)

  if (missing(data)) {
    stop("a dataset is needed to run the fulldata() function.")
  }

  csemListCall <- list(cSEM::csem, .model = model, .data = data)
  csemListCall <- c(csemListCall, csemArgs)
  csemListCall <- csemListCall[!duplicated(names(csemListCall))]
  fit <- list()
  fit$FitList <- list(suppressWarnings(eval(as.call(csemListCall))))
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
  fit$DataList <- data
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
