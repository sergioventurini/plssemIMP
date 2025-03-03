run_sims <- function(runs = 10, argsCD = list(n = 1e3, method = "model", pkg = "cSEM.DGP"),
  argsMM = list(prop = .5, mech = "MCAR", method = "ampute"),
  argsMI = list(m = 5, methods = c("pmm", "norm"), pkg = "mice"),
  argscSEM = list()) {
  nCD <- argsCD$n
  methodCD <- argsCD$method
  pkgCD <- argsCD$pkg
  modelCD <- argsCD$model
  methodsMI <- argsMI$methods
  propMM <- argsMM$prop
  mechMM <- argsMM$mech
  methodMM <- argsMM$method
  mMI <- argsMI$m
  pkgMI <- argsMI$pkg

  res_all <- res <- list()
  for (run in 1:runs) {
    cat(paste0("Simulation run ", run, " of ", runs, "\n"))
    ## STEP 1: generate complete data
    dat <- create_data(n = nCD, run = run, method = methodCD, pkg = pkgCD,
                       args = list(model = modelCD))
    
    ## STEP 2: generate missing data (using the ampute() function from the mice package)
    dat_miss <- make_missing(dat, prop = .5, mech = "MCAR", method = "ampute",
      missArgs = argsMM[setdiff(names(argsMM), c("prop", "mech", "method"))],
      seed = run)
    dat_orig <- dat_miss$orig
    dat <- as.data.frame(dat_miss$amputed)
    
    ## STEP 3: perform multiple imputation
    res <- list()
    for (m in 1:length(methodsMI)) {
      res[[m]] <- plssemMI(model = modelCD, data = dat, m = mMI, miArgs = list(method = methodsMI[m]),
                           miPackage = pkgMI, seed = run, csemArgs = argscSEM)
    }
    names(res) <- methodsMI
    res$dat_orig <- dat_orig
    res$dat_miss <- dat
    res_all[[run]] <- res
  }

  res_all
}
