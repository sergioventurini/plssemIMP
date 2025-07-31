library(plssemIMP)
library(tidyverse)

path_to_save <- "/Users/Sergio/Documents/Dati_VENTURINI/2_Research/1_Methods/PLS-SEM_missing/paper/results"

nruns <- 2
nimp <- 2
nboot <- 3
conflev <- 0.95
ngb <- c(5, 9, 15)

# set gloabl random seed
global_seed <- 1809
set.seed(global_seed, "L'Ecuyer-CMRG")

source(file.path(path_to_save, "sims_models.R"))
nsample <- 100
consistent <- FALSE
boot_mi <- "weighted_bootmi"
miss_mech <- "MCAR"
miss_prop <- 0.05

# mice::ampute() arguments
nvars <- lapply(models_dgp, function(m) length(cSEM::parseModel(m)$indicators))
mice_patt <- lapply(nvars,
                    function(x) rbind(generate_patterns(x, 2)))
# mice_patt <- lapply(nvars,
#                     function(x) rbind(generate_patterns(x, 1),
#                                       generate_patterns(x, 2)))
mice_freq <- mapply(
  function(pt, nv)
    generate_freqs(pt,
                   reps = choose(nv, 2),
                   weights = 1),
  mice_patt, nvars, SIMPLIFY = FALSE)
# mice_freq <- mapply(
#   function(pt, nv)
#     generate_freqs(pt,
#                    reps = c(choose(nv, 1), choose(nv, 2)),
#                    weights = c(0.7, 0.3)),
#   mice_patt, nvars, SIMPLIFY = FALSE)
mice_weights <- NULL
mice_cont <- TRUE
mice_odds <- NULL
mice_type <- "TAIL"
mice_methods <- c("norm", "pmm", "rf")

argsBOOT <- list(parallel = ifelse(.Platform$OS.type == "unix", "multicore", "snow"),
                 ncpus = parallel::detectCores())

res_all <- list()
allscenarios <- 1
for (md in 1:length(models_dgp)) {
  argsMI <- list(m = nimp, pkg = "mice",
                 methods = mice_methods,
                 model = models_csem[[md]],  # same model as the DGP
                 maxit = 10,
                 blocks = make_blocks(models_csem[[md]]))
  for (n in nsample) {
    argsCD <- list(method = "model", pkg = "cSEM.DGP", model = models_dgp[[md]], n = n)
    for (mmech in miss_mech) {
      for (mprop in miss_prop) {
        argsMM <- list(prop = mprop, mech = mmech, method = "ampute",
                       patterns = mice_patt[[md]], freq = mice_freq[[md]],
                       weights = mice_weights, cont = mice_cont,
                       type = mice_type, odds = mice_odds)
        
        # generation of all the data within each simulation scenario
        alldata <- run_sims(runs = nruns,
                            argsCD = argsCD,
                            argsMM = argsMM,
                            argsMI = list(),
                            argscSEM = list(),
                            argsBOOT = list(),
                            verbose = FALSE,
                            global_seed = NULL,
                            replic_seeds = NULL,
                            boot_mi = "miboot",  # this is irrelevant when generating the data sets
                            level = conflev,
                            meanimp = FALSE,
                            knnimp = FALSE,
                            listwise = FALSE,
                            fulloriginal = FALSE,
                            wgtType = "rows",
                            datalist = NULL,
                            datamisslist = NULL,
                            runALL = FALSE)
        datafulllist <- extract_data(alldata, full = TRUE)
        datamisslist <- extract_data(alldata, full = FALSE)

        # each data is now analyzed using the different estimation methods (PLS-SEM/PLSc, bootstrap/MI)
        for (PLSc in consistent) {
          argscSEM <- list(.disattenuate = PLSc,
                           .R = nboot,
                           .tolerance = 1e-07,
                           .resample_method = "bootstrap",
                           .handle_inadmissibles = "replace",
                           .eval_plan = ifelse(.Platform$OS.type == "unix", "multicore", "multisession"))
          for (bs_mi in boot_mi) {
            scenario <- paste0("SCENARIO ", allscenarios, ": MODEL ", md,
                                                          " - N = ", n,
                                                          " - PLSc = ", ifelse(PLSc, "yes", "no"),
                                                          " - BOOT/MI = ", bs_mi,
                                                          " - MECH = ", mmech,
                                                          " - PROP = ", mprop,
                                                          "\n")
            cat(paste0(strrep("#", nchar(scenario)), "\n"))
            cat(scenario)
            cat(paste0(strrep("#", nchar(scenario)), "\n"))
            res <- run_sims(runs = nruns,
                            argsCD = argsCD,
                            argsMM = argsMM,
                            argsMI = argsMI,
                            argscSEM = argscSEM,
                            argsBOOT = argsBOOT,
                            verbose = TRUE,
                            global_seed = NULL,
                            replic_seeds = NULL,
                            boot_mi = bs_mi,
                            level = conflev,
                            meanimp = FALSE,
                            knnimp = FALSE, argsKNN = list(k = ngb),
                            listwise = FALSE,
                            fulloriginal = FALSE,
                            wgtType = "rows",
                            datalist = datafulllist,
                            datamisslist = datamisslist,
                            runALL = TRUE)

            mod <- paste0("m", md)
            cnst <- ifelse(PLSc, "PLSc", "PLS-SEM")
            # save(res, file = file.path(path_to_save,
            #      paste0("res_", mod, "_", n, "_", cnst, "_", bs_mi, "_", mmech, "_", mprop, ".RData")))
            res_all[[paste0("mod", mod, "_", n, "_", cnst, "_", bs_mi, "_", mmech, "_", mprop)]] <- res
            allscenarios <- allscenarios + 1
            cat("\n")
          }
        }
      }
    }
  }
}

# res1 <- res_all
# res2 <- res_all
# res3 <- res_all
# all.equal(res1$`modm2_100_PLS-SEM_weighted_bootmi_MCAR_0.05`$run_2$rf$path_sd,
#           res2$`modm2_100_PLS-SEM_weighted_bootmi_MCAR_0.05`$run_2$rf$path_sd)
# all.equal(res1$`modm2_100_PLS-SEM_weighted_bootmi_MCAR_0.05`$run_2$rf$path_sd,
#           res3$`modm2_100_PLS-SEM_weighted_bootmi_MCAR_0.05`$run_2$rf$path_sd)
# all.equal(res1$`modm3_100_PLS-SEM_weighted_bootmi_MCAR_0.05`$run_1$pmm$load_vcov,
#           res2$`modm3_100_PLS-SEM_weighted_bootmi_MCAR_0.05`$run_1$pmm$load_vcov)
# all.equal(res1$`modm3_100_PLS-SEM_weighted_bootmi_MCAR_0.05`$run_1$pmm$load_vcov,
#           res3$`modm3_100_PLS-SEM_weighted_bootmi_MCAR_0.05`$run_1$pmm$load_vcov)
