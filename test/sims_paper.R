library(plssemIMP)
library(future)
library(future.apply)

options(future.globals.maxSize = 600 * 1024^2)  # 600 MB

# function for logging the simulation in each scenario
log_msg <- function(msg, log_file) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg),
      file = log_file, append = TRUE)
}

# path_to_save <- "/Users/Sergio/Documents/Dati_VENTURINI/2_Research/1_Methods/PLS-SEM_missing/paper/results"
path_to_save <- "/home/ubuntu/plssem_simulation/results"

nruns <- 500
nimp <- 20
nboot <- 200
conflev <- 0.95
ngb <- c(5, 9, 13)

# set global random seed
global_seed <- 3110
set.seed(global_seed, kind = "L'Ecuyer-CMRG")

# source(file.path(path_to_save, "sims_models.R"))
source(file.path("/home/ubuntu/plssem_simulation", "sims_models.R"))
nsample <- c(200, 1000)  #c(200, 500, 1000)
consistent <- TRUE  #c(FALSE, TRUE)
boot_mi <- c("miboot", "bootmi", "weighted_bootmi")
miss_mech <- c("MCAR", "MAR")
miss_prop <- c(0.1, 0.8)  #c(0.1, 0.3, 0.8)  # prop of incomplete cases (not overall prop of missings)

# scenario table
scenario_grid <- expand.grid(
  md        = seq_along(models_dgp),
  n         = nsample,
  PLSc      = consistent,
  bs_mi     = boot_mi,
  mmech     = miss_mech,
  mprop     = miss_prop,
  stringsAsFactors = FALSE
)

# amputation arguments
nvars <- lapply(models_dgp, function(m) length(cSEM::parseModel(m)$indicators))
# mice_patt <- lapply(nvars,
#                     function(x) rbind(generate_patterns(x, 1)))
# mice_freq <- mapply(
#   function(pt, nv)
#     generate_freqs(pt,
#                    reps = choose(nv, 1),
#                    weights = 1),
#   mice_patt, nvars, SIMPLIFY = FALSE)
mice_patt <- lapply(nvars,
                    function(x) rbind(generate_patterns(x, 1),
                                      generate_patterns(x, 2)))
mice_freq <- mapply(
  function(pt, nv)
    generate_freqs(pt,
                   reps = c(choose(nv, 1), choose(nv, 2)),
                   weights = c(0.7, 0.3)),
  mice_patt, nvars, SIMPLIFY = FALSE)
mice_weights <- NULL
mice_cont <- TRUE
mice_odds <- NULL
mice_type <- "TAIL"
mice_methods <- c("norm", "pmm")  #c("norm", "pmm", "rf")

# single scenario simulation function
run_one_scenario <- function(i, scenario_grid, global_seed) {
  sc <- scenario_grid[i, ]

  # log the scenario
  log_file <- file.path(path_to_save, "logs", sprintf("scenario_%03d.log", i))
  start_time <- Sys.time()

  log_msg(
    paste0(
      "START scenario ", i,
      " | model=", sc$md,
      " n=", sc$n,
      " PLSc=", sc$PLSc,
      " boot_mi=", sc$bs_mi,
      " mech=", sc$mmech,
      " prop=", sc$mprop
    ),
    log_file
  )
  
  # set scenario-specific random seed
  # set.seed(global_seed + i, kind = "L'Ecuyer-CMRG")  # this potentially overwrites the RNG seed propagation
                                                       # set by future

  md    <- sc$md
  n     <- sc$n
  PLSc  <- sc$PLSc
  bs_mi <- sc$bs_mi
  mmech <- sc$mmech
  mprop <- sc$mprop

  outfile <- file.path(
    path_to_save,
    paste0("res_",
           "m", md, "_", n, "_",
           ifelse(PLSc, "PLSc", "PLS-SEM"), "_",
           bs_mi, "_", mmech, "_", mprop, ".RData")
  )
  if (file.exists(outfile)) {
    log_msg("Result file already exists — skipping scenario.", log_file)
    return(invisible(TRUE))
  }
  
  tryCatch({
    argsCD <- list(
      method = "model",
      pkg = "cSEM.DGP",
      model = models_dgp[[md]],
      n = n
    )

    argsMM <- list(
      prop = mprop,
      mech = mmech,
      method = "ampute",
      patterns = mice_patt[[md]],
      freq = mice_freq[[md]],
      weights = mice_weights,
      cont = mice_cont,
      type = mice_type,
      odds = mice_odds
    )

    argsMI <- list(
      m = nimp,
      pkg = "mice",
      methods = mice_methods,
      model = models_csem[[md]],
      maxit = 10,
      blocks = make_blocks(models_csem[[md]])
    )

    argscSEM <- list(
      .disattenuate = PLSc,
      .R = nboot,
      .tolerance = 1e-07,
      .resample_method = "bootstrap",
      .handle_inadmissibles = "replace",
      # .handle_inadmissibles = "drop",
      .eval_plan = "sequential"
    )

    argsBOOT <- list(parallel = "no", ncpus = 1)

    log_msg(
      paste0(
        "Starting estimation: nruns=", nruns,
        ", nboot=", nboot,
        ", nimp=", nimp
      ),
      log_file
    )

    res <- run_sims(
      runs = nruns,
      argsCD = argsCD,
      argsMM = argsMM,
      argsMI = argsMI,
      argscSEM = argscSEM,
      argsBOOT = argsBOOT,
      verbose = FALSE,
      boot_mi = bs_mi,
      level = conflev,
      meanimp = TRUE,
      knnimp = TRUE, argsKNN = list(k = ngb),
      listwise = TRUE,
      fulloriginal = TRUE,
      wgtType = "rows",
      datalist = NULL,
      datamisslist = NULL,
      store_data = FALSE,     # the data sets are NOT stored/saved
      runALL = TRUE,
      log_file = log_file
    )

    mod  <- paste0("m", md)
    cnst <- ifelse(PLSc, "PLSc", "PLS-SEM")

    log_msg("Estimation completed. Saving results.", log_file)

    save(
      res,
      file = file.path(
        path_to_save,
        paste0("res_", mod, "_", n, "_", cnst, "_", bs_mi, "_", mmech, "_", mprop, ".RData"))
      )

    elapsed <- difftime(Sys.time(), start_time, units = "mins")

    log_msg(
      sprintf(
        "END scenario %d | elapsed time: %.1f minutes",
        i, as.numeric(elapsed)
      ),
      log_file
    )
  }, error = function(e) {
    log_msg(paste("ERROR:", conditionMessage(e)), log_file)
    stop(e)
  })

  invisible(TRUE)
}

# run_one_scenario(1, scenario_grid, global_seed)  # miboot
# run_one_scenario(5, scenario_grid, global_seed)  # bootmi
# run_one_scenario(9, scenario_grid, global_seed)  # weighted_bootmi

# parallel simulation of all scenarios
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

plan(sequential)   # reset any existing plan
# nworkers <- 4      # it does not work with more then 4 cores on my laptop :(
nworkers <- max(1, parallelly::availableCores() - 2)
plan(multisession, workers = nworkers)

future_lapply(
  seq_len(nrow(scenario_grid)),
  function(i) {
    run_one_scenario(
      i = i,
      scenario_grid = scenario_grid,
      global_seed = global_seed
    )
  },
  future.seed = TRUE
)

plan(sequential)   # reset the existing plan
