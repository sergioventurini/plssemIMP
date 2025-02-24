# res <- run_sims(10)

seed <- 1406

## STEP 1: generate complete data

n <- 1e4

#  (1) multivariate normal data
# mean <- rep(0, 9)
# sigma <- diag(9)
# sigma1 <- matrix(c(1, .4, .3, .4, 1, .5, .3, .5, 1), nrow = 3)
# sigma2 <- matrix(c(1, .1, .2, .1, 1, .2, .2, .2, 1), nrow = 3)
# sigma3 <- matrix(c(1, .7, .5, .7, 1, .6, .5, .6, 1), nrow = 3)
# sigma[1:3, 1:3] <- sigma1
# sigma[4:6, 4:6] <- sigma2
# sigma[7:9, 7:9] <- sigma3
# dat <- create_data(n = n, run = seed, method = "norm",
#                    args = list(mean = mean, sigma = sigma))
# pairs(dat)

#  (2) multivariate non-normal data (Vale-Maurelli approach)
# mean <- rep(0, 9)
# sigma <- diag(9)
# sigma1 <- matrix(c(1, .4, .3, .4, 1, .5, .3, .5, 1), nrow = 3)
# sigma2 <- matrix(c(1, .1, .2, .1, 1, .2, .2, .2, 1), nrow = 3)
# sigma3 <- matrix(c(1, .7, .5, .7, 1, .6, .5, .6, 1), nrow = 3)
# sigma[1:3, 1:3] <- sigma1
# sigma[4:6, 4:6] <- sigma2
# sigma[7:9, 7:9] <- sigma3
# skew <- rep(c(1.5, 1.5, 0.5), 3)
# kurt <- rep(c(3.75, 3.5, 0.5), 3)
# dat <- create_data(n = n, run = seed, method = "vm",
#                    args = list(mean = mean, sigma = sigma,
#                                skew = skew, kurt = kurt))
# pairs(dat)

# (3) using a specified SEM model through the simstandard package
# model <- '
# # Structural model
# eta3 ~ 0.4*eta1 + 0.35*eta2
# eta4 ~ 0.7*eta3
# 
# # Measurement model
# eta1 =~ 0.8*y11 + 0.9*y12 + 0.8*y13
# eta2 =~ 0.7*y21 + 0.7*y22 + 0.9*y23
# eta3 =~ 0.7*y31 + 0.8*y32 + 0.7*y33
# eta4 =~ 0.7*y41 + 0.7*y42 + 0.6*y43
# 
# # Within block indicator correlation of eta1
# y11 ~~ 0.2*y12
# y11 ~~ 0.1*y13
# y12 ~~ 0.2*y13
# 
# # Within block indicator correlation of eta2
# y21 ~~ 0.3*y22
# y21 ~~ 0.1*y23
# y22 ~~ 0.2*y23
# 
# # Within block indicator correlation of eta3
# y31 ~~ -0.2*y32
# y31 ~~ -0.3*y33
# y32 ~~ -0.5*y33
# 
# # Within block indicator correlation of eta4
# y41 ~~ 0.1*y42
# y41 ~~ 0.4*y43
# y42 ~~ 0.2*y43
# '
# 
# dat <- create_data(n = n, run = seed, method = "model", pkg = "simstandard",
#                    args = list(model = model,
#                                cn = c("y11", "y12", "y13",
#                                       "y21", "y22", "y23",
#                                       "y31", "y32", "y33",
#                                       "y41", "y42", "y43")))
# pairs(dat)

# (4) using a specified SEM model through the cSEM.DGP package
model <- '
# Structural model
eta3 ~ 0.4*eta1 + 0.35*eta2
eta4 ~ 0.7*eta3

# Measurement model
eta1 =~ 0.8*y11 + 0.9*y12 + 0.8*y13
eta2 =~ 0.7*y21 + 0.7*y22 + 0.9*y23
eta3 =~ 0.7*y31 + 0.8*y32 + 0.7*y33
eta4 =~ 0.7*y41 + 0.7*y42 + 0.6*y43
'

dat <- create_data(n = n, run = seed, method = "model", pkg = "cSEM.DGP",
                   args = list(model = model,
                               cn = c("y11", "y12", "y13",
                                      "y21", "y22", "y23",
                                      "y31", "y32", "y33",
                                      "y41", "y42", "y43")))
# pairs(dat)

## STEP 2: generate missing data (using the ampute() function from the mice package)

dat <- make_missing(dat, prop = .5)
dat_orig <- dat$orig
dat <- dat$amputed
# sapply(dat, function(x) sum(is.na(x)))
# sum(complete.cases(dat)) # or nrow(na.omit(dat))
# pairs(dat)

## STEP 3: perform multiple imputation
res <- plssemMI(model = model, data = dat, m = 5, miArgs = list(),
                miPackage = "mice", seed = seed,
                # .approach_2ndorder            = "2stage", # c("2stage", "mixed"),
                # .approach_cor_robust          = "none", # c("none", "mcd", "spearman"),
                # .approach_nl                  = "sequential", # c("sequential", "replace"),
                # .approach_paths               = "OLS",
                # .approach_weights             = "PLS-PM",
                # .conv_criterion               = "diff_absolute",
                .disattenuate                 = TRUE,
                # .dominant_indicators          = NULL,
                # .estimate_structural          = TRUE,
                # .id                           = NULL,
                # .instruments                  = NULL,
                # .iter_max                     = 100,
                # .normality                    = FALSE,
                # .PLS_approach_cf              = "dist_squared_euclid",
                # .PLS_ignore_structural_model  = FALSE,
                # .PLS_modes                    = NULL,
                # .PLS_weight_scheme_inner      = "path",
                # .reliabilities                = NULL,
                # .starting_values              = NULL,
                .resample_method              = "bootstrap", # c("none", "bootstrap", "jackknife"),
                # .resample_method2             = "none",
                .R                            = 99,
                # .handle_inadmissibles         = "drop",
                # .user_funs                    = NULL,
                # .eval_plan                    = "sequential",
                # .seed                         = NULL,
                # .sign_change_option           = "none",
                .tolerance                    = 1e-07)

## STEP 4: apply Rubin's combination rule
parest_mi <- rubin_parest(res$ParamList)
vcov_mi <- rubin_vcov(res$ParamList, res$VCOVList)
sd_mi <- sqrt(diag(vcov_mi))

## STEP 5: analysis on original complete data
opts <- res$FitList$Data_1$Information$Arguments
opts$.data <- opts$.model <- NULL
opts$.resample_method <- "bootstrap"
opts$.R <- 99
res_orig_call <- list(cSEM::csem, .model = model, .data = dat_orig)
res_orig_call <- c(res_orig_call, opts)
res_orig <- eval(as.call(res_orig_call))
colSums(res_orig$Estimates$Path_estimates)
colSums(res_orig$Estimates$Loading_estimates)
sqrt(diag(cov(res_orig$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled)))
sqrt(diag(cov(res_orig$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled)))
