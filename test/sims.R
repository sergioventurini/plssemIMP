library(plssemMI)

#  (1) multivariate normal data
# mean <- rep(0, 9)
# sigma <- diag(9)
# sigma1 <- matrix(c(1, .4, .3, .4, 1, .5, .3, .5, 1), nrow = 3)
# sigma2 <- matrix(c(1, .1, .2, .1, 1, .2, .2, .2, 1), nrow = 3)
# sigma3 <- matrix(c(1, .7, .5, .7, 1, .6, .5, .6, 1), nrow = 3)
# sigma[1:3, 1:3] <- sigma1
# sigma[4:6, 4:6] <- sigma2
# sigma[7:9, 7:9] <- sigma3
# argsCD <- list(mean = mean, sigma = sigma)

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
# argsCD <- list(mean = mean, sigma = sigma,
#                skew = skew, kurt = kurt)

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
# argsCD <- list(model = model)

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

## multiple imputation analysis
res <- run_sims(runs = 50,
                argsCD = list(n = 5e2, method = "model", pkg = "cSEM.DGP",
                              model = model),
                argsMM = list(prop = .5, mech = "MCAR", method = "ampute"),
                argsMI = list(m = 5, methods = c("pmm", "norm"), pkg = "mice"),
                argscSEM = list(.disattenuate = TRUE,
                                .R = 100,
                                .tolerance = 1e-07,
                                .resample_method = "bootstrap",
                                .handle_inadmissibles = "replace"))

## analysis on original complete data
opts <- list(.disattenuate = TRUE,
             .R = 100,
             .tolerance = 1e-07,
             .resample_method = "bootstrap",
             .handle_inadmissibles = "replace")
dat_orig <- res[[1]]$dat_orig
res_orig_call <- list(cSEM::csem, .model = model, .data = dat_orig)
res_orig_call <- c(res_orig_call, opts)
res_orig <- eval(as.call(res_orig_call))
res_orig <- csem_combine(res_orig)

## complete-case (listwise deletion) analysis
opts <- list(.disattenuate = TRUE,
             .R = 100,
             .tolerance = 1e-07,
             .resample_method = "bootstrap",
             .handle_inadmissibles = "replace")
dat_miss <- na.omit(res[[1]]$dat_miss)
res_miss_call <- list(cSEM::csem, .model = model, .data = dat_miss)
res_miss_call <- c(res_miss_call, opts)
res_miss <- eval(as.call(res_miss_call))
res_miss <- csem_combine(res_miss)
