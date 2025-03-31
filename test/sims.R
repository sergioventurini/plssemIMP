library(plssemMI)
library(tidyverse)

## setup the data creation
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
argsCD <- list(n = nsample, method = "model", pkg = "cSEM.DGP", model = model)

true_model <- cSEM::parseModel(model)
true_path <- colSums(true_model$structural2)
true_path <- true_path[true_path != 0]
true_load <- colSums(true_model$measurement2)
true_coefs <- c(true_path, true_load)

## perform imputation analysis
# set.seed(101)
nruns <- 100
nsample <- 1e3
nimp <- 5
nboot <- 200
conflev <- 0.95
argscSEM <- list(.disattenuate = TRUE,
                 .R = nboot,
                 .tolerance = 1e-07,
                 .resample_method = "bootstrap",
                 .handle_inadmissibles = "replace",
                 .eval_plan = ifelse(.Platform$OS.type == "unix", "multicore", "multisession"))
res <- run_sims(runs = nruns,
                argsCD = argsCD,
                argsMM = list(prop = .5, mech = "MCAR", method = "ampute"),
                argsMI = list(m = nimp, methods = c("pmm", "norm"), pkg = "mice",
                              model = model),  # WE ARE USING THE SAME MODEL AS IN THE DGP!
                argscSEM = argscSEM,
                argsBOOT = list(parallel = ifelse(.Platform$OS.type == "unix", "multicore", "snow"),
                                ncpus = parallel::detectCores()),
                verbose = TRUE, boot_mi = "miboot", level = conflev,
                meanimp = TRUE, knnimp = TRUE, argsKNN = list(k = c(5, 7)),
                listwise = TRUE, fulloriginal = TRUE,
                seed = 1406)

## aggregate all results
res_df <- aggregate_results(res, true_coefs = true_coefs,
                            methods = c("pmm", "listwise", "fulloriginal"),
                            qual_meas = c("PB", "CR"))

## represent the results graphically
res_path_est <- extract_results(res, approach = "pmm", type = "path", what = "est")
res_load_est <- extract_results(res, approach = "pmm", type = "load", what = "est")
mi_est <- rbind(res_path_est, res_load_est)
mi_est <- data.frame(param = rownames(mi_est), mi_est, true_coef)
df_long <- mi_est %>%
  pivot_longer(cols = starts_with("run_"),
               names_to = "simulation",
               values_to = "estimate")

constants <- data.frame(param = rownames(mi_est), true_coefs = true_coefs)
constants <- constants %>%
  mutate(x = as.numeric(factor(param)))

ggplot(df_long, aes(x = param, y = estimate, fill = param)) +
  geom_boxplot(alpha = 0.7) +
  geom_segment(data = constants,
               aes(x = x - 0.4, xend = x + 0.4,
                   y = true_coefs, yend = true_coefs),
               linewidth = 2, linetype = "solid",
               alpha = 0.4) +
  theme_minimal() +
  labs(title = "Boxplot of MI estimates by parameter",
       x = "parameter",
       y = "estimate")
