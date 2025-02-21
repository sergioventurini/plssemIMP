test_impute <- function(data, m = 5, method = "norm", ...) {
  imp <- mice(data, method = method, m = m, print = FALSE, ...)
  fit <- with(imp, lm(y ~ x))
  tab <- summary(pool(fit), "all", conf.int = TRUE)
  as.numeric(tab["x", c("estimate", "2.5 %", "97.5 %")])
}
