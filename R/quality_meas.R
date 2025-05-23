rawbias <- function(estimates, true) {
  rowMeans(estimates) - true
}

percentbias <- function(estimates, true) {
  100*rawbias(estimates, true)/true
}

coveragerate <- function(ciestimates, true) {
  lower <- as.matrix(ciestimates[, , 1])
  upper <- as.matrix(ciestimates[, , 2])
  rowMeans(lower < true & true < upper)
}

averagewidth <- function(ciestimates) {
  rowMeans(as.matrix(ciestimates[, , 2]) - as.matrix(ciestimates[, , 1]))
}

rmse <- function(estimates, true) {
  sqrt(rowMeans((estimates - true)^2))
}

qualitymeasures <- function(estimates, true) {
  data.frame(RB = rawbias(estimates, true),
    PB = percentbias(estimates, true),
    CR = coveragerate(estimates, true),
    AW = averagewidth(estimates),
    RMSE = rmse(estimates, true))
}
