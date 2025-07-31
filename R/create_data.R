create_data <- function(n = 50, method = "norm", pkg = NULL, args = list(),
  seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed = seed, "L'Ecuyer-CMRG")
  }

  if (method == "reg") {
    beta <- args$beta
    sigma2 <- args$sigma2
    p <- length(beta)
    x <- rnorm(n)
    y <- list()
    for (j in 1:p) {
      y[[j]] <- beta[j]*x + rnorm(n, sd = sqrt(sigma2[j]))
    }
    data <- cbind(x = x, y = y)
    cn <- c("x", paste("y", 1:p))
  }
  else if (method == "norm") {
    mean <- args$mean
    sigma <- args$sigma
    p <- length(mean)
    data <- mvtnorm::rmvnorm(n = n, mean = mean, sigma = sigma)
    cn <- paste0("y", 1:p)
  }
  else if (method == "vm") {
    mean <- args$mean
    sigma <- args$sigma
    skew <- args$skew
    kurt <- args$kurt
    p <- length(mean)
    data <- SimDesign::rValeMaurelli(n = n, mean = mean, sigma = sigma,
      skew = skew, kurt = kurt)
    cn <- paste0("y", 1:p)
  }
  else if (method == "model") {
    model <- args$model
    if (pkg == "cSEM.DGP") {
      skew <- args$skew
      kurt <- args$kurt
      data <- cSEM.DGP::generateData(.model = model, .N = n,
        .skewness = skew, .kurtosis = kurt)
    }
    else if (pkg == "simstandard") {
      data <- simstandard::sim_standardized(m = model, n = n,
        observed = TRUE, latent = FALSE, errors = FALSE,
        factor_scores = FALSE, composites = FALSE, matrices = FALSE)
    }
    else {
      stop("the specified data generation package is not available.")
    }
  }
  else {
    stop("the specified data creation method is not available.")
  }

  if (!is.null(args$cn)) {
    colnames(data) <- args$cn
  }
  else {
    if (is.null(colnames(data))) {
      colnames(data) <- cn
    }
  }

  data <- as.data.frame(data)
  data
}
