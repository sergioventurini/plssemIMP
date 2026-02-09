model_coef_vec <- function(res, estimates = FALSE, what = "all") {
  if (!estimates) {
    parsed_model <- cSEM::parseModel(res)

    # path coefficients
    path_vec <- matrix2vec(parsed_model$structural2, names = TRUE)
    path_vec <- path_vec[path_vec != 0, ]

    # loadings
    load_vec <- matrix2vec(parsed_model$measurement2, names = TRUE)
    load_vec <- load_vec[load_vec != 0, ]
  }
  else {
    # path coefficients
    path_vec <- matrix2vec(res$Estimates$Path_estimates, names = TRUE)
    path_vec <- path_vec[path_vec != 0, ]

    # loadings
    load_vec <- matrix2vec(res$Estimates$Loading_estimates, names = TRUE)
    load_vec <- load_vec[load_vec != 0, ]
  }
  names(path_vec) <- gsub("eta", "gamma", gsub("_eta", "", names(path_vec), fixed = TRUE), fixed = TRUE)
  names(load_vec) <- gsub(".*_", "", names(load_vec))

  # combine
  if (what == "all") {
    coef_vec <- c(path_vec, load_vec)
  }
  else if (what == "path") {
    coef_vec <- path_vec
  }
  else if (what == "load") {
    coef_vec <- load_vec
  }
  else {
    stop("the objects to combine are not available.")
  }

  coef_vec
}

csem_model_coef <- function(model) {
  model_parsed <- cSEM::parseModel(model)
  eta_endo <- model_parsed$cons_endo
  eta_exo <- model_parsed$cons_exo

  B <- model_parsed$structural2[eta_endo, eta_endo]
  Gamma <- model_parsed$structural2[eta_endo, eta_exo]

  Lambda_x <- model_parsed$measurement2[eta_exo, ]
  Lambda_x <- t(Lambda_x[, which(colSums(Lambda_x) > 0)])
  Lambda_y <- model_parsed$measurement2[eta_endo, ]
  Lambda_y <- t(Lambda_y[, which(colSums(Lambda_y) > 0)])

  list(B = B, Gamma = Gamma, Lambda_x = Lambda_x, Lambda_y = Lambda_y)
}

csem_to_mean <- function(model, kappa = NULL, alpha = NULL, tau_x = NULL, tau_y = NULL) {
  csem_model <- csem_model_coef(model)
  B <- csem_model$B
  Gamma <- csem_model$Gamma
  Lambda_x <- csem_model$Lambda_x
  Lambda_y <- csem_model$Lambda_y
  
  if (is.null(kappa)) {
    kappa <- matrix(0, nrow = nrow(Gamma), ncol = 1)
  }
  if (is.null(alpha)) {
    alpha <- matrix(0, nrow = nrow(Gamma), ncol = 1)
  }
  if (is.null(tau_x)) {
    tau_x <- matrix(0, nrow = nrow(Lambda_x), ncol = 1)
  }
  if (is.null(tau_y)) {
    tau_y <- matrix(0, nrow = nrow(Lambda_y), ncol = 1)
  }

  mu_x <- numeric(nrow(Lambda_x))
  mu_y <- numeric(nrow(Lambda_y))

  B_star <- solve(diag(nrow(B)) - B)
  mu_y <- tau_y + Lambda_y %*% B_star %*% (alpha + Gamma %*% kappa)
  mu_x <- tau_x + Lambda_x %*% kappa

  rbind(mu_y, mu_x)
}

csem_to_var <- function(model, Phi = NULL, Psi = NULL, Theta_e = NULL, Theta_d = NULL, Theta_de = NULL) {
  csem_model <- csem_model_coef(model)
  B <- csem_model$B
  Gamma <- csem_model$Gamma
  Lambda_x <- csem_model$Lambda_x
  Lambda_y <- csem_model$Lambda_y

  if (is.null(Psi)) {
    Psi <- diag(nrow(B))
  }
  if (is.null(Phi)) {
    Phi <- diag(nrow(Gamma))
  }
  if (is.null(Theta_d)) {
    Theta_d <- diag(nrow(Lambda_x))
  }
  if (is.null(Theta_e)) {
    Theta_e <- diag(nrow(Lambda_y))
  }
  if (is.null(Theta_de)) {
    Theta_de <- matrix(0, nrow(Lambda_x), ncol(Lambda_y))
  }
  
  B_star <- solve(diag(nrow(B)) - B)
  Sigma_11 <- Lambda_y %*% B_star %*% (Gamma %*% Phi %*% t(Gamma) + Psi) %*% t(B_star) %*% t(Lambda_y) + Theta_e
  Sigma_21 <- Lambda_y %*% B_star %*% Gamma %*% Phi %*% t(Lambda_x) + t(Theta_de)
  Sigma_12 <- Lambda_x %*% Phi %*% t(Gamma) %*% t(B_star) %*% t(Lambda_y) + Theta_de
  Sigma_22 <- Lambda_x %*% Phi %*% t(Lambda_x) + Theta_d

  rbind(cbind(Sigma_11, Sigma_12), cbind(Sigma_21, Sigma_22))
}

#' Generate the indicator correlation matrix
#'
#' Generates the indicator correlation matrix based on the parameters of a
#' structural equation model in [lavaan model syntax][lavaan::model.syntax].
#'
#' @usage generateSigma(
#'  .model                    = NULL,
#'  .handle_negative_definite = c("stop", "drop", "set_NA")
#'  )
#'
#' @param .model A model generated using [generate_model()].
#' @param .handle_negative_definite Character string. How should negative definite
#'   indicator correlation matrices be handled? One of `"stop"`, `"drop"` or `"set_NA"`
#'   in which case an `NA` is produced. Defaults to `"stop"`.
#
#' @return A K by K matrix of indicator correlations. K is the number of indicators
#'
#' @export
generate_vcov <- function(.model = NULL,
  .handle_negative_definite = c("stop", "drop", "set_NA")) {
  ## Match arguments
  .handle_negative_definite <- match.arg(.handle_negative_definite)

  .model <- generate_model(.model)

  ## Get relevant objects
  model       <- .model
  con_type    <- model$construct_type
  # Loadings
  Lambda      <- t(model$measurement)
  # Measurement errors
  Theta       <- model$error_cor
  # Correlation between exogenous (Phi) (if correlation between all constructs
  # is supplied all constructs are "exogenous)
  Phi         <- model$phi

  if (!all(model$structural == 0)) {
    vars_exo    <- model$cons_exo
    vars_endo   <- model$cons_endo
    path_matrix <- model$structural2
    # Path from exogenous to endogenous
    Gamma       <- path_matrix[vars_endo, vars_exo, drop = FALSE]
    # Path from endogenous to endogenous
    B           <- path_matrix[vars_endo, vars_endo, drop = FALSE]

    ### Checks and errors --------------------------------------------------------
    #   A maximum of 5 exogenous constructs is allowed:
    #     1. If there is 1 exogenous construct  : a maximum of 7 endogenous constructs is allowed
    #     2. If there are 2 exogenous constructs: a maximum of 6 endogenous constructs is allowed
    #     3. If there are 3 exogenous constructs: a maximum of 5 endogenous constructs is allowed
    #     4. If there are 4 exogenous constructs: a maximum of 4 endogenous constructs is allowed
    #     5. If there are 5 exogenous constructs: a maximum of 4 endogenous constructs is allowed

    if (length(vars_exo) > 5) {
      stop("Models containing more than 5 exogenous constructs are not supported.",
           call. = FALSE)
    }
    if (length(vars_endo) > 7) {
      stop("Models containing more than 7 endogenous constructs are not supported.",
           call. = FALSE)
    }
    if (length(vars_exo) > 2 && length(vars_endo) > 6) {
      stop("Models containing more than 2 exogenous AND more than 6 endogenous constructs are not supported.",
           call. = FALSE)
    }
    if (length(vars_exo) > 3 && length(vars_endo) > 5) {
      stop("Models containing more than 3 exogenous AND more than 5 endogenous constructs are not supported.",
           call. = FALSE)
    }
  }

  ## Modify and fill Lambda and Theta ------------------------------------------
  for (j in colnames(Lambda)) {
    if (con_type[j] == "Composite") {
      indicators <- colnames(model$measurement2[j, model$measurement2[j, ] != 0, drop = FALSE])
      # If j is a composite, the values of measurement2 are interpreted as weights
      w_j  <- model$measurement2[j, indicators]

      # If weights are given, the within-block indicator correlation matrix
      # must be given as well. If j is a composite values in indicator_cor
      # are interpreted as indicator correlations
      Sigma_jj <- as.matrix(model$indicator_cor[indicators, indicators])
      diag(Sigma_jj) <- 1

      if (nrow(Sigma_jj) > 1 && sum(Sigma_jj[lower.tri(Sigma_jj)]) == 0) {
        stop("Indicator correlation matrix of indicators: ",
             paste0("`", indicators, "`", collapse = ","),
             " (construct `", j, "`) is zero.\n",
             "Please specify the correlation using e.g., `",
             indicators[1],  " ~~ 0.4*", indicators[2], "`", call. = FALSE)
      }

      # Scale weigths
      ws_j <- w_j / c(sqrt(w_j %*% Sigma_jj %*% w_j))

      # Compute lambda_j
      lambda_j <- c(Sigma_jj %*% ws_j)

      # Replace corresponding elements in Lambda
      Lambda[indicators, j] <- lambda_j

      # Compute theta
      Theta[indicators, indicators] <-  Sigma_jj - lambda_j %*% t(lambda_j)


    } else {# Common factor
      indicators <- colnames(model$measurement2[j, model$measurement2[j, ] != 0, drop = FALSE])

      # Replace corresponding elements in Lambda
      lambda_j <- model$measurement2[j, indicators]
      Lambda[indicators, j] <- lambda_j

      # Get measurement error correlation
      # If j is a common factor values in indicator_cor
      # are interpreted as measurement error correlations
      Theta[indicators, indicators] <- as.matrix(model$indicator_cor[indicators, indicators])

      # Compute Theta
      if (!is.null(dim(Theta[indicators, indicators]))) {
        diag(Theta[indicators, indicators]) <-  1 - diag(lambda_j %*% t(lambda_j))
      } else {
        Theta[indicators, indicators] <- 1 - diag(lambda_j %*% t(lambda_j))
      }
    }
  }

  if (any(model$construct_order == "Second order")) {
    ## Second order model
    vars_2nd <- model$vars_2nd
    vars_attached_to_2nd <- model$vars_attached_to_2nd
    vars_not_attached_to_2nd <- model$vars_not_attached_to_2nd

    ## Lambda matrix for then "inner" model
    Lambda_2nd <- Lambda[vars_attached_to_2nd, vars_2nd, drop = FALSE]
    I <- diag(length(vars_not_attached_to_2nd))
    Lambda_inner <- cbind(
      rbind(Lambda_2nd, matrix(0, nrow = nrow(I), ncol = ncol(Lambda_2nd))),
      rbind(matrix(0, nrow = nrow(Lambda_2nd), ncol = ncol(I)), I)
    )
    rownames(Lambda_inner) <- c(vars_attached_to_2nd, vars_not_attached_to_2nd)
    colnames(Lambda_inner) <- c(vars_2nd, vars_not_attached_to_2nd)

    ## Theta matrix for the "inner" model
    Theta_2nd <- Theta[vars_attached_to_2nd, vars_attached_to_2nd, drop = FALSE]
    Theta_inner <- cbind(
      rbind(Theta_2nd, matrix(0, nrow = nrow(I), ncol = ncol(Theta_2nd))),
      rbind(matrix(0, nrow = nrow(Theta_2nd), ncol = ncol(I)), matrix(0, nrow = nrow(I), ncol = ncol(I)))
    )

    rownames(Theta_inner) <- c(vars_attached_to_2nd, vars_not_attached_to_2nd)
    colnames(Theta_inner) <- c(vars_attached_to_2nd, vars_not_attached_to_2nd)

    ## "Clean" Lambda and Theta
    selector1 <- !(rownames(Lambda) %in% vars_attached_to_2nd)
    Lambda <- Lambda[selector1, !(colnames(Lambda) %in% vars_2nd), drop = FALSE]

    Theta  <- Theta[selector1, selector1, drop = FALSE]

    vcv_construct_inner <- generate_cor(.Gamma = Gamma, .B = B, .Phi = Phi)

    vcv_construct <- Lambda_inner[colnames(Lambda), rownames(vcv_construct_inner)] %*%
      vcv_construct_inner %*% t(Lambda_inner[colnames(Lambda), rownames(vcv_construct_inner)]) +
      Theta_inner[colnames(Lambda), colnames(Lambda)]

  } else {
    # Compute the construct correlation matrix if a structural model has been
    # supplied
    if (all(model$structural == 0)) {
      vcv_construct <- Phi
    } else {
      vcv_construct <- generate_cor(.Gamma = Gamma, .B = B, .Phi = Phi)
    }
  }

  # Compute the indicator correlation matrix (Sigma)
  Sigma <- Lambda %*% vcv_construct %*% t(Lambda) + Theta

  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]

  # Check if semi-positve definite
  if (!matrixcalc::is.positive.semi.definite(Sigma)) {
    if (.handle_negative_definite %in% c("drop", "set_NA")) {
      NA
    } else if (.handle_negative_definite == "stop") {
      stop("Indicator correlation matrix is not semi-positive definite.",
           call. = FALSE)
    }
  } else {
    Sigma
  }
}

#' Internal: Compute the construct correlation matrix
#'
#' Calculate the construct correlation matrix.
#'
#' @param .Gamma A matrix containing the path coefficients from the exogenous on
#'   the endogenous constructs.
#' @param .B A matrix containing the path coefficients from the endogenous on
#'   the endogenous constructs.
#' @param .Phi A symmetric matrix of correlations between exogenous constructs
#'
#' @return A matrix of construct correlations.
#'
#' @export
generate_cor <- function(.Gamma = NULL, .B = NULL, .Phi = NULL) {
  # Number of endogenous
  k <- nrow(.B)

  # Set up empty matrix of variances of the structural error terms
  Psi <- matrix(0, nrow = nrow(.B), ncol = ncol(.B), dimnames = dimnames(.B))

 # Define gamma, beta and phi for the following calculations
  gamma_temp  <- matrix(0, nrow = 7, ncol = 5)
  beta_temp   <- matrix(0, nrow = 7, ncol = 7)
  phi_temp    <- matrix(0, nrow = 5, ncol = 5)

  beta_temp[1:nrow(.B), 1:ncol(.B)]          <- .B
  gamma_temp[1:nrow(.Gamma), 1:ncol(.Gamma)] <- .Gamma
  phi_temp[1:nrow(.Phi), 1:ncol(.Phi)]       <- .Phi

  Psi[1, 1] <- varzeta1(beta_temp, gamma_temp, phi_temp)

  if (k >= 2){
    Psi[2, 2] <- varzeta2(beta_temp, gamma_temp, phi_temp)
  }
  if (k >= 3){
    Psi[3, 3] <- varzeta3(beta_temp, gamma_temp, phi_temp)
  }
  if (k >= 4){
    Psi[4, 4] <- varzeta4(beta_temp, gamma_temp, phi_temp)
  }
  if (k >= 5){
    Psi[5, 5] <- varzeta5(beta_temp, gamma_temp, phi_temp)
  }
  if (k >= 6){
    Psi[6, 6] <- varzeta6(beta_temp, gamma_temp, phi_temp)
  }
  if (k >= 7){
    Psi[7, 7] <- varzeta7(beta_temp, gamma_temp, phi_temp)
  }

  I            <- diag(nrow(.B))
  Pi           <- solve(I - .B) %*% .Gamma
  VCV_endo_exo <- Pi %*% .Phi
  VCV_endo     <- Pi %*% .Phi %*% t(Pi) + solve(I - .B) %*% Psi %*% t(solve(I - .B))

  VCV          <- cbind(rbind(.Phi, VCV_endo_exo),rbind(t(VCV_endo_exo),VCV_endo))
  return(VCV)
}

#' Internal: generate cSEMModels
#'
#' Generate all possible [cSEMModel][cSEM::csem_model]s.
#'
#' @usage generate_model(.model, ...)
#'
#' @param .model A model in [lavaan model syntax][lavaan::model.syntax] possibly
#'   containing labels.
#' @param ... `"name" = values` pairs. `"name"` is a character string giving the
#'   label used for the parameter of interest. `values` is a numeric vector of
#'   values to use for the paramter given by `"name"`.
#
#' @return A list of [cSEMModel][cSEM::csem_model]s
#'
#' @keywords internal
generate_model <- function(.model, ...) {
  ## Collect dotdotdot (...) arguments
  params <- as.list(list(...))
  param_names <- names(params)

  xx  <- cSEM::parseModel(.model)

  if (is.null(xx$measurement2)) {
    stop("No population values given. Please specify all population values", call. = FALSE)
  }

  if (xx$model_type == "Nonlinear") {
    stop("Currently, models containing nonlinear terms are not supported.", call. = FALSE)
  }

  ss  <- xx$structural2
  m   <- xx$measurement2
  e   <- xx$indicator_cor
  cc  <- xx$construct_cor

  ## Only the Phi matrix is required (correlation matrix between exogenous constructs)
  Phi <- matrix(0,
                nrow = length(xx$cons_exo),
                ncol = length(xx$cons_exo),
                dimnames = list(xx$cons_exo, xx$cons_exo)
  )

  # Get row and column names for constructs
  row_index <- intersect(rownames(Phi), rownames(cc))
  col_index <- intersect(colnames(Phi), colnames(cc))

  Phi[row_index, col_index] <- cc[row_index, col_index, drop = FALSE]

  # Set diagonal elements to 1
  diag(Phi) <- 1

  ## Structural model ----------------------------------------------------------
  # Which elements in param_names match elements in ss?
  param_names_path <- intersect(param_names, c(ss))

  # Get the array indices for these matches
  indices <- which(matrix(ss %in% param_names_path, dim(ss)), arr.ind = TRUE)

  # Compute all combinations of the variables and create new data structural
  # models
  path_coefs <- NULL
  if (nrow(indices) > 0) {
    path_coefs <- expand.grid(params[param_names_path])

    sl <- lapply(1:nrow(path_coefs), function(x) {
      # Its crucial to order path_coef (using ss[indices])!!
      ss[indices] <- unlist(path_coefs[x, ss[indices]])
      class(ss) <- "numeric"
      ss
    })
  } else {
    class(ss) <- "numeric"
    sl <- list(ss)
    sl
  }

  ## Measurement/composite model -----------------------------------------------
  # Which elements in param_names match elements in m?
  param_names_measurement <- intersect(param_names, c(m))

  # Get the array indices for these matches
  indices <- which(matrix(m %in% param_names_measurement, dim(m)), arr.ind = TRUE)

  #  Compute all combinations of the variables and create new data structural
  # models
  measurement_coefs <- NULL
  if (nrow(indices) > 0) {
    measurement_coefs <- expand.grid(params[param_names_measurement])

    ml <- lapply(1:nrow(measurement_coefs), function(x) {
      m[indices] <- unlist(measurement_coefs[x, m[indices]])
      class(m) <- "numeric"
      m
    })
  } else {
    class(m) <- "numeric"
    ml <- list(m)
    ml
  }

  ## Indicator correlation (if composite)/ Measurement error correlation
  ## (if common factor) ---------------------------------------------
  # Which elements in param_names match elements in e?
  param_names_error <- intersect(param_names, c(e))

  # Get the array indices for these matches
  indices <- which(matrix(e %in% param_names_error, dim(e)), arr.ind = TRUE)

  #  Compute all combinations of the variables and create new data structural
  # models
  error_coefs <- NULL
  if (nrow(indices) > 0) {
    error_coefs <- expand.grid(params[param_names_error])

    el <- lapply(1:nrow(error_coefs), function(x) {
      e[indices] <- unlist(error_coefs[x, e[indices]])
      class(e) <- "numeric"
      e
    })
  } else {
    class(e) <- "numeric"
    el <- list(e)
    el
  }

  ## Structural error correlation ----------------------------------------------
  # Which elements in param_names match elements in Phi?
  param_names_Phi <- intersect(param_names, c(Phi))

  # Get the array indices for these matches
  indices <- which(matrix(Phi %in% param_names_Phi, dim(Phi)), arr.ind = TRUE)

  # Compute all combinations of the variables and create new data structural
  # models
  Phi_coefs <- NULL
  if (nrow(indices) > 0) {
    Phi_coefs <- expand.grid(params[param_names_Phi])

    Phil <- lapply(1:nrow(Phi_coefs), function(x) {
      Phi[indices] <- unlist(Phi_coefs[x, Phi[indices]])
      class(Phi) <- "numeric"
      Phi
    })
  } else {
    class(Phi) <- "numeric"
    Phil <- list(Phi)
    Phil
  }

  ## Merge
  coef_df <- NULL
  if (!is.null(path_coefs)) {
    coef_df <- path_coefs
    if (!is.null(measurement_coefs)) {
      coef_df <- merge(coef_df, measurement_coefs, sort = FALSE)
    }
    if (!is.null(error_coefs)) {
      coef_df <- merge(coef_df, error_coefs, sort = FALSE)
    }
    if (!is.null(Phi_coefs)) {
      coef_df <- merge(coef_df, Phi_coefs, sort = FALSE)
    }
  } else if (!is.null(measurement_coefs)) {
    coef_df <- measurement_coefs
    if (!is.null(error_coefs)) {
      coef_df <- merge(coef_df, error_coefs, sort = FALSE)
    }
    if (!is.null(Phi_coefs)) {
      coef_df <- merge(coef_df, Phi_coefs, sort = FALSE)
    }
  } else if (!is.null(error_coefs)) {
    coef_df <- error_coefs
    if (!is.null(Phi_coefs)) {
      coef_df <- merge(coef_df, Phi_coefs, sort = FALSE)
    }
  } else if (!is.null(Phi_coefs)) {
    coef_df <- Phi_coefs
  }

  # Add rownames as column
  if (!is.null(coef_df)) {
    coef_df <- data.frame("Id" = 1:nrow(coef_df), coef_df)
  }

  ## Combine Structural, measurement/composite and error correlation matrices
  ## and add the rest of the information
  sme <-
    unlist(lapply(Phil, function(Phi){
      unlist(lapply(el, function(e) {
        unlist(lapply(ml, function(m) {
          lapply(sl, function(s){
            l <- list(
              "structural"  = xx$structural,
              "measurement" = xx$measurement,
              "error_cor"   = xx$error_cor,
              "construct_type"  = xx$construct_type,
              "construct_order" = xx$construct_order,
              "model_type"      = xx$model_type,
              "cons_exo"        = xx$cons_exo,
              "cons_endo"       = xx$cons_endo,
              "vars_2nd"        = xx$vars_2nd,
              "vars_attached_to_2nd" = xx$vars_attached_to_2nd,
              "vars_not_attached_to_2nd" = xx$vars_not_attached_to_2nd,
              "structural2"   = s,
              "measurement2"  = m,
              "indicator_cor" = e,
              "phi"           = Phi
            )
            class(l) <- "cSEMModel"
            l
          })
        }), recursive = FALSE)
      }), recursive = FALSE)
    }), recursive = FALSE)

  # Return
  # list("Models" = sme, "Coef_df" = coef_df)
  sme[[1]]
}
