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
