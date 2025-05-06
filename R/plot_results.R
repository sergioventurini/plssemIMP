plot_results <- function(res, true_coefs = NULL, methods = "all", values = "est", ylim = NULL) {
  if (is.null(true_coefs) | missing(true_coefs)) {
    true_coefs <- NA
  }

  nruns <- length(res) - 4
  nmeth <- length(res[[1]]) - 2
  methnm <- names(res[[1]])[1:(length(res[[1]]) - 2)]
  if (is.null(methods) | any(methods == "all")) {
    methods <- methnm
  }
  else {
    if (!all(methods %in% methnm))
      stop("the specified methods are not available.")
  }
  valuenm <- c("est", "sd")
  if (is.null(values) | any(values == "all")) {
    values <- valuenm
  }
  else {
    if (!all(values %in% valuenm))
      stop("the specified values to graph are not available.")
  }

  for (meth in methods) {
    for (val in values) {
      res_path_toplot <- extract_results(res, approach = meth, type = "path", what = val)
      rownames(res_path_toplot) <- paste0("eta", 1:nrow(res_path_toplot))
      res_load_toplot <- extract_results(res, approach = meth, type = "load", what = val)
      rownames(res_load_toplot) <- paste0("lambda", 1:nrow(res_load_toplot))
      res_toplot <- rbind(res_path_toplot, res_load_toplot)
      if (!is.null(true_coefs) & !missing(true_coefs) & !any(is.na(true_coefs))) {
        res_toplot <- data.frame("param" = rownames(res_toplot), res_toplot, true_coefs)
      }
      else {
        res_toplot <- data.frame("param" = rownames(res_toplot), res_toplot)
      }
      df_long <- res_toplot %>%
        tidyr::pivot_longer(cols = tidyselect::starts_with("run_"),
                            names_to = "simulation",
                            values_to = "estimate")

      if (val == "est") {
        constants <- data.frame("param" = rownames(res_toplot), true_coefs = true_coefs)
        constants <- constants %>%
          dplyr::mutate(x = as.numeric(factor(param)))
      }

      out <- ggplot2::ggplot(df_long, ggplot2::aes(x = param, y = estimate, fill = param)) +
        ggplot2::geom_boxplot(alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Boxplot of MI estimates by parameter",
                      x = "parameter",
                      y = ifelse(val == "est", "estimate", "std. dev."))
      if (!is.null(ylim) & length(ylim) == 2) {
        if (ylim[1] < ylim[2]) {
          out <- out + ggplot2::coord_cartesian(ylim = ylim)
        }
        else {
          stop("the vertical axis lower limit must be smaller than the upper one.")
        }
      }
      if (!is.null(true_coefs) & !missing(true_coefs) & !any(is.na(true_coefs)) & val == "est") {
        out <- out +
          ggplot2::geom_segment(data = constants,
                                ggplot2::aes(x = x - 0.4, xend = x + 0.4,
                                    y = true_coefs, yend = true_coefs),
                                linewidth = 2, linetype = "solid",
                                alpha = 0.4)
      }

      plot(out)
    }
  }
}
