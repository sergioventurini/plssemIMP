###############################################################################
# creating the master file with all results
###############################################################################

library(plssemIMP)
library(tidyverse)

path_to_save <- "/Users/Sergio/Documents/Dati_VENTURINI/2_Research/1_Methods/PLS-SEM_missing/paper/results"
path_figure <- "/Users/Sergio/Documents/Dati_VENTURINI/2_Research/1_Methods/PLS-SEM_missing/paper/Figs"

nruns <- 500
meth_nm <- c("norm", "pmm", "mean", "knn_5", "knn_9", "knn_13", "listwise")
pars_nm <- c("path_est", "path_sd", "path_lower", "path_upper", "load_est", "load_sd", "load_lower", "load_upper")
path_nm <- c("gamma21", "gamma31", "gamma32")
load_nm <- c("x11", "x12", "x13", "x21", "x22", "x23", "x31", "x32", "x33")
mech_nm <- c("MCAR", "MAR")
miboot_nm <- c("bootmi", "miboot")
param_type_levels <- c("path", "loading")
model_nm <- c("m1", "m2")
nsample <- c(200, 1000)
missprop <- c(.1, .8)
metrics <- c("RB", "PB", "CR", "RMSE", "AW")

# convert (one run × one method) into a data frame
extract_run_method <- function(
  res,
  run_id,
  method,
  scenario_id,
  md,
  nsize,
  miboot,
  mech,
  missprop
) {
  block <- res[[run_id]][[method]]
  
  ## paths
  paths <- tibble(
    scenario   = scenario_id,
    model      = md,
    nsize      = nsize,
    miboot     = miboot,
    mech       = mech,
    mprop      = missprop,
    run        = run_id,
    method     = method,
    param_type = "path",
    parameter  = names(block$path_est),
    est        = as.numeric(block$path_est),
    se         = as.numeric(block$path_sd),
    lower      = as.numeric(block$path_lower),
    upper      = as.numeric(block$path_upper)
  )
  
  ## loadings
  loads <- tibble(
    scenario   = scenario_id,
    model      = md,
    nsize      = nsize,
    miboot     = miboot,
    mech       = mech,
    mprop      = missprop,
    run        = run_id,
    method     = method,
    param_type = "loading",
    parameter  = names(block$load_est),
    est        = as.numeric(block$load_est),
    se         = as.numeric(block$load_sd),
    lower      = as.numeric(block$load_lower),
    upper      = as.numeric(block$load_upper)
  )
  
  bind_rows(paths, loads)
}

# scale up inside one scenario file
extract_scenario <- function(res, scenario_id) {
  runs    <- seq_along(res)
  methods <- names(res[[1]])[match(meth_nm, names(res[[1]]))]
  idx <- gregexpr("_", scenario_id, fixed = TRUE)[[1]]
  
  map_dfr(
    runs,
    function(r) {
      map_dfr(
        methods,
        ~ extract_run_method(
          res         = res,
          run_id      = r,
          method      = .x,
          scenario_id = scenario_id,
          md          = substr(scenario_id, start = idx[1] + 1, stop = idx[2] - 1),
          nsize       = as.integer(substr(scenario_id, start = idx[2] + 1, stop = idx[3] - 1)),
          miboot      = substr(scenario_id, start = idx[4] + 1, stop = idx[5] - 1),
          mech        = substr(scenario_id, start = idx[5] + 1, stop = idx[6] - 1),
          missprop    = as.numeric(substr(scenario_id, start = idx[6] + 1, stop = nchar(scenario_id)))
        )
      )
    }
  )
}

# stack all 32 scenarios
files <- list.files(
  path = path_to_save,
  pattern = "res_m.*\\.RData$",
  full.names = TRUE
)

master_df <- map_dfr(
  files,
  function(f) {
    load(f)                 # loads `res`
    print(f)
    res_tmp <- res[1:nruns]
    scenario_id <- tools::file_path_sans_ext(basename(f))
    extract_scenario(res_tmp, scenario_id)
  }
)

# make parameter names consistent
param_map <- tribble(
  ~param_type, ~raw,          ~canonical,
  "path",     "eta2_eta1",    "gamma21",
  "path",     "eta3_eta1",    "gamma31",
  "path",     "eta3_eta2",    "gamma32",
  
  "loading",  "eta1_x11",     "x11",
  "loading",  "eta1_x12",     "x12",
  "loading",  "eta1_x13",     "x13",
  "loading",  "eta2_x21",     "x21",
  "loading",  "eta2_x22",     "x22",
  "loading",  "eta2_x23",     "x23",
  "loading",  "eta3_x31",     "x31",
  "loading",  "eta3_x32",     "x32",
  "loading",  "eta3_x33",     "x33"
)

master_df <- master_df %>%
  left_join(
    param_map,
    by = c("param_type", "parameter" = "raw")
  ) %>%
  mutate(
    parameter = if_else(
      is.na(canonical),
      parameter,   # keep original if already canonical
      canonical
    )
  ) %>%
  select(-canonical)

# master_df_std %>%
#   distinct(param_type, parameter)

# add true parameter values
true_path <- tibble(
  model      = NA_character_,
  param_type = "path",
  parameter  = c("gamma21", "gamma31", "gamma32"),
  true_value = c(0.5, 0.3, 0.7)
)
true_load_m1 <- tibble(
  model = "m1",
  param_type = "loading",
  parameter  = c("x11", "x12", "x13", "x21", "x22", "x23", "x31", "x32", "x33"),
  true_value = c(.9, .8, .7, .7, .7, .7, .8, .8, .7)
)
true_load_m2 <- tibble(
  model = "m2",
  param_type = "loading",
  parameter  = c("x11", "x12", "x13", "x21", "x22", "x23", "x31", "x32", "x33"),
  true_value = c(.37, .24, .35, .13, .39, .25, .15, .35, .13)
)
true_table <- bind_rows(
  true_path,
  true_load_m1,
  true_load_m2
)

master_df <- master_df %>%
  left_join(
    true_table,
    by = c("model", "param_type", "parameter")
  ) %>%
  left_join(
    true_table %>% filter(is.na(model)),
    by = c("param_type", "parameter"),
    suffix = c("", "_path")
  ) %>%
  mutate(
    true_value = coalesce(true_value, true_value_path)
  ) %>%
  select(-model_path) %>%
  select(-true_value_path)

master_df <- master_df %>%
  mutate(
    mech       = factor(mech, levels = mech_nm),
    miboot     = factor(miboot, levels = miboot_nm),
    model      = factor(model, levels = model_nm),
    method     = factor(method, levels = meth_nm),
    param_type = factor(param_type, levels = param_type_levels),
    parameter  = factor(parameter, levels = param_map$canonical)
  )

## preliminary computations
master_df <- master_df %>%
  mutate(
    bias      = est - true_value,
    perc_bias = 100 * bias / true_value,
    sq_error  = (est - true_value)^2,
    covered   = lower <= true_value & upper >= true_value,
    ci_width  = upper - lower
  )

# ## raw bias
# rb_df <- master_df %>%
#   group_by(
#     model, nsize, mech, mprop, miboot,
#     method, param_type, parameter
#   ) %>%
#   summarise(
#     RB = mean(bias),
#     .groups = "drop"
#   )
# 
# ## percent bias
# pb_df <- master_df %>%
#   filter(true_value != 0) %>%
#   group_by(
#     model, nsize, mech, mprop, miboot,
#     method, param_type, parameter
#   ) %>%
#   summarise(
#     PB = mean(perc_bias),
#     .groups = "drop"
#   )
# 
# ## root mean squared error
# rmse_df <- master_df %>%
#   group_by(
#     model, nsize, mech, mprop, miboot,
#     method, param_type, parameter
#   ) %>%
#   summarise(
#     RMSE = sqrt(mean(sq_error)),
#     .groups = "drop"
#   )
# 
# ## coverage rate
# cr_df <- master_df %>%
#   group_by(
#     model, nsize, mech, mprop, miboot,
#     method, param_type, parameter
#   ) %>%
#   summarise(
#     CR = mean(covered),
#     .groups = "drop"
#   )
# 
# ## average width
# aw_df <- master_df %>%
#   group_by(
#     model, nsize, mech, mprop, miboot,
#     method, param_type, parameter
#   ) %>%
#   summarise(
#     AW = mean(ci_width),
#     .groups = "drop"
#   )

## combine all
perf_df <- master_df %>%
  group_by(
    model, nsize, mech, mprop, miboot,
    method, param_type, parameter
  ) %>%
  summarise(
    RB   = mean(bias),
    PB   = mean(perc_bias),
    RMSE = sqrt(mean(sq_error)),
    CR   = mean(covered),
    AW   = mean(ci_width),
    .groups = "drop"
  )

# ## parameter aggregated version
# perf_paramtype <- perf_df %>%
#   group_by(
#     model, nsize, mech, mprop, miboot,
#     method, param_type, parameter
#   ) %>%
#   summarise(
#     RB   = mean(RB),
#     PB   = mean(PB),
#     RMSE = mean(RMSE),
#     CR   = mean(CR),
#     AW   = mean(AW),
#     .groups = "drop"
#   )

# glimpse(master_df)
save(master_df, perf_df, file = file.path(path_to_save, "master_results.RData"))

###############################################################################

load(file.path(path_to_save, "master_results.RData"))

###############################################################################

###############################################################################
# tables
###############################################################################

library(knitr)
library(kableExtra)

tab_pb <- perf_df %>%
  filter(
    model == "m1",
    nsize == 1000,
    mech  == "MAR",
    mprop == 0.1
  ) %>%
  select(
    miboot, method, parameter, PB
  ) %>%
  pivot_wider(
    names_from  = parameter,
    values_from = PB
  ) %>%
  arrange(miboot, method) %>%
  ungroup()
tab_pb <- tab_pb %>%
  mutate(
    method = recode(
      method,
      "knn_5"  = "kNN(5)",
      "knn_9"  = "kNN(9)",
      "knn_13" = "kNN(13)"
    )
  )
tab_pb <- tab_pb %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
kbl <- tab_pb %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    caption =
      "Percent bias (PB) by method and MI--bootstrap approach -- model 1 (Figure~1), $n = 1000$, MAR, missing proportion = 10\\%",
    align = c("l", "l", rep("r", ncol(tab_pb) - 2)),
    escape = FALSE
  ) %>%
  kable_styling(
    latex_options = c("repeat_header"),
    font_size = 7
  )
block_sizes <- tab_pb %>%
  count(miboot) %>%
  pull(n)
kbl <- kbl %>%
  pack_rows(
    index = setNames(
      block_sizes,
      unique(tab_pb$miboot)
    ),
    bold = FALSE
  )
kbl

###############################################################################
# graphs
###############################################################################

library(ggplot2)
library(patchwork)

make_perf_plot <- function(df, model_id, nsize_id, mech_id, mprop_id,
                           metric, title_tag = NULL, ylims = NULL) {
  plot_df <- df %>%
    filter(
      model == model_id,
      nsize == nsize_id,
      mech  == mech_id,
      mprop == mprop_id
    )
  
  boot_mi_methods <- c("norm", "pmm")
  mi_boot_methods <- c("norm", "pmm", "mean", "knn_5", "knn_9", "knn_13", "listwise")
  
  plot_df <- bind_rows(
    plot_df %>%
      filter(miboot == "bootmi", method %in% boot_mi_methods),
    plot_df %>%
      filter(miboot == "miboot", method %in% mi_boot_methods)
  )
  
  plot_df <- plot_df %>%
    mutate(
      method = recode(
        method,
        "knn_5"  = "kNN(5)",
        "knn_9"  = "kNN(9)",
        "knn_13" = "kNN(13)"
      ),
      method = factor(
        method,
        levels = c("norm", "pmm", "mean", "kNN(5)", "kNN(9)", "kNN(13)", "listwise")
      )
    )
  
  plot_df <- plot_df %>%
    mutate(
      vis_group = case_when(
        miboot == "bootmi" ~ "Boot MI",
        miboot == "miboot" & method %in% c("norm", "pmm") ~ "MI Boot",
        miboot == "miboot" ~ "Single imputation"
      )
    )

  plot_df <- plot_df %>%
    mutate(
      parameter_math = recode(
        parameter,
        gamma21  = "gamma[21]",
        gamma31  = "gamma[31]",
        gamma32  = "gamma[32]"
      )
    )
  
  y_label <- switch(
    metric,
    "PB"   = "Percent bias",
    "CR"   = "Coverage rate",
    "AW"   = "Average CI width",
    "RMSE" = "RMSE",
    metric
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = method,
      y = .data[[metric]],
      color = vis_group,
      shape = vis_group
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.55),
      size = 2.2
    ) +
    facet_wrap(
      ~ parameter_math,
      # scales = "free_y",
      labeller = label_parsed
    ) +
    scale_color_manual(
      values = c(
        "Boot MI"                = "#333333",
        "MI Boot"                = "#7F7F7F",
        "Single imputation"      = "#BFBFBF"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Boot MI"                = 16,  # solid circle
        "MI Boot"                = 17,  # solid triangle
        "Single imputation"      = 15   # solid square
      )
    ) +
    labs(
      x = "MI method",
      y = y_label,
      color = NULL,
      shape = NULL,
      title = title_tag
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",  # remove legend here
    )
  
  if (metric %in% c("RB", "PB", "CR")) {
    hval <- if (metric %in% c("RB", "PB")) 0 else 0.95
    p <- p +
      geom_hline(yintercept = hval, linetype = "dashed", linewidth = 0.4)
  }
  if (!is.null(ylims)) {
    p <- p + coord_cartesian(ylim = ylims)
  }
  
  p
}

get_common_ylims <- function(df, model_id, mech_id, metric) {
  vals <- df %>%
    filter(
      model == model_id,
      mech  == mech_id
    ) %>%
    pull(.data[[metric]])
  
  range(vals, na.rm = TRUE)
}

expand_range <- function(x, frac = 0.05) {
  rng <- range(x, na.rm = TRUE)
  d <- diff(rng)
  c(rng[1] - frac * d, rng[2] + frac * d)
}

perf_df_path <- perf_df %>%
  filter(param_type == "path")
metrics <- c("PB", "CR") #, "RB", "RMSE", "AW")
for (mtr in metrics) {
  for (mod in model_nm) {
    for (mch in mech_nm) {
      plts <- list()
      idx <- 1
      for (nsz in nsample) {
        for (mp in missprop) {
          ylims <- expand_range(
            get_common_ylims(perf_df_path, mod, mch, mtr)
          )
          plts[[idx]] <- make_perf_plot(perf_df_path, mod, nsz, mch, mp, mtr,
            title_tag = paste0("n = ", nsz, " — ", "missing prop = ", mp),
            ylims = ylims)
          idx <- idx + 1
        }
      }
      full_plot <- (plts[[1]] + plts[[2]]) / (plts[[3]] + plts[[4]])
      legend_plot <- plts[[1]] + theme(legend.position = "bottom")
      full_plot <- full_plot +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
      
      plot_name <- paste0("fig_", mtr, "_", mod, "_", mch, ".png")
      ggsave(filename = file.path(path_figure, plot_name),
        plot = full_plot, width = 12, height = 10,
        units = "in", dpi  = 300, bg = "white")
    }
  }
}
