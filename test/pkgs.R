options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c(
  "boot",
  "cSEM",
  "cSEM.DGP",
  "gamlss",
  "ggplot2",
  "graphics",
  "Hmisc",
  "ImputeRobust",
  "magrittr",
  "mice",
  "mvtnorm",
  "parallel",
  "SimDesign",
  "simstandard",
  "stats4",
  "tools",
  "tidyr",
  "tidyselect",
  "knitr",
  "testthat",
  "methods",
  "stats",
  "utils"
)

# install helper packages if needed
if (!requireNamespace("gh", quietly = TRUE)) {
  install.packages("gh")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

for (pkg in pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(pkg, " is already installed")
    next
  }

  message("Trying CRAN for ", pkg, " ...")

  cran_ok <- tryCatch({
    install.packages(pkg, dependencies = TRUE)
    requireNamespace(pkg, quietly = TRUE)
  }, error = function(e) {
    FALSE
  })

  if (cran_ok) {
    message(pkg, " installed from CRAN")
    next
  }

  message("Trying GitHub for ", pkg, " ...")

  gh_ok <- tryCatch({
    res <- gh::gh(
      "/search/repositories",
      q = paste0(pkg, " in:name language:R")
    )

    if (length(res$items) == 0) {
      message(pkg, " not found on GitHub")
      FALSE
    } else {
      repo_names <- vapply(res$items, function(x) x$name, character(1))
      full_names <- vapply(res$items, function(x) x$full_name, character(1))

      # prefer exact repo name match
      exact <- which(tolower(repo_names) == tolower(pkg))

      repo_to_install <- if (length(exact) >= 1) {
        full_names[exact[1]]
      } else if (length(full_names) == 1) {
        full_names[1]
      } else {
        message(
          pkg, ": multiple GitHub matches found: ",
          paste(full_names, collapse = ", ")
        )
        FALSE
      }

      if (isFALSE(repo_to_install)) {
        FALSE
      } else {
        remotes::install_github(repo_to_install, dependencies = TRUE)
        requireNamespace(pkg, quietly = TRUE)
      }
    }
  }, error = function(e) {
    message("GitHub install failed for ", pkg, ": ", e$message)
    FALSE
  })

  if (gh_ok) {
    message(pkg, " installed from GitHub")
  } else {
    message("Could not install ", pkg)
  }
}

###

pkgs <- c(
  "boot",
  "cSEM",
  "cSEM.DGP",
  "gamlss",
  "ggplot2",
  "graphics",
  "Hmisc",
  "ImputeRobust",
  "magrittr",
  "mice",
  "mvtnorm",
  "parallel",
  "SimDesign",
  "simstandard",
  "stats4",
  "tools",
  "tidyr",
  "tidyselect",
  "knitr",
  "testthat",
  "methods",
  "stats",
  "utils"
)
all(pkgs %in% rownames(installed.packages()))
