# Package-wide global variables
.plssemMIEnv <- new.env()

#' @exportPattern "^[[:alpha:]]+"
#'
#' @importFrom utils globalVariables
#' @importFrom tools file_path_as_absolute
.onLoad <- function(lib, pkg){
  # quiets concerns of R CMD check re: the .'s that appear in pipelines
  # (see https://github.com/tidyverse/magrittr/issues/29)
  if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(".", "", "estimate", "param", "x"))
  }
  .plssemMIEnv$path.to.me <- tools::file_path_as_absolute(lib)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf("Package %s (%s) loaded.\nTo cite, type citation(\"%s\")",
    pkgname, utils::packageDescription(pkgname)$Version, pkgname))
  # packageStartupMessage("Caution: ", pkgname, " is under active developement! Breaking changes may occur in the future.")
    packageStartupMessage("Please report improvements and bugs to: https://github.com/sergioventurini/plssemMI/issues")
}
