# plssemMI

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/plssemMI)](https://cran.r-project.org/package=plssemMI)
[![](http://cranlogs.r-pkg.org/badges/grand-total/plssemMI?color=blue)](https://cran.r-project.org/package=plssemMI)
[![Travis build status](https://travis-ci.org/sergioventurini/plssemMI.svg?branch=master)](https://travis-ci.org/sergioventurini/plssemMI)

<!-- badges: end -->

## Overview

###### Current release: 0.1.0
###### R version required: at least 3.6.0
`R` package that implements different missing data imputation approaches
for partial least squares structural equation models (PLS-SEM).

## Installation

Since the package requires some code to be compiled, you need a working C++
compiler. To get it:

- On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
- On Mac, install Xcode from the app store.
- On Linux, `sudo apt-get install r-base-dev` or similar.

Then, the easiest way to get the package is to install it from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("sergioventurini/plssemMI")
```

See the main help page of the package by executing `?plssemMI` or run any of
the demos available in the package by executing the code
`demo(demo_name, package = "plssemMI")`.

## Authors
Sergio Venturini, Department of Economic and Social Sciences, Università Cattolica del Sacro Cuore, Cremona, Italy

E-mail: sergio.venturini@unicatt.it

Mehmet Mehmetoglu, Department of Psychology, NTNU, Norwegian University of Science and Technology, Trondheim, Norway

E-mail: mehmet.mehmetoglu@ntnu.no
