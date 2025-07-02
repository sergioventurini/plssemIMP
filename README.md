# plssemIMP

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/plssemIMP)](https://cran.r-project.org/package=plssemIMP)
[![](http://cranlogs.r-pkg.org/badges/grand-total/plssemIMP?color=blue)](https://cran.r-project.org/package=plssemIMP)
[![Travis build status](https://travis-ci.org/sergioventurini/plssemIMP.svg?branch=master)](https://travis-ci.org/sergioventurini/plssemIMP)

<!-- badges: end -->

## Overview

###### Current release: 0.1.11
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
devtools::install_github("sergioventurini/plssemIMP")
```

See the main help page of the package by executing `?plssemIMP` or run any of
the demos available in the package by executing the code
`demo(demo_name, package = "plssemIMP")`.

## Authors
Sergio Venturini, Department of Economic and Social Sciences, Università Cattolica del Sacro Cuore, Cremona, Italy

E-mail: sergio.venturini@unicatt.it

Mehmet Mehmetoglu, Department of Psychology, NTNU, Norwegian University of Science and Technology, Trondheim, Norway

E-mail: mehmet.mehmetoglu@ntnu.no
