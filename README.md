
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biPOD <a href="caravagnalab.github.io/biPOD"><img src="man/figures/logo.png" align="right" height="69" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/caravagnalab/biPOD/workflows/R-CMD-check/badge.svg)](https://github.com/caravagnalab/biPOD/actions)
[![pkgdown](https://github.com/caravagnalab/biPOD/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/caravagnalab/biPOD/actions/workflows/pkgdown.yaml)
[![R-CMD-check](https://github.com/caravagnalab/biPOD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/caravagnalab/biPOD/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/caravagnalab/biPOD/branch/main/graph/badge.svg)](https://app.codecov.io/gh/caravagnalab/biPOD?branch=main)
<!-- badges: end -->

biPOD is a package to infer kinetic parameters of an evolving population
whose size is observed at discrete intervals of time. The tool is able
to:

1.  Infer growth rates and instant of birth of a population of cells
    when time windows and change points are known. This includes:
    - Calculating exponential growth rates for different time periods
      (e.g., treatment phase, relapse phase)
    - Determining the upper limit for the population’s birth time if
      growth is observed in the initial window
    - Handling multiple time windows with predefined change points
2.  Identify optimal change points in population dynamics when these
    transition times are not known, by:
    - Analyzing population size data to detect significant changes in
      growth patterns
    - Computing the most likely moments when population behavior shifts
    - Optimizing the placement of change points to best fit the observed
      data
3.  Analyze dual-population scenarios typical in treatment-to-relapse
    cases, characterized by:
    - Inferring separate growth rates for resistant and sensitive cell
      populations
    - Determining the birth time of the resistant population
    - Calculating the death time of the sensitive population
    - Modeling the characteristic U-shaped dynamics that emerge from the
      interplay between these populations

#### Help and support

## [![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/biPOD/-yellow.svg)](https://caravagnalab.github.io/biPOD)

## Installation

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/biPOD")
```

------------------------------------------------------------------------

#### Copyright and contacts

Cancer Data Science (CDS) Laboratory, University of Trieste, Italy.

[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
