
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fishLengthAssess

<!-- badges: start -->
<!-- badges: end -->

This repository contains the code and documentation for the
fishLengthAssess R package. fishLengthAssess acts as an aggregator of
length-based stock assessment models. The package imposes a standardized
formatting and metadata structure on length data sets, which can then be
used in a variety of length-based assessment methods that are aggregated
within the package. It also provides wrapper functions to link
standardized formatting of length data sets with external R packages
that contain length-based assessment methods (e.g. LBSPR package).

## Installation

You can install the development version of fishLengthAssess from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("natureanalytics-ca/fishLengthAssess")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fishLengthAssess)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
