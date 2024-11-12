
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

The main functions in the package rely on life histories and length
composition objects. The life history object gathers information about
life history traits of each species (length-weight relationship, von
Bertalanffy growth parameters, length at 50% maturity, natural
mortality, and others). The life history object can be created in the R
package fishSimGTG, see more details
[here](https://natureanalytics-ca.github.io/fishSimGTG/om-pop.html#life-history-object).

``` r
library(fishLengthAssess)
## basic example code
```
