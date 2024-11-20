
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

The length composition object is a length data set that can be
structured either as raw data (i.e. individual length measurements) or
length frequency data (i.e. numbers of fish per length bin). For the
length frequency data, the first column of the data set must contain the
mid-points of the length bins, and the remaining columns contain the
counts for each length class.

``` r
# Load the package
library(fishLengthAssess)

# Example data set with length composition defined as length measurements by year
knitr::kable(head(LengthCompExampleFreq@dt))
```

| LMids | 2015 | 2016 | 2017 | 2018 | 2019 |
|------:|-----:|-----:|-----:|-----:|-----:|
|    11 |    0 |    0 |    0 |    0 |    0 |
|    13 |    0 |    0 |    0 |    0 |    0 |
|    15 |    0 |    0 |    0 |    0 |    0 |
|    17 |    0 |    0 |    0 |    0 |    0 |
|    19 |    0 |    0 |    0 |    0 |    0 |
|    21 |    0 |    0 |    0 |    0 |    0 |

``` r

# Example data set with length composition defined as length frequencies by year
knitr::kable(head(LengthCompExampleLength@dt))
```

| 2015 | 2016 | 2017 | 2018 | 2019 |
|-----:|-----:|-----:|-----:|-----:|
|   55 |   49 |   51 |   43 |   29 |
|   55 |   53 |   53 |   43 |   31 |
|   55 |   53 |   53 |   45 |   31 |
|   57 |   57 |   55 |   45 |   33 |
|   57 |   57 |   57 |   45 |   33 |
|   59 |   57 |   59 |   47 |   33 |

## Vignette

Documentation for the fishLengthAssess package can be found here as a
vignette:
