% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methodsAssess.R
\name{PmatFunc}
\alias{PmatFunc}
\title{Proportion of the catch that is above length at maturity}
\usage{
PmatFunc(LifeHistoryObj, LengthCompObj, byGroup = FALSE)
}
\arguments{
\item{LifeHistoryObj}{A LifeHistory object from fishSimGTG.}

\item{LengthCompObj}{A LengthComp object}

\item{byGroup}{A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.}
}
\description{
Froese, R., 2004. Keep it simple: three indicators to deal with overfishing. Fish. Fish 5, 86-91. Proportion of the length observations that are equal to or greater than the length at 50\% maturity parameter from a fishSimGTG LifeHistory object.
}
\examples{
library(fishSimGTG)
PmatFunc(fishSimGTG::LifeHistoryExample, fishLengthAssess::LengthCompExampleFreq, byGroup = FALSE)
}
