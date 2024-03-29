% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{poolLengthComp}
\alias{poolLengthComp}
\title{Utility that pools length data when multiple columns exist}
\usage{
poolLengthComp(LengthCompObj, byGroup = FALSE)
}
\arguments{
\item{LengthCompObj}{A life history object.}

\item{byGroup}{A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.}
}
\description{
Utility that pools length data when multiple columns exist
}
\examples{
poolLengthComp(fishLengthAssess::LengthCompExampleFreq)
}
