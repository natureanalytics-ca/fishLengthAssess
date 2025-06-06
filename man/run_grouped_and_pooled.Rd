% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtg_dome_lbspr.R
\name{run_grouped_and_pooled}
\alias{run_grouped_and_pooled}
\title{Run Both Grouped and Pooled Optimization}
\usage{
run_grouped_and_pooled(
  lifeHistoryObj,
  fixedFleetPars,
  LengthCompObj,
  SizeBins = NULL,
  Lc = 0,
  mod = "GTG"
)
}
\arguments{
\item{lifeHistoryObj}{S4 life history object (same format as GTGDomeLBSPRSim2)}

\item{fixedFleetPars}{List of fixed fleet parameters including selectivity type}

\item{LengthCompObj}{S4 LengthComp object containing length composition data}

\item{SizeBins}{Optional list of size bin specifications with Linc and ToSize elements}

\item{Lc}{Length at first capture for removing small fish from analysis (default 0).
Only applies to fishery-independent data}

\item{mod}{Character string specifying model type (default "GTG").
Currently only "GTG" is supported}
}
\value{
List containing:
\itemize{
\item grouped: Results from DoOptDome.LengthComp() with byGroup = TRUE
\item pooled: Results from DoOptDome.aggregated()
}
}
\description{
Convenience function that runs both grouped and pooled optimization analyses
on the same dataset, allowing comparison between approaches.
}
\details{
This function simplifies comparative analysis by automatically running:
\enumerate{
\item Grouped analysis (each group fitted separately)
\item Pooled analysis (all groups combined)
}

This allows users to compare whether pooling data affects parameter estimates
and model fit.
}
\examples{
\donttest{
# Compare grouped vs pooled optimization results
# Create LifeHistory object
LifeHistoryObj <- new("LifeHistory")
LifeHistoryObj@Linf <- 120
LifeHistoryObj@K <- 0.2
LifeHistoryObj@L50 <- 60
LifeHistoryObj@L95delta <- 2
LifeHistoryObj@MK <- 1.5
LifeHistoryObj@LW_A <- 0.01
LifeHistoryObj@LW_B <- 3
LifeHistoryObj@Steep <- 0.7
LifeHistoryObj@R0 <- 1E6
attr(LifeHistoryObj, "NGTG") <- 13
attr(LifeHistoryObj, "CVLinf") <- 0.1
attr(LifeHistoryObj, "MaxSD") <- 2

# Create LengthComp object
data(gtg_catch_frequency)
LengthCompObj <- new("LengthComp",
                     dt = gtg_catch_frequency,
                     dataType = "Frequency",
                     L_source = "FD",
                     header = TRUE)

# Fleet parameters
fixedFleetPars <- list(
  selectivityCurve = "Normal.sca",
  SL1 = 55,
  SL2 = 18,
  SLmesh = c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19),
  use_aggregated = TRUE
)

# Run both analyses
results <- run_grouped_and_pooled(LifeHistoryObj, fixedFleetPars, LengthCompObj, mod="GTG")
 results$grouped$group_results$Catch_1$lbPars
results$pooled$lbPars
}
}
\seealso{
\code{\link{DoOptDome.LengthComp}}, \code{\link{DoOptDome.aggregated}}
}
