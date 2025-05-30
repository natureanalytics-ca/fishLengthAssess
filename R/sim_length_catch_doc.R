#' Example Catch Frequency Data for GTG dome-shaped LBSPR Analysis
#'
#' Simulated catch-at-length frequency data for 5 samples, generated using
#' GTG dome selectivity simulation with Normal.sca selectivity curve.
#'
#' @format A data frame with length bins (rows) and catch samples (columns):
#' \describe{
#'   \item{Length}{Length bin midpoints (cm)}
#'   \item{Catch_1}{Catch frequencies for sample 1}
#'   \item{Catch_2}{Catch frequencies for sample 2}
#'   \item{Catch_3}{Catch frequencies for sample 3}
#'   \item{Catch_4}{Catch frequencies for sample 4}
#'   \item{Catch_5}{Catch frequencies for sample 5}
#' }
#'
#' @details
#' Data simulated using the following parameters:
#' - Linf = 120 cm, K = 0.2, M/K = 1.5
#' - L50 = 60 cm, L95 = 62 cm  
#' - Normal.sca selectivity: SL1 = 55.28, SL2 = 18.87
#' - Mesh sizes: 13.5-19 cm
#' - Sample size: 500 fish per sample (using multinomial distribution)
#'
#' @source Generated using \code{\link{GTGDomeLBSPRSim2}} with example parameters
#' @seealso \code{\link{gtg_catch_lengths}}, \code{\link{DoOptDome}}, \code{\link{DoOptDome.LengthComp}}
"gtg_catch_frequency"

#' Example Individual Length Data for GTG Analysis
#'
#' Simulated individual fish length measurements for 5 samples, corresponding
#' to the frequency data in gtg_catch_frequency. Generated using GTG dome 
#' selectivity simulation with Normal.sca selectivity curve.
#'
#' @format A data frame with individual fish (rows) and samples (columns):
#' \describe{
#'   \item{Year_1}{Individual length measurements (cm) for sample 1}
#'   \item{Year_2}{Individual length measurements (cm) for sample 2}
#'   \item{Year_3}{Individual length measurements (cm) for sample 3}
#'   \item{Year_4}{Individual length measurements (cm) for sample 4}
#'   \item{Year_5}{Individual length measurements (cm) for sample 5}
#' }
#' Data represents realistic catch from a Normal.sca dome-shaped selectivity fishery.
#'
#' @details
#' Data simulated using the same parameters as gtg_catch_frequency:
#' - Linf = 120 cm, K = 0.2, M/K = 1.5
#' - L50 = 60 cm, L95 = 62 cm  
#' - Normal.sca selectivity: SL1 = 55.28, SL2 = 18.87
#' - Mesh sizes: 13.5-19 cm
#' - Sample size: 500 fish per sample (using multinomial distribution)
#'
#' @source Generated using \code{\link{GTGDomeLBSPRSim2}} with example parameters
#' @seealso \code{\link{gtg_catch_frequency}}, \code{\link{GTGDomeLBSPRSim2}}
"gtg_catch_lengths"