#' Example Catch Frequency Data for GTG dome-shaped LBSPR Analysis (5 groups)
#'
#' Simulated catch-at-length frequency data for 5 catch samples, generated using
#' the GTG dome-shaped LBSPR simulation with Normal.sca selectivity curve.
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
#' - FM = 1.0
#' - Normal.sca selectivity: SL1 = 55.28, SL2 = 4.343513^2
#' - Reference mesh: 13.5 cm, fishery mesh: 15.9 cm
#' - Sample size: 500 fish per sample (multinomial distribution)
#'
#' @source Generated using \code{\link{GTGDomeLBSPRSim2}} with example parameters
#' @seealso \code{\link{gtg_catch_lengths}}, \code{\link{gtg_catch_frequency_1group}},
#'   \code{\link{run_grouped_and_pooled}}
"gtg_catch_frequency"


#' Example Individual Length Data for GTG dome-shaped LBSPR Analysis (5 groups)
#'
#' Simulated individual fish length measurements for 5 catch samples, corresponding
#' to the frequency data in \code{\link{gtg_catch_frequency}}. Generated using the
#' GTG dome-shaped LBSPR simulation with Normal.sca selectivity curve.
#'
#' @format A data frame with individual fish (rows) and samples (columns):
#' \describe{
#'   \item{Catch_1}{Individual length measurements (cm) for sample 1}
#'   \item{Catch_2}{Individual length measurements (cm) for sample 2}
#'   \item{Catch_3}{Individual length measurements (cm) for sample 3}
#'   \item{Catch_4}{Individual length measurements (cm) for sample 4}
#'   \item{Catch_5}{Individual length measurements (cm) for sample 5}
#' }
#'
#' @details
#' Data simulated using the same parameters as \code{\link{gtg_catch_frequency}}:
#' - Linf = 120 cm, K = 0.2, M/K = 1.5
#' - L50 = 60 cm, L95 = 62 cm
#' - FM = 1.0
#' - Normal.sca selectivity: SL1 = 55.28, SL2 = 4.343513^2
#' - Reference mesh: 13.5 cm, fishery mesh: 15.9 cm
#' - Sample size: 500 fish per sample (multinomial distribution)
#'
#' @source Generated using \code{\link{GTGDomeLBSPRSim2}} with example parameters
#' @seealso \code{\link{gtg_catch_frequency}}, \code{\link{gtg_catch_lengths_1group}},
#'   \code{\link{run_grouped_and_pooled}}
"gtg_catch_lengths"


#' Example Catch Frequency Data for GTG dome-shaped LBSPR Analysis (1 group)
#'
#' Simulated catch-at-length frequency data for a single catch sample, generated
#' using the GTG dome-shaped LBSPR simulation with Normal.sca selectivity curve.
#' This dataset is the single-group version of \code{\link{gtg_catch_frequency}}.
#'
#' @format A data frame with length bins (rows) and one catch sample (column):
#' \describe{
#'   \item{Length}{Length bin midpoints (cm)}
#'   \item{Catch_1}{Catch frequencies for the single sample}
#' }
#'
#' @details
#' Data simulated using the same parameters as \code{\link{gtg_catch_frequency}}:
#' - Linf = 120 cm, K = 0.2, M/K = 1.5
#' - L50 = 60 cm, L95 = 62 cm
#' - FM = 1.0
#' - Normal.sca selectivity: SL1 = 55.28, SL2 = 4.343513^2
#' - Reference mesh: 13.5 cm, fishery mesh: 15.9 cm
#' - Sample size: 500 fish (multinomial distribution)
#'
#' @source Generated using \code{\link{GTGDomeLBSPRSim2}} with example parameters
#' @seealso \code{\link{gtg_catch_frequency}}, \code{\link{gtg_catch_lengths_1group}},
#'   \code{\link{run_grouped_and_pooled}}
"gtg_catch_frequency_1group"


#' Example Individual Length Data for GTG dome-shaped LBSPR Analysis (1 group)
#'
#' Simulated individual fish length measurements for a single catch sample,
#' generated using the GTG dome-shaped LBSPR simulation with Normal.sca
#' selectivity curve. This dataset is the single-group version of
#' \code{\link{gtg_catch_lengths}}.
#'
#' @format A data frame with individual fish (rows) and one sample (column):
#' \describe{
#'   \item{Catch_1}{Individual length measurements (cm) for the single sample}
#' }
#'
#' @details
#' Data simulated using the same parameters as \code{\link{gtg_catch_frequency}}:
#' - Linf = 120 cm, K = 0.2, M/K = 1.5
#' - L50 = 60 cm, L95 = 62 cm
#' - FM = 1.0
#' - Normal.sca selectivity: SL1 = 55.28, SL2 = 4.343513^2
#' - Reference mesh: 13.5 cm, fishery mesh: 15.9 cm
#' - Sample size: 500 fish (multinomial distribution)
#'
#' @source Generated using \code{\link{GTGDomeLBSPRSim2}} with example parameters
#' @seealso \code{\link{gtg_catch_lengths}}, \code{\link{gtg_catch_frequency_1group}},
#'   \code{\link{run_grouped_and_pooled}}
"gtg_catch_lengths_1group"
