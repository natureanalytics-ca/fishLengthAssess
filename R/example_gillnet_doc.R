
#' Example gillnet catch data
#'
#' This dataset contains example gillnet catch data across multiple mesh sizes (eight mesh sizes).
#' It is adapted from SELECT and TropFishR and used to demonstrate dome-shaped selectivity fitting.
#' First column: MidLength.
#' Second - ninth columns: Mesh sizes.
#'
#' @format A data frame with 11 rows and 9 variables:
#' \describe{
#'   \item{MidLength}{Midpoint of each length bin (e.g., cm)}
#'   \item{Mesh1–Mesh8}{Catch counts for each mesh size}
#' }
#'
#' @source Adapted from SELECT and TropFishR
"raw_data_gillnet"