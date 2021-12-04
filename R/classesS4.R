
#----------------------
#Length composition object
#----------------------

#Roxygen header
#'Length composition object
#'
#'An S4 object that holds length composition data in a specified format.
#'
#'The length composition object is contains a data frame of length observations.
#'These observations are organized according to 'dataType'. When 'dataType' is set to 'lengths', each row contains
#'a length observation and columns represent discrete observation units, such as year.
#'When 'dataType' is set to 'freq', each row is the number of length observations corresponding to a length bin.
#'The first column must contain the mid-points of each length bin, with columns 2 and above representing discrete observation units
#' @param title A title for the object, useful for displaying the contents of the object
#' @param description A longer description of the object
#' @param L_units Units of measure for object, e.g., cm
#' @param L_type Options are TL, FL, or SL (Total length, fork length, and standard length, respectively)
#' @param L_source Options are FD or FI. FD is lengths or frequencies in the catch. FI is for fishery independent lengths or frequencies, such as a survey or research collection.
#' @param dataType Options are Length or Frequency (See details for additional information)
#' @param observationGroup Provide a description for interpreting column header, e.g., Year
#' @param header Logical. Whether dataframe contains a descriptive  header. FALSE will cause the header to be hidden
#' @param dt The length composition data frame.
#' @import LBSPR
#' @import fishSimGTG
#' @importFrom methods new

setClass("LengthComp",
         representation(
           title = "character",
           description = "character",
           L_units = "character",
           L_type = "character",
           L_source = "character",
           dataType = "character",
           observationGroup = "character",
           header = "logical",
           dt = "data.frame")
)

