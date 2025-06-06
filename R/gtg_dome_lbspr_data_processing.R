#' Process Length Composition Data
#'
#' @description Processes S4 LengthComp objects for use in GTG dome-shaped LBSPR model
#' @param LengthCompObj S4 LengthComp object containing length data
#' @param byGroup Logical, whether to process by group (default FALSE). If TRUE, 
#'   analyzes each group separately; if FALSE, pools all data together.
#' @param SizeBins List with Linc and ToSize elements (optional). If NULL, defaults 
#'   to 1cm bins up to reasonable maximum size.
#' @param Lc Length at first capture (default 0). For fishery-independent data, 
#'   must be >= 0. Fish smaller than Lc are removed from analysis.
#'   
#' @return List containing processed length data with elements:
#' \describe{
#'   \item{LenDat}{Frequency data matrix or vector}
#'   \item{LenMids}{Length bin midpoints}
#'   \item{group_names}{Names of groups analyzed}
#'   \item{n_groups}{Number of groups}
#'   \item{was_pooled}{Logical indicating if data was pooled}
#'   \item{original_dataType}{Original data format ("Frequency" or "Length")}
#'   \item{L_source}{Data source ("FI" or "FD")}
#'   \item{pooled_data}{Raw pooled data frame}
#'   \item{binWidth}{Width of length bins used}
#' }
#' @examples
#' # Example 1: Process frequency data
#' data(gtg_catch_frequency)
#' freq_obj <- new("LengthComp", 
#'                dt = gtg_catch_frequency, 
#'                dataType = "Frequency",
#'                L_source = "FD",
#'                header = TRUE)
#'
#' # Process pooled data
#' result1 <- processLengthCompData(freq_obj, byGroup = FALSE)
#' print(result1$group_names)
#' # Process by groups
#' result2 <- processLengthCompData(freq_obj, byGroup = TRUE)
#' print(result2$group_names)
#' # Example 2: Process raw length data
#' data(gtg_catch_lengths)
#' length_obj <- new("LengthComp",
#'                  dt = gtg_catch_lengths,
#'                  dataType = "Length", 
#'                  L_source = "FD",
#'                  header = TRUE)
#'
#' # Process pooled
#' result3 <- processLengthCompData(length_obj, byGroup = FALSE)
#' print(result3$group_names)
#' # Process by groups
#' result4 <- processLengthCompData(length_obj, byGroup = TRUE)
#' print(result4$group_names)
#' @export
processLengthCompData <- function(LengthCompObj, byGroup = FALSE, SizeBins = NULL, Lc = 0) {
  
# Input Validation (similar to lbsprWrapper validation)
# If any of these conditions are TRUE the fucntion stops and provude the message
  
  if(!is(LengthCompObj, "LengthComp") ||
     length(LengthCompObj@dt) == 0 ||
     length(LengthCompObj@dataType) != 1 ||
     !(LengthCompObj@dataType %in% c("Frequency", "Length")) ||
     !(LengthCompObj@L_source %in% c("FI", "FD")) ||
     (LengthCompObj@L_source == "FI" & Lc < 0)) {
    stop("Invalid LengthComp object or parameters")
  }
  
# Data processing: poolLengthComp(): it prepares the data, calls the fishLengthAssess package, handles the S4 object structure,
# applies byGroup logic (separate vs pooled),returns a clean data frame
  
# Frequency data - poolLengthComp() returned a data frame of pooled frequency data
  if(LengthCompObj@dataType == "Frequency") {
    dt <- poolLengthComp(LengthCompObj, byGroup)
    if(is.null(dt)) {
      stop("poolLengthComp returned NULL")
    }
    if(LengthCompObj@L_source == "FI") {
      dt <- dt[which(dt[,1] >= Lc),]    #Remove fish smaller than Lc
    }
  }
  
# Length data (raw measurements)- when byGroup=FALSE, LengthCompObj() removes grouping and returns a single sorted vector of all measured lengths
  if(LengthCompObj@dataType == "Length") {
    dt <- poolLengthComp(LengthCompObj, byGroup)
    if(is.null(dt)) {
      stop("poolLengthComp returned NULL")
    }
    
# Setting up the bins: setting binWidth (following lbspr wrapper). 
# The following code uses bin width from SizeBins, or defaults to 1 cm
# Basicaly, if SizeBins exists and  it has a Linc value, use that value as bin width. 
# Otherwise 1 cm  
    binWidth <- if(!is.null(SizeBins) && !is.null(SizeBins$Linc)) SizeBins$Linc else 1
    if(binWidth <= 0) binWidth <- 1
    
# Create bins similar to lbsprWrapper but using SizeBins if provided
# if SizeBins and SizeBins$ToSize both exist, it defines new bin edges using these.
    if(!is.null(SizeBins) && !is.null(SizeBins$ToSize)) {
      Lbins_new <- seq(from = 0, to = SizeBins$ToSize, by = binWidth)
    } else {
      
# If no SizeBins are not provided, falls back to an automatic method based on the data (like lbspr)
# figure out the range from the actual data
    max_length <- max(max(dt, na.rm = TRUE) + binWidth, 150)  # Use reasonable default - Use whichever is larger: the biggest fish + 1 bin, OR 150 cm
  
# define the bins using the provided data range    
    Lbins_new <- seq(from = min(dt, na.rm = TRUE), to = max_length, by = binWidth)
    }
    
# Convert raw lengths to frequencies (If raw lengths are provided, these are converted into frequency counts) 
# Following lpspr logic    
# Convert raw lengths to frequencies for each group/column
# For each column (e.g., Year_1, Year_2, etc. takes all the raw length measurements
# Uses hist() to count how many fish fall in each length bin
# Return the frequency counts for each bin
# tmp contains the frequency table
    tmp <- sapply(X = 1:NCOL(dt), function(X) {
      if(LengthCompObj@L_source == "FI") {
        dtIn <- dt[which(dt[,X] >= Lc), X]     #for FI data, remove small fish 
      } else {
        dtIn <- dt[,X]
      }
      dtIn <- dtIn[!is.na(dtIn)]
      hist(dtIn, Lbins_new, plot = FALSE)$counts # extract only the number of fish in each bin
    })
    
# Create midpoints (calculate midpoints of each bin) and combine
# Lbins_new: contains the edges of the length bins (e.g., 0,1,2,3,4,5).
# that defines 5 bins: 0-1, 1-2, 2-3, 3-5, 4-5   
    Lmids_new <- Lbins_new[1:(length(Lbins_new)-1)] + binWidth/2  #Lbins_new[1:5] + 1/2 = c(0.5, 1.5, 2.5, 3.5, 4.5)
    dt <- cbind(Lmids_new, tmp)
  }
  
# Group handling (after binning the data) - Following lpspr logic  
# Extract length data matrix
  if(byGroup && NCOL(dt) > 2) {

# Multiple groups - analyze each separately
    LenDat <- as.matrix(dt[,-1])  # Remove length column (keep data column)
    LenMids <- dt[,1]             # keep midpoints
    
# Get group names from the data
    if(LengthCompObj@header) {
      if(LengthCompObj@dataType == "Frequency") {
        group_names <- colnames(LengthCompObj@dt[,-1]) # Skip "Length" column
      } else {
        group_names <- colnames(LengthCompObj@dt)      # For length data - All columns are groups
      }
    } else {

# if not headers - create generic names
      if(LengthCompObj@dataType == "Frequency") {
        group_names <- paste0("Group", 1:NCOL(LengthCompObj@dt[,-1]))
      } else {
        group_names <- paste0("Group", 1:NCOL(LengthCompObj@dt))
      }
    }
    
  } else {

# Either only one group or byGroup = FALSE (pooled case)
# takes all the data columns and adds them together and creates one combined frequency distribution

# two cols (midpoints + one group)    
       if(NCOL(dt) == 2) {
      LenDat <- dt[,2]  # Single frequency column
      LenMids <- dt[,1]
      group_names <- "Pooled"
    } else {

# Multiple columns but pooled - sum across groups
      LenDat <- rowSums(dt[,-1], na.rm = TRUE)
      LenMids <- dt[,1]
      group_names <- "Pooled"
    }
  }
  
# Clean the data and create a clean data list (ensure all numeric)
  if(is.matrix(LenDat)) {
    LenDat[is.na(LenDat)] <- 0              # Replace NAs with zeros
    LenDat <- apply(LenDat, 2, as.numeric)  # Ensure numeric
    n_groups <- NCOL(LenDat)
  } else {
# when just one group or byGroup=FALSE (pooled)    
    LenDat <- as.numeric(LenDat)
    LenDat[is.na(LenDat)] <- 0
    n_groups <- 1
  }

#final output: the cleanded LendDat
  return(list(
    LenDat = LenDat,                  # The actual data used in the optimization
    LenMids = LenMids,                # Length bin midpoints
    group_names = group_names,
    n_groups = n_groups,
    was_pooled = !byGroup,            # Whether data was pooled (if byGroup = FALSE, so was_pooled = TRUE)
    original_dataType = LengthCompObj@dataType, # Original format
    L_source = LengthCompObj@L_source,
    pooled_data = dt,
    binWidth = if(LengthCompObj@dataType == "Length") binWidth else SizeBins$Linc
  ))
}
