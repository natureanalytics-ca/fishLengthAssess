# This script generates example data for fishLengthAssess package, specifically for
# the GTG dome-shaped LBSPR analysis
# Run this script from within the package development environment
# Functions GTGDomeLBSPRSim2, new("LifeHistory"), etc. are available from the package

# Create LifeHistory S4 object using direct slot assignment
LifeHistoryObj <- new("LifeHistory")
LifeHistoryObj@title <- "Fish Stock"
LifeHistoryObj@speciesName <- "Generic Fish"
LifeHistoryObj@shortDescription <- "Stock for LBSPR Analysis"
LifeHistoryObj@L_type <- "TL"
LifeHistoryObj@L_units <- "cm"
LifeHistoryObj@Walpha_units <- "g"

# Core life history parameters
LifeHistoryObj@Linf <- 120
LifeHistoryObj@K <- 0.2  # Calculated from your MK ratio
LifeHistoryObj@t0 <- 0
LifeHistoryObj@L50 <- 60
LifeHistoryObj@L95delta <- 2  # L95 - L50 = 62 - 60
LifeHistoryObj@M <- 0.3  # Calculated from MK * K = 1.5 * 0.2
LifeHistoryObj@MK <- 1.5  # Your original M/K ratio

# Length-weight and reproduction parameters
LifeHistoryObj@LW_A <- 0.01
LifeHistoryObj@LW_B <- 3
LifeHistoryObj@Steep <- 0.7
LifeHistoryObj@R0 <- 1E6
LifeHistoryObj@Tmax <- -log(0.01) / LifeHistoryObj@M

# Additional parameters for spawning (set defaults)
LifeHistoryObj@recSD <- 0.6
LifeHistoryObj@recRho <- 0
LifeHistoryObj@isHermaph <- FALSE
LifeHistoryObj@H50 <- 0
LifeHistoryObj@H95delta <- 0

# GTG-specific parameters as attributes
attr(LifeHistoryObj, "CVLinf") <- 0.1
attr(LifeHistoryObj, "MaxSD") <- 2
attr(LifeHistoryObj, "NGTG") <- 13
attr(LifeHistoryObj, "Mpow") <- 0
attr(LifeHistoryObj, "FecB") <- 3

#define fleet parameters
FleetPars <- NULL
FleetPars$FM <- 1

# Normal.sca (normal with proportional spread)
FleetPars$selectivityCurve <- "Normal.sca"
FleetPars$SL1 <- 55.28228    # Modal length parameter
FleetPars$SL2 <- 4.343513^2  # Spread parameter (squared)
FleetPars$SLmesh <- c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19)
FleetPars$SLMin <- NA
FleetPars$use_aggregated <- TRUE

# FleetPars$selectivityCurve <- "Logistic"
# FleetPars$SL1 = 50
# FleetPars$SL2 = 60

SizeBins <- list()
SizeBins$Linc <- 1
SDLinf <- attr(LifeHistoryObj, "CVLinf") * LifeHistoryObj@Linf
SizeBins$ToSize <- LifeHistoryObj@Linf + SDLinf * attr(LifeHistoryObj, "MaxSD")

#Simulate the dynamic based on life history and fleet parameters
simGTG <- GTGDomeLBSPRSim2(LifeHistoryObj, FleetPars, SizeBins)

# Helper function to convert frequency data to individual lengths
convert_freq_to_individual_lengths <- function(lengths, frequencies) {
  unlist(mapply(rep, lengths, frequencies))
}

# Generate frequency data and individual lengths
simulate_frequency_and_lengths <- function(simGTG, n_samples = 5, sample_size = 500, seed = 123) {
  set.seed(seed)
  len_bins <- simGTG$LenMids
  sel_probs <- simGTG$LCatchFished / sum(simGTG$LCatchFished)
  
  # Remove any zero or negative probabilities
  valid_indices <- sel_probs > 0
  len_bins <- len_bins[valid_indices]
  sel_probs <- sel_probs[valid_indices]
  
  # Renormalize probabilities
  sel_probs <- sel_probs / sum(sel_probs)
  
  freq_df <- data.frame(Length = len_bins)
  ind_list <- list()
  
  for (i in 1:n_samples) {
    freq_vec <- as.numeric(rmultinom(1, sample_size, sel_probs))
    freq_df[[paste0("Catch_", i)]] <- freq_vec
    ind_list[[paste0("Year_", i)]] <- convert_freq_to_individual_lengths(len_bins, freq_vec)
  }
  
  # Convert individual lengths to wide format
  max_n <- max(sapply(ind_list, length))
  len_matrix <- data.frame(matrix(NA, nrow = max_n, ncol = n_samples))
  colnames(len_matrix) <- names(ind_list)
  for (i in seq_along(ind_list)) {
    len_matrix[1:length(ind_list[[i]]), i] <- ind_list[[i]]
  }
  
  return(list(catch_matrix = freq_df, length_matrix = len_matrix))
}

# Run simulation
cat("Generating simulated catch data...\n")
simulated_data <- simulate_frequency_and_lengths(simGTG, n_samples = 5, sample_size = 500)

# Create the final 2 datasets for the package
gtg_catch_frequency <- simulated_data$catch_matrix      # Frequency data
gtg_catch_lengths <- simulated_data$length_matrix    # Individual lengths

# Check the data
cat("Frequency data (first few rows):\n")
print(head(gtg_catch_frequency))

cat("\nLength data (first few rows):\n")
print(head(gtg_catch_lengths))

cat("\nData ranges:\n")
cat("Frequency data length range:", range(gtg_catch_frequency$Length), "\n")
cat("Individual length data range:", range(unlist(gtg_catch_lengths), na.rm = TRUE), "\n")

# Save as .rda files in the data/ directory instead of CSV
usethis::use_data(gtg_catch_frequency, overwrite = TRUE)
usethis::use_data(gtg_catch_lengths, overwrite = TRUE)