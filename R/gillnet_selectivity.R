#' @importFrom TropFishR select_Millar
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap labs theme_minimal geom_line geom_point scale_color_manual geom_vline theme guides guide_legend scale_size scale_fill_gradient2 geom_hline element_text xlim ggsave scale_linetype_manual scale_size_manual
#' @importFrom reshape2 melt
#' @importFrom pracma findpeaks
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr filter select arrange desc mutate across where %>%
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @importFrom stats qlogis quantile sd setNames
#' @importFrom grDevices rainbow
#' @importFrom utils write.csv
NULL

# Suppress R CMD check notes about ggplot2/dplyr variables
utils::globalVariables(c("MidLength", "Catch", "Mesh", "TotalCatch", "Length", 
                         "Selectivity", "MeshSize", "Residuals", "Curve", "Model", 
                         "LogLikelihood", "Deviance", "Mode1", "StdDev1", "Mode2", 
                         "StdDev2", "P_Mode1"))

#' Fit Gillnet Dome-Shaped Selectivity Models
#'
#' This function fits multiple selectivity models to gillnet catch data, including
#' both unimodal and bimodal selectivity curves. It automatically estimates starting
#' values, fits multiple models, and provides diagnostic plots.
#'
#' @param input_data A data frame with the first column as mid-length values and
#'   subsequent columns as catch data for each mesh size
#' @param mesh_sizes Numeric vector of mesh sizes corresponding to the catch data columns
#' @param run_bimodal Logical. Whether to include bimodal selectivity models in fitting.
#'   Default is TRUE
#' @param manual_x0_list Optional list of manual starting values for model parameters.
#'   Should be named list with elements corresponding to model types
#' @param length_seq Numeric sequence of lengths for plotting selectivity curves.
#'   Default is seq(40, 100, 0.1)
#' @param output_dir Character string specifying directory for saving plots.
#'   Default is "model_plots"
#' @param criterion Character string specifying model selection criterion.
#'   Default is "Deviance"
#' @param sd_spread Numeric value for spread around modes when detecting peaks.
#'   Default is 7
#' @param rel.power Optional numeric vector of relative fishing power for each mesh size.
#'   Default is NULL (equal power)
#' @param verbose Logical. Whether to print progress messages. Default is TRUE
#'
#' @return A list containing:
#' \itemize{
#'   \item results: Model fitting results from TropFishR::select_Millar
#'   \item summary_table: Data frame comparing all fitted models
#'   \item selectivity_curves: List of selectivity curves for each model
#' }
#' 
#' @details This function implements the following workflow:
#' \enumerate{
#'   \item Exploratory data analysis and automatic detection of modes
#'   \item Calculation of starting values for model parameters
#'   \item Fitting of multiple selectivity models (norm.loc, norm.sca, lognorm,
#'         and optionally binorm.sca, bilognorm)
#'   \item Generation of diagnostic plots and model comparison tables
#' }
#' @seealso \code{\link{compare_stats}}, \code{\link{plot_mesh_curves}}, \code{\link{get_composite_curve}}
#' @examples
#' data(raw_data_gillnet)
#' mesh_sizes <- c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19)
#' result <- fit_gillnet_dome(
#'   input_data = raw_data_gillnet,
#'   mesh_sizes = mesh_sizes,
#'   output_dir = tempdir(),
#'   length_seq = seq(40, 100, 1)
#' )
#' names(result$results)
#' 
#' @export
fit_gillnet_dome <- function(input_data,
                        mesh_sizes,
                        run_bimodal=TRUE,#new argument
                        manual_x0_list = list(),
                        length_seq = seq(40, 100, 0.1),
                        output_dir = "model_plots",
                        criterion = "Deviance",
                        sd_spread = 7,
                        rel.power = NULL, #added rel.power argument
                        verbose = TRUE) {
  
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  


  
  # This first line of code is used to tell the function to use manual staring values if available
  # or calculated staring values when manual are not provided 
  # this is useful for bimodal models only
  # basically it says if the first thing is not null, use it. Otherwise use the second one.
  # `%||%`  the 'backticks' are needed when creating the shortcut - once it's defined it works as normal operator
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b #Use a if it is not NULL, otherwise use b (allows to set default values dynamically, only if nothing is provided)
  
  #=================================================================#
  # Step 1: Exploratory plots and calculation of starting values X0 #
  #=================================================================#
  if (verbose) cat("\n[Step 1] Exploratory plots and auto x0 estimation...\n") # just print a message indicating that step 1 is beginning
  
  # Create a copy of the data for plotting
  
  plot_data <- input_data
  colnames(plot_data) <- c("MidLength", paste0("Mesh_", seq_along(mesh_sizes)))
  input_long <- melt(plot_data, id.vars = "MidLength", variable.name = "Mesh", value.name = "Catch")
  input_long$Mesh <- factor(rep(mesh_sizes, each = nrow(plot_data)))
  
  # Exploratory plots
  # Histogram per mesh size
  p1 <- ggplot(input_long, aes(x = MidLength, y = Catch, fill = Mesh)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    facet_wrap(~ Mesh, scales = "free_y") +
    labs(title = "Length Distribution by Mesh Size", x = "Fish Length (cm)", y = "Catch Count", fill = "Mesh Size (cm)") +
    theme_minimal()
  
  # Line plot across meshes
  p2 <- ggplot(input_long, aes(x = MidLength, y = Catch, color = Mesh, group = Mesh)) +
    geom_line(linewidth = 1.2) +
    geom_point() +
    scale_color_manual(values = rainbow(length(mesh_sizes))) +
    labs(title = "Catch Distribution Across Length Classes", x = "Fish Length (cm)", y = "Catch Count", color = "Mesh Size (cm)") +
    theme_minimal()
  
  # adds a new col that sums the catch across all mesh sizes for each length class.
  plot_data$TotalCatch <- rowSums(plot_data[, -1]) 
  
  # Bar plot of total catch across length classes (used to identify peaks)
  p3 <- ggplot(plot_data, aes(x = MidLength, y = TotalCatch)) +
    geom_bar(stat = "identity", fill = "grey", alpha = 0.7) +
    labs(title = "Aggregated Length Distribution", x = "Fish Length (cm)", y = "Total Catch Count") +
    theme_minimal()
  
  total_counts <- plot_data$TotalCatch
  
  # Finding peaks (mode) for bimodal models
  
  # Use findpeaks() to identify peaks in the vector total_counts, which represents total catch per length class 
  # findpeaks() finds local maxima in a numeric vector (from pracma package)
  # findpeaks() will return a matrix, where the number of rows is the number of peaks
  peaks <- findpeaks(total_counts, nups = 2, ndowns = 2, minpeakheight = max(total_counts) * 0.1) #nups=require at least 2 increasing points before the peak;ndowns= require at least 2 decreasing points after the peak;minpeakheight= ignore small peaks (less than 10% of max) 
  
  # If peaks is not NULL, and it has at least 2 rows, means  2 or more clear peaks in total catch across length classes
  if (!is.null(peaks) && nrow(peaks) >= 2) {
  
  # get the lengths where the peaks were found
  # the output is a vector of positions that correspond to the peaks
  # If two peaks are found, use them as Mode1 and Mode2  
    peak_lengths <- plot_data$MidLength[peaks[, 2]]#column 2 = index of peak (position)
    Mode1 <- peak_lengths[1]
    Mode2 <- peak_lengths[2]
    
  # If only one peak is found, estimate Mode2 based on other high-catch areas.
  } else if (!is.null(peaks) && nrow(peaks) == 1) {
    Mode1 <- plot_data$MidLength[peaks[1, 2]] #takes only first peak
  
  # It filters out lengths that are too close to Mode1 (within 5 cm) -Keeps only those lengths that are at least 5 cm away from Mode1 on either side.  
    remaining_lengths <- plot_data$MidLength[abs(plot_data$MidLength - Mode1) > 5] 
  
  # find the next best candidate peak  
    Mode2 <- if (length(remaining_lengths) > 0) {
    remaining_lengths[which.max(total_counts[abs(plot_data$MidLength - Mode1) > 5])] #If there are any valid remaining lengths, finds the length with the highest catch.
  
  # If we could not find any good candidates far from Mode1
  # we just use the 75th percentile length instead - to ensure Mode2 is on the higher side    
    }else {
      quantile(plot_data$MidLength, 0.75) #If no valid lengths exist (e.g., very short length range), falls back to using the 75th percentile (pick a high-end value based on percentiles) of all length classes as Mode2.   
    }
    if (verbose) message("Only one peak found. Mode2 estimated dynamically.")
  } 
  
  # If no clear peaks, fall back to using 30th and 70th percentiles (This ensures some spread across the length range)
    else {
    Mode1 <- quantile(plot_data$MidLength, 0.3)
    Mode2 <- quantile(plot_data$MidLength, 0.7)
    if (verbose) message("No clear peaks found. Using quantiles.")
    }
  
  # Ensure Mode2 is not too close to Mode1 (for stability in model fitting).
  # I need to improve this, I could use diff. values depending on the length range
  if (abs(Mode2 - Mode1) < 8) Mode2 <- Mode1 + 8 
  
  #It keeps all lengths within a window of ±sd_spread cm around Mode1 and Mode2. 
  #e.g. if Mode1 = 40 and sd_spread = 7, it keeps lengths from 33 to 47 cm
  subset1 <- plot_data$MidLength[plot_data$MidLength > (Mode1 - sd_spread) & plot_data$MidLength < (Mode1 + sd_spread)]#selects a subset of fish lengths around Mode1
  subset2 <- plot_data$MidLength[plot_data$MidLength > (Mode2 - sd_spread) & plot_data$MidLength < (Mode2 + sd_spread)]#selects a subset of fish lengths around Mode2
  
  #Now if If the SD is missing (NA), or too small (less than 3 cm) (too narrow curve). Then it defaults to 3.5 cm as a safe minimum value 
  # Also, I need to improve this. I could use diff. values
  StdDev1 <- ifelse(is.na(sd(subset1)) || sd(subset1) < 3, 3.5, sd(subset1))
  StdDev2 <- ifelse(is.na(sd(subset2)) || sd(subset2) < 3, 4.5, sd(subset2))
  
  # Convert raw SD to log-scale SD for the bilognorm model using lognormal theory
  # This is the correct mathematical conversion from raw-space standard deviation to lognormal scale standard deviation
  # log-scale standard deviations (for the bilognorm model) based on calc. Mode and Standard Deviation values.
  # The bilognorm model needs its inputs on a log scale
  log_sd_from_raw <- function(mean_val, sd_val) sqrt(log(1 + (sd_val / mean_val)^2))#standard deviation of a lognormal distribution
  LogStdDev1 <- min(max(log_sd_from_raw(Mode1, StdDev1), 0.1), 0.4) #Calculates the log-scale SD around Mode1 (bounded to 0.1 (avoid too narrow dome) to 0.4 (avoid too flat) - It can be modified)
  LogStdDev2 <- min(max(log_sd_from_raw(Mode2, StdDev2), 0.1), 0.4) #max: ensures the value is at least 0.1; min: ensures the value is at most 0.4
  
  #sum catch under each mode to calc.later the proportion of catch in Mode 1 (Pmode1)
  Catch_Mode1 <- sum(total_counts[plot_data$MidLength > (Mode1 - sd_spread) & plot_data$MidLength < (Mode1 + sd_spread)], na.rm = TRUE)# sum the catch around mode1, from Mode1 - sd_spread to Mode1 + sd_spread (total catch under the first mode)
  Catch_Mode2 <- sum(total_counts[plot_data$MidLength > (Mode2 - sd_spread) & plot_data$MidLength < (Mode2 + sd_spread)], na.rm = TRUE)# value >x and <y (the range around mode)
  
  #If no peaks were detected (peaks is null or fewer than 2), it just sets P_Mode1 = 0.95 by default
  P_Mode1 <- if (!is.null(peaks) && nrow(peaks) >= 2) min(max(Catch_Mode1 / (Catch_Mode1 + Catch_Mode2), 0.75), 0.95) else 0.95 # catch under M1/(catch under M1+M2) (min 0.75, and max 0.95 to avoid unrealistic value). Not peak detected, set to 0.95
  P_Mode1_logit <- qlogis(P_Mode1)#convert to logit scale
  
  #Starting values that will be passed to select_Millar()
  x0_binorm_sca <- c(Mode1, StdDev1, Mode2, StdDev2, P_Mode1_logit)
  x0_bilognorm <- c(log(Mode1), LogStdDev1, log(Mode2), LogStdDev2, P_Mode1_logit)
  
  # barplot with detected or calculated peaks (Mode1 and Mode2)
  p4 <- ggplot(plot_data, aes(x = MidLength, y = TotalCatch)) +
    geom_bar(stat = "identity", fill = "grey", alpha = 0.7) +
    geom_vline(xintercept = Mode1, color = "red", linetype = "dashed", size = 1.2) +
    geom_vline(xintercept = Mode2, color = "blue", linetype = "dashed", size = 1.2) +
    labs(title = "Detected or calculated Peaks (Mode1 & Mode2)", x = "Fish Length (cm)", y = "Total Catch Count") +
    theme_minimal()
  
  # just printing for visualization purposes, but these are also stored in a folder
    print(p1)
    print(p2)
    print(p3)
    print(p4)
  
  # combine the plots and save them
  combined_exploratory_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  exploratory_plot_path <- file.path(output_dir, "exploratory_plots_peaks.jpeg")
  
  ggsave(filename = exploratory_plot_path,plot = combined_exploratory_plot,
    width = 12,
    height = 10,
    dpi = 300,
    device = "jpeg")
  
  if (verbose) cat("Saved combined exploratory plots to:", exploratory_plot_path, "\n")
  
  #=================================================================#
  #                 Step 2: Now -Models fitting                     #
  #=================================================================#
  
  if (verbose) cat("\n[Step 2] Fitting selectivity models...\n")
  
  # just reset the original input data
  input_data <- input_data[, 1:(ncol(input_data))]  
  colnames(input_data) <- NULL
  
  # used when manual stating values are provided
  full_x0_list <- manual_x0_list 
  
  # here I am using the first line of this code to provide option of manual calc or automatic calc of starting values for bimodal models.
  full_x0_list$binorm.sca <- full_x0_list$binorm.sca %||% x0_binorm_sca #use the manual starting values provided (full_x0_list) or use the calculated starting values (x0_binorm_sca)
  full_x0_list$bilognorm <- full_x0_list$bilognorm %||% x0_bilognorm    #same for the other model
  
  # Default rel.power to 1 if NULL, otherwise define a vector length= Nmesh sizes
  if (is.null(rel.power)) rel.power <- rep(1, length(mesh_sizes))
  
  # prepare the data to be used by the selectivity models
  midLengths <- input_data[, 1] #extract mid lengths
  
  # extracts everything except the first column, convert into a matrix
  CatchPerNet_mat <- as.matrix(input_data[, -1]) 
  colnames(CatchPerNet_mat) <- NULL
  
  # create the data list
  data_list <- list(midLengths = midLengths, meshSizes = mesh_sizes, CatchPerNet_mat = CatchPerNet_mat, rel.power = rel.power)#need to be stored in a data list
  
  # Always include unimodal models
  models <- c("norm.loc", "norm.sca", "lognorm")
  
  # Conditionally include bimodal models if requested, the next line tests whether the value of run_bimodal is TRUE. If so, include bimodal models
  if (run_bimodal==TRUE) {
    manual_x0_list <- manual_x0_list %||% list()    # use starting values (the 2 options)
    models <- c(models, "binorm.sca", "bilognorm")  # and add the bimodal models
  }
  results <- list()# empty list to store each model outputs (one result oer model)
  
  # Fitting all the sel models (loop over every model)
  all_selectivity_curves <- list() # creates an empty list that will store all the selectivity curves for each model
  for (model in models) { #loops through each of the model types
    cat("\nFitting model:", model, "...\n") #prints a status message to show which model is being fitted
  
  # The following line does: 
  # Checking if the current model name (e.g., "norm.loc", "binorm.sca", etc.) is one of 
  # the names in the full_x0_list. 
  # If it is, the code retrieves the starting values for that model from full_x0_list. 
  # If not, it sets x0 to NULL  
  x0 <- if (model %in% names(full_x0_list)) full_x0_list[[model]] else NULL
  
  # using tryCatch to catch any error during the fitting
  # If an error occurs, the error = function(e) {...} part is executed
  # It prints an error message, returns NULL, which means this model will be marked as failed.
  # The loop then continues with the next model.
  # If the model fitting succeeds, the result is stored in results[[model]]
  
  results[[model]] <- tryCatch({
    
  # the next line does the model fitting  
      select_Millar(data_list, x0 = x0, rtype = model, rel.power = rel.power, plot = FALSE)
    }, error = function(e) {
      cat("Error fitting model:", model, "- Skipping.\n")
      return(NULL)
    })
  }
  
  # cleaning up the results - removes failed models from the result list (return True for returned NUll (failed models), false (fitted models)
  results <- results[!sapply(results, is.null)] # which one is not Null (reversed !)
  if (length(results) == 0) {
    cat("No models fitted successfully.\n")
    return(NULL)
  }
  # for example:
  # norm.loc    norm.sca    lognorm    binorm.sca    bilognorm
  # TRUE         FALSE       FALSE       TRUE         FALSE
  # ! to flip and keep only False entries
  
  # Summary: Extract key statistics- extract the model name from the result list
  # x is just a placeholder that represents each element of results
  # e.g., in sapply(results, function(x) x$out["model.l", 1]): sapply goes through each element of results (each fitted model)
  # For each model, it extracts x$out["model.l", 1] (the log-likelihood value)
  # It collects all these log-likelihood values into a vector
  summary_table <- data.frame(
    Model = names(results),
    LogLikelihood = sapply(results, function(x) x$out["model.l", 1]),
    Deviance = sapply(results, function(x) x$out["Deviance", 1]),
    Mode1 = sapply(results, function(x) x$estimates[1, "par"]),
    StdDev1 = sapply(results, function(x) x$estimates[2, "par"]),
    Mode2 = sapply(results, function(x) ifelse(nrow(x$estimates) > 2, x$estimates[3, "par"], NA)),
    StdDev2 = sapply(results, function(x) ifelse(nrow(x$estimates) > 3, x$estimates[4, "par"], NA)),
    P_Mode1 = sapply(results, function(x) ifelse(nrow(x$estimates) > 4, x$estimates[5, "par"], NA))
  )
  summary_table <- summary_table[order(summary_table$LogLikelihood, decreasing = TRUE), ]
  print(format(summary_table, digits = 3, nsmall = 3))
  
  #=================================================================#
  #               Step 3: Plot selectivity and residuals            #
  #=================================================================#
  
  # Load internal (non-exported) TropFishR fucntion
  rtypes_Millar <- get("rtypes_Millar", envir = getNamespace("TropFishR"))
  #if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (model in names(results)) {
    res <- results[[model]]
    meshSizes <- res$meshSizes
    nmeshes <- length(meshSizes)
    plotlens <- length_seq
    
  # In the next line, outer: applies the selectivity function for:
  # all combinations of length (plotlens) × mesh (meshSizes);
  # res$rtype: tells which model to use
  # res$par: take the estimated pars for that model
  # outer is an R function outer(rows, cols, fun,additional args). In this case outer()calculates the selectivity value for each combination of fish length (plotlens) and mesh size (meshSizes), using the appropriate selectivity function (rtypes_Millar(res$rtype)) and the estimated parameters(res$par).
    rmatrix <- outer(plotlens, meshSizes, rtypes_Millar(res$rtype), res$par)#matrix output: rows:lengths and cols:mesh sizes
    rmatrix <- t(t(rmatrix) * res$rel.power) # transpose to multiply each row by one number (mesh power), then t() again to bring matrix back to original shape
    
  # Data format format for ggplot (data frame for selectiviy)
  # creates a df with the first column named "Length" = plotlens
  # adds the selectivity matrix as additional columns (one column per mesh size)
    selectivity_df <- data.frame(Length = plotlens, rmatrix)
  
  # sets the column names for all columns except the first one which is Length (-1)
    colnames(selectivity_df)[-1] <- as.character(meshSizes)
    all_selectivity_curves[[model]] <- selectivity_df# stores the completed df in all_selectivity_curves
    # output ex: all_selectivity_curves$norm.loc
    
  # For ggplot format
    selectivity_long <- melt(selectivity_df, id.vars = "Length", variable.name = "MeshSize", value.name = "Selectivity")
    
  # gsub("Mesh_", "", selectivity_long$MeshSize) removes any "Mesh_" prefix from the mesh size values
    selectivity_long$MeshSize <- factor(gsub("Mesh_", "", selectivity_long$MeshSize))
    
  # setNames() assigns names to the elements of the rainbow color vector
    mesh_colors <- setNames(rainbow(nmeshes), levels(selectivity_long$MeshSize))
  
    
  # Plotting the estimated selectivity curve for each mesh size  
    p1 <- ggplot(selectivity_long, aes(x = Length, y = Selectivity, color = MeshSize, group = MeshSize)) +
      geom_line(linewidth = 1.5) +
      scale_color_manual(values = mesh_colors) +
      labs(title = paste("Selectivity Curve -", model), x = "Fish Length (cm)", y = "Relative Retention") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = "Mesh size [cm]"))
    
  # data frame for residuals
    dev_res_df <- data.frame(
      Length = rep(res$midLengths, times = nmeshes),
      MeshSize = factor(rep(meshSizes, each = length(res$midLengths))),
      Residuals = as.vector(res$Dev.resids)
    )
    
  # bubble plot - plotting the model residuals 
    p2 <- ggplot(dev_res_df, aes(x = Length, y = MeshSize, size = abs(Residuals), fill = Residuals)) +
      geom_point(shape = 21, color = "black", alpha = 0.8) + #creates circle points with
      scale_size(range = c(1, 8), guide = "none") + #controls the bubble sizes. Smallest =1, largest =8. guide = "none" hides the size legend
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + #scale_fill_gradient2(...) creates a diverging color scale: - blue, + red, near zero are white. midpoint = 0 centers the color scale at zero
      labs(title = paste("Deviance Residuals -", model), x = "Fish Length (cm)", y = "Mesh Size (cm)", fill = "Residual Value") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "right") +
      geom_hline(yintercept = meshSizes, linetype = "dotted", color = "gray60")#just for better visualization
 
  # combine the plots and save them   
    combined_plot <- grid.arrange(p1, p2, ncol = 1, heights = c(1, 1))

    jpeg_filename <- file.path(output_dir, paste0("Combined_", model, ".jpeg"))
    
    ggsave(filename = jpeg_filename,
               plot = combined_plot,
              width = 8, 
             height = 10,
              units = "in",
                dpi = 300)
    
    cat("Saved combined selectivity & residuals plot:", jpeg_filename, "\n")
  }
  
  return(list(results = results, summary_table = summary_table,selectivity_curves = all_selectivity_curves))
} 

#' Compare Selectivity Model Statistics
#'
#' Creates a formatted comparison table of fitted selectivity models with
#' key statistics and parameter estimates.
#'
#' @param result_object Output from fit_gillnet_dome function
#' @param include_bimodal Logical. Whether to include bimodal models in comparison.
#'   Default is FALSE
#' @param caption Character string for table caption.
#'   Default is "Comparison of Selectivity Models"
#' @param save_csv Logical. Whether to save results as CSV file. Default is FALSE
#' @param filename Character string for CSV filename. If NULL, generates automatic name
#'
#' @return A formatted kable table comparing model performance
#'
#' @details This function extracts key statistics from the fitted models and creates
#'   a comparison table ordered by model performance (highest log-likelihood first).
#'   For unimodal models, it shows Mode and Standard Deviation. For bimodal models,
#'   it additionally shows Mode2, StdDev2, and the proportion parameter P_Mode1.
#' @seealso \code{\link{fit_gillnet_dome}}, \code{\link{plot_mesh_curves}}
#' @examples
#' # First fit the models
#' data(raw_data_gillnet)
#' mesh_sizes <- c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19)
#' result <- fit_gillnet_dome(
#'   input_data = raw_data_gillnet,
#'   mesh_sizes = mesh_sizes,
#'   output_dir = tempdir(),
#'   length_seq = seq(40, 100, 1)
#' )
#' 
#' # Compare unimodal models only
#' summary_stats <- compare_stats(
#'   result_object = result,
#'   include_bimodal = FALSE,
#'   caption = "Comparison of Selectivity Models",
#'   save_csv = FALSE, 
#'   filename = NULL
#' )
#' print(summary_stats)
#' 
#' # Compare all models including bimodal
#' all_models <- compare_stats(
#'   result_object = result,
#'   include_bimodal = TRUE
#' )
#' print(all_models)
#' 
#' @export
compare_stats <- function(result_object, include_bimodal=FALSE, caption=
                              "Comparison of Selectivity Models", save_csv=FALSE, filename=NULL){
  
  
  # First get the summary table from the result object
  model_summary <- result_object$summary_table
  
  # Filter models
  if (!include_bimodal) {
    # Only include unimodal models
    models_to_include <- c("norm.loc", "norm.sca", "lognorm")
    model_comparison <- model_summary %>%
      filter(Model %in% models_to_include) %>%
      select(Model, LogLikelihood, Deviance, Mode1, StdDev1)
  } else {
    # Include all models
    model_comparison <- model_summary %>%
      select(Model, LogLikelihood, Deviance, Mode1, StdDev1, Mode2, StdDev2, P_Mode1)
  }
    
  # Order by best fit (highest LogLikelihood)
  model_comparison <- model_comparison %>% 
    arrange(desc(LogLikelihood))
  
  # Format for better readability
  formatted_comparison <- model_comparison %>%
    mutate(across(where(is.numeric), round, 2))
  
  # Create column names based on whether bimodal models are included
  if (!include_bimodal) {
    col_names <- c("Model", "Log-Likelihood", "Deviance", "Mode (cm)", "Std Dev (cm)")
  } else {
    col_names <- c("Model", "Log-Likelihood", "Deviance", "Mode1 (cm)", 
                   "StdDev1 (cm)", "Mode2 (cm)", "StdDev2 (cm)", "P_Mode1")
    colnames(formatted_comparison) <- col_names
  }
  
  # Save as CSV if requested
  if (save_csv) {
    # Generate default filename if none provided
    if (is.null(filename)) {
      if (!include_bimodal) {
        filename <- "unimodal_comparison.csv"
      } else {
        filename <- "all_models_comparison.csv"
      }
    }
    write.csv(formatted_comparison, filename, row.names = FALSE)
    cat("Table saved to:", filename, "\n")
  }
  
  # Create nice table using kableExtra
  formatted_comparison %>%
    kable(format = "markdown", 
          caption = caption,
          col.names = col_names) %>%
    kable_styling(full_width = FALSE)
  }
  
#' Calculate Composite Selectivity Curve
#' 
#' Calculates the aggregated selectivity curve across all mesh sizes...
#' 
#' @param model_results Model results object from fit_gillnet_dome results list
#' @param length_seq Numeric sequence of lengths for calculation
#' @return Data frame with Length and Selectivity columns containing the composite selectivity curve normalized to maximum value of 1
#' 
#' @details This function takes the selectivity parameters from a fitted model
#'   and calculates the combined selectivity across all mesh sizes. The resulting
#'   curve represents the overall selectivity pattern when using multiple mesh sizes
#'   simultaneously.
#' @seealso \code{\link{fit_gillnet_dome}}, \code{\link{plot_mesh_curves}}
#' @examples
#' # First fit the models
#' data(raw_data_gillnet)
#' mesh_sizes <- c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19)
#' result <- fit_gillnet_dome(
#'   input_data = raw_data_gillnet,
#'   mesh_sizes = mesh_sizes,
#'   output_dir = tempdir(),
#'   length_seq = seq(40, 100, 1)
#' )
#' 
#' # Get composite curve for the best fitting model
#' best_model <- result$results[[1]]  # First model in results (ordered by fit)
#' composite_curve <- get_composite_curve(
#'   model_results = best_model,
#'   length_seq = seq(40, 100, 0.5)
#' )
#' 
#' # View the composite selectivity curve
#' head(composite_curve)
#' 
#' # Plot the composite curve
#' library(ggplot2)
#' ggplot(composite_curve, aes(x = Length, y = Selectivity)) +
#'   geom_line(linewidth = 1.2, color = "blue") +
#'   labs(title = "Composite Selectivity Curve",
#'        x = "Fish Length (cm)", 
#'        y = "Relative Selectivity") +
#'   theme_minimal()
#' 
#' # Compare composite curves from different models
#' norm_loc_curve <- get_composite_curve(result$results$norm.loc, seq(40, 100, 0.5))
#' lognorm_curve <- get_composite_curve(result$results$lognorm, seq(40, 100, 0.5))
#'  
#' @export  
  get_composite_curve <- function(model_results, length_seq) {
  
  # Get the model type and parameters
  rtype <- model_results$rtype
  params <- model_results$par
  mesh_sizes <- model_results$meshSizes
  
  # Get the appropriate selectivity function
  sel_fun <- get("rtypes_Millar", envir = asNamespace("TropFishR"))(rtype)
  
  # Calculate selectivity matrix (rows = lengths, cols = mesh sizes)
  rmatrix <- outer(length_seq, mesh_sizes, sel_fun, params)
  
  # Sum across all mesh sizes
  comp <- rowSums(rmatrix)
  
  # Normalize to maximum of 1
  comp <- comp / max(comp, na.rm = TRUE)
  return(data.frame(Length = length_seq, Selectivity = comp))
  }

  
  
  #' Plot Individual and Composite Selectivity Curves
  #'
  #' Creates plots showing individual mesh selectivity curves and the composite
  #' curve for a specific model. Individual curves are shown in colors while the
  #' composite curve is shown as a black dashed line.
  #'
  #' @param result_object Output list from fit_gillnet_dome function containing
  #'   results and selectivity_curves elements
  #' @param model_name Character string specifying which model to plot 
  #'   (e.g., "norm.loc", "binorm.sca")
  #' @param length_seq Numeric sequence for plotting x-axis. Default is seq(40, 100, 1)
  #' @param save_plot Logical. Whether to save plot to file. Default is FALSE
  #' @param output_dir Character string for output directory when saving plot.
  #'   Required if save_plot = TRUE
  #'
  #' @return A ggplot object that can be displayed or further customized
  #' @seealso \code{\link{fit_gillnet_dome}}, \code{\link{get_composite_curve}}
  #' @examples
  #' 
  #' # First fit the models
  #' data(raw_data_gillnet)
  #' mesh_sizes <- c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19)
  #' result <- fit_gillnet_dome(
  #'   input_data = raw_data_gillnet,
  #'   mesh_sizes = mesh_sizes,
  #'   output_dir = tempdir(),
  #'   length_seq = seq(40, 100, 1)
  #' )
  #' 
  #' # Plot selectivity curves for the normal location model
  #' plot_mesh_curves(
  #'   result_object = result,
  #'   model_name = "norm.loc",
  #'   length_seq = seq(40, 100, 0.5)
  #' )
  #' 
  #' # Plot and save to file
  #' plot_mesh_curves(
  #'   result_object = result,
  #'   model_name = "lognorm",
  #'   length_seq = seq(40, 100, 1),
  #'   save_plot = TRUE,
  #'   output_dir = tempdir()
  #' )
  #' 
  #' # Compare different models by plotting multiple
  #' library(gridExtra)
  #' p1 <- plot_mesh_curves(result, "norm.loc")
  #' p2 <- plot_mesh_curves(result, "norm.sca") 
  #' grid.arrange(p1, p2, ncol = 1)
  #' 
  #' @export
    plot_mesh_curves <- function(result_object, model_name, length_seq = seq(40, 100, 1), 
                               save_plot = FALSE, output_dir = NULL) {
    # get the model results
    model_data <- result_object$selectivity_curves[[model_name]]
    
    # get mesh sizes (all column names except "Length")
    mesh_sizes <- setdiff(colnames(model_data), "Length") 
    
    # create empty data frame to store plotting data
    plot_data <- data.frame()
    
    # Add each individual mesh curve to the plot data frame
    for (mesh in mesh_sizes) {
      temp_data <- data.frame(
        Length = model_data$Length,
        Selectivity = model_data[[mesh]],
        Curve = paste0("Mesh ", mesh, " cm")
      )
      plot_data <- rbind(plot_data, temp_data)
    }
    
    # Calculate and add the composite curve using the same length sequence
    composite_curve <- get_composite_curve(result_object$results[[model_name]], length_seq)
    composite_data <- data.frame(
      Length = composite_curve$Length,
      Selectivity = composite_curve$Selectivity,
      Curve = "Composite"
    )
    plot_data <- rbind(plot_data, composite_data)
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = Length, y = Selectivity, color = Curve, 
                               linetype = (Curve == "Composite"), 
                               size = (Curve == "Composite"))) +
      geom_line() +
      scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"), guide = "none") +
      scale_size_manual(values = c("TRUE" = 1.5, "FALSE" = 0.8), guide = "none") +
      scale_color_manual(values = c(setNames(rainbow(length(mesh_sizes)), 
                                             paste0("Mesh ", mesh_sizes, " cm")),
                                    "Composite" = "black")) +
      labs(title = paste("Selectivity Curves for", model_name, "Model"),
           x = "Fish Length (cm)",
           y = "Relative Selectivity",
           color = "Curves") +
      theme_minimal() +
      theme(
        # Increase title font
        plot.title = element_text(size = 16, face = "bold"),
        
        # Increase axis title font
        axis.title = element_text(size = 14),
        
        # Increase axis text (numbers)
        axis.text = element_text(size = 12),
        
        # Increase legend title
        legend.title = element_text(size = 14),
        
        # Increase legend text
        legend.text = element_text(size = 12),
        
        # Keep legend position
        legend.position = "right"
      ) +
      xlim(min(length_seq), max(length_seq))
    
    # Save the plot if requested
    if (save_plot) {
      # Check if output_dir is provided
      if (is.null(output_dir)) {
        stop("When save_plot is TRUE, you must specify an output_dir")
      }
      
      # Make sure the directory exists
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("Created directory:", output_dir, "\n")
      }
      
      # Create filename
      filename <- file.path(output_dir, paste0("Curves_", model_name, ".jpeg"))
      
      # Save the plot
      ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 300, bg = "white")
      cat("Plot saved to:", filename, "\n")
    }
    
    # Return the plot
    return(p)
  }
  
  