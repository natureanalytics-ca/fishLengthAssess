

#---------------------------------
#Pool length data across groups
#----------------------------------

#Roxygen header
#'Utility that pools length data when multiple columns exist
#'
#' @param LengthCompObj  A life history object.
#' @param byGroup A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.
#' @export
#' @examples
#' poolLengthComp(fishLengthAssess::LengthCompExampleFreq)

poolLengthComp<-function(LengthCompObj, byGroup = FALSE) {

  #Check for errors in data entry
  if(length(LengthCompObj@dt) == 0 || length(LengthCompObj@dataType) != 1 ||  !(LengthCompObj@dataType %in%  c("Frequency", "Length"))) {
    return(NULL)
  } else {

    #Frequency
    if(LengthCompObj@dataType == "Frequency") {
      #Check there are at least two columns as there should be for Frequency data
      if(NCOL(LengthCompObj@dt) > 1) {
        #By group (or pooled?)
        if(byGroup){
          return(LengthCompObj@dt)
        } else {
          #Multiple groups to pool?
          if(NCOL(LengthCompObj@dt[,-1]) > 1) {
            return(data.frame(
              Lmids = LengthCompObj@dt[,1],
              Freq = rowSums(LengthCompObj@dt[,-1]))
            )
          } else {
            return(LengthCompObj@dt)
          }
        }
      } else {
        return(NULL)
      }
    }

    #Length
    if(LengthCompObj@dataType == "Length") {
      #Check there is at least 1 column as there should be for Frequency data
      if(NCOL(LengthCompObj@dt) > 0) {
        #By group (or pooled?)
        if(byGroup){
          return(LengthCompObj@dt)
        } else {
          return(data.frame(sort(as.vector(t(LengthCompObj@dt)))))
        }
      } else {
        return(NULL)
      }
    }
  }
}


#---------------------------------
#Optimum harvest length Beverton
#----------------------------------

#Roxygen header
#'Optimum harvest length Beverton (1992)
#'
#'Beverton, R.J.H., 1992. Patterns of reproductive strategy parameters in some marine teleost fishes. J. Fish. Biol. 41, 137-160.
#'
#' @param LifeHistoryObj  A life history object.
#' @import fishSimGTG
#' @importFrom methods is
#' @export
#' @examples
#' library(fishSimGTG)
#' LoptFunc(fishSimGTG::LifeHistoryExample)

LoptFunc<-function(LifeHistoryObj) {
  if(!is(LifeHistoryObj, "LifeHistory") || length(LifeHistoryObj@M) != 1 || length(LifeHistoryObj@Linf) != 1 || length(LifeHistoryObj@K) != 1) {
    return(NULL)
  } else {
    3*LifeHistoryObj@Linf/(3+LifeHistoryObj@M/LifeHistoryObj@K)
  }
}


#------------------------------------------
#Estimate length at full selectivity (mode)
#------------------------------------------

#Roxygen header
#'Estimate length at full selectivity using the mode of the length-frequency distribution
#'#'
#' @param LengthCompObj  A LengthComp object.
#' @param byGroup A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.
#' @importFrom stats loess
#' @importFrom stats predict
#' @export
#' @examples
#' library(fishSimGTG)
#' LcFunc(fishLengthAssess::LengthCompExampleFreq, byGroup = FALSE)

LcFunc<-function(LengthCompObj, byGroup = FALSE) {

  if(!is(LengthCompObj, "LengthComp") || length(LengthCompObj@dt) == 0 || length(LengthCompObj@dataType) != 1 ||  !(LengthCompObj@dataType %in%  c("Frequency", "Length"))) {
    return(NULL)
  } else {

    show_condition <- function(code) {
      tryCatch({
        x<-code
        c(x)
      },  error = function(c) NULL
      )
    }

    #Frequency data
    if(LengthCompObj@dataType == "Frequency") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      a<-as.numeric(dt[order(dt[,1]), 1])
      return(sapply(X=2:NCOL(dt), function(X){
        show_condition({
          tmp<-dt[order(dt[,1]), X]
          z1=cumsum(tmp)
          z1=z1/sum(tmp)
          d=loess(z1~a)
          a1=seq(min(a)+1, max(a),by=1)
          d1=predict(d,newdata=a1)
          d2=d1[2:length(a1)]-d1[1:(length(a1)-1)]
          d3=a1[1:length(d2)][d2==max(d2)]
          d3
        })
      }))
    }

    #Length data
    if(LengthCompObj@dataType == "Length") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      return(sapply(X=1:NCOL(dt), function(X){
        show_condition({
          z=table(dt[,X])
          z1=cumsum(z)
          z1=z1/sum(z)
          a=as.numeric(names(z))
          a1=seq(trunc(min(a))+1,max(a),by=1)
          d=loess(z1~a)
          d1=predict(d,newdata=a1)
          d2=d1[2:length(a1)]-d1[1:(length(a1)-1)]
          d3=a1[1:length(d2)][d2==max(d2)]
          d3
        })
      }))
    }
  }
}


#-----------------------------------------------------
#Estimate length at full selectivity (Kernel smoother)
#-----------------------------------------------------

#Roxygen header
#'Estimate length at full selectivity using a Kernel smoother
#'#'
#' @param LengthCompObj  A LengthComp object.
#' @param byGroup A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.
#' @importFrom stats density
#' @export
#' @examples
#' library(fishSimGTG)
#' modeKernelFunc(fishLengthAssess::LengthCompExampleFreq, byGroup = FALSE)

modeKernelFunc<-function(LengthCompObj, byGroup = FALSE) {

  if(!is(LengthCompObj, "LengthComp") || length(LengthCompObj@dt) == 0 || length(LengthCompObj@dataType) != 1 ||  !(LengthCompObj@dataType %in%  c("Frequency", "Length"))) {
    return(NULL)
  } else {

    show_condition <- function(code) {
      tryCatch({
        x<-code
        c(x)
      },  error = function(c) NULL
      )
    }

    #Frequency data
    if(LengthCompObj@dataType == "Frequency") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      LMids <- dt[,1]
      return(sapply(X=2:NCOL(dt), function(X){
        show_condition({
          frequencies <- dt[,X]
          freq_smooth <- density(x = rep(LMids,frequencies), bw="nrd0", na.rm=TRUE)
          round(freq_smooth$x[freq_smooth$y == max(freq_smooth$y)],0)
        })
      }))
    }

    #Length data
    if(LengthCompObj@dataType == "Length") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      return(sapply(X=1:NCOL(dt), function(X){
        show_condition({
          length_smooth <- density(dt[,X], bw="nrd0", na.rm = TRUE)
          round(length_smooth$x[length_smooth$y == max(length_smooth$y)],0)
        })
      }))
    }
  }
}

#-----------------------------------------------------
#Build SizeBins from LifeHistory
#-----------------------------------------------------

#Roxygen header
#'Build SizeBins list from a LifeHistory object
#'Creates the SizeBins list required by GTG dome-shaped LBSPR functions.
#'Uses CVLinf (0.1) and MaxSD (2) if not set as attributes on the LifeHistory object
#' CVLinf(0.1) and maxSD(2) ensure the length vector always extends beyond Linf to capture the full length distribution
#'
#' @param lifeHistoryObj  A LifeHistory S4 object with Linf defined
#' @param Linc Numeric. Length bin width in cm (default = 1).
#' @importFrom methods is
#' @export
#' @examples
#' library(fishSimGTG)
#' lh <- fishSimGTG::LifeHistoryExample
#' SizeBins <- make_size_bins(lh)

make_size_bins <- function(lifeHistoryObj, Linc = 1) {
  CVLinf <- if (!is.null(attr(lifeHistoryObj, "CVLinf"))) attr(lifeHistoryObj, "CVLinf") else 0.1
  MaxSD  <- if (!is.null(attr(lifeHistoryObj, "MaxSD")))  attr(lifeHistoryObj, "MaxSD")  else 2
  list(
    Linc   = Linc,
    ToSize = lifeHistoryObj@Linf + CVLinf * lifeHistoryObj@Linf * MaxSD
  )
}

#-----------------------------------------------------
#Get observed frequencies for plotting
#-----------------------------------------------------

#Roxygen header
#' Convert LengthComp data to frequency for plotting
#' Converts either Frequency or Length data from a LengthComp object into
#' a frequency distribution using the model length bins.
#' Allows observed vs predicted plots to be created regardless of the original(lengths or frecuency)
#' data type
#'
#' @param lengthCompObj A LengthComp S4 object with dataType "Frequency" or "Length".
#' @param SizeBins A list with Linc and ToSize elements (from make_size_bins()).
#' @return A data frame with columns Length and Observed.
#' @export
#' @examples
#' library(fishSimGTG)
#' SizeBins <- make_size_bins(fishSimGTG::LifeHistoryExample)
#' obs_df <- get_observed_freq(fishLengthAssess::LengthCompExampleFreq, SizeBins)

get_observed_freq <- function(lengthCompObj, SizeBins) {
  len_bins <- seq(0.5, SizeBins$ToSize - 0.5, by = SizeBins$Linc)

  if (lengthCompObj@dataType == "Frequency") {
    dt       <- lengthCompObj@dt
    freq_cols <- dt[, -1, drop = FALSE]
    observed  <- if (ncol(freq_cols) > 1) rowSums(freq_cols) else freq_cols[, 1]
    lengths  <- dt[, 1]

  } else if (lengthCompObj@dataType == "Length") {
    all_lengths <- unlist(lengthCompObj@dt)
    all_lengths <- all_lengths[!is.na(all_lengths)]
    observed    <- hist(all_lengths,
                        breaks = c(len_bins - 0.5, max(len_bins) + 0.5),
                        plot   = FALSE)$counts
    lengths     <- len_bins
  }

  data.frame(Length = lengths, Observed = observed)
}

#---------------------------------
# Plot pooled observed vs predicted
#---------------------------------

#' Plot pooled observed vs predicted catch length distribution
#'
#' Creates a bar and line plot comparing observed catch (grey bars) to
#' predicted catch (red line) from the pooled GTG dome LBSPR model fit.
#' Works for both Frequency and Length data types.
#'
#' @param results_all Output from \code{run_grouped_and_pooled()}.
#' @param lengthCompObj A LengthComp S4 object.
#' @param SizeBins A list with Linc and ToSize (from \code{make_size_bins()}).
#' @param output_dir Character. Folder to save the plot. Created if it does not exist.
#' @param filename Character. File name for the saved plot (default = "fit_pooled.jpeg").
#' @export
#' @examples
#' \dontrun{
#' plot_fit_pooled(results_all, lengthComp, SizeBins, output_dir = "plots")
#' }

plot_fit_pooled <- function(results_all, lengthCompObj, SizeBins,
                            output_dir = "plots", filename = "fit_pooled.jpeg") {
  if (!dir.exists(output_dir)) dir.create(output_dir)

  obs_df         <- get_observed_freq(lengthCompObj, SizeBins)
  length_vector  <- obs_df$Length
  observed_data  <- obs_df$Observed
  predicted_data <- results_all$pooled$PredLen

  compare_df <- data.frame(
    Length    = length_vector,
    Observed  = observed_data,
    Predicted = predicted_data
  )

  p <- ggplot(compare_df, aes(x = Length)) +
    geom_bar(aes(y = Observed), stat = "identity", fill = "grey70", alpha = 0.8) +
    geom_line(aes(y = Predicted), color = "red", linewidth = 1.2) +
    labs(
      title    = "Observed vs predicted catch length distribution",
      subtitle = paste0("F/M = ",  round(results_all$pooled$lbPars["F/M"],  3),
                        "   SPR = ", round(results_all$pooled$lbPars["SPR"], 3)),
      x = "Fish length (cm)",
      y = "Number of fish"
    ) +
    theme_bw()

  ggsave(file.path(output_dir, filename), plot = p, width = 7, height = 5, dpi = 300)
  cat("Plot saved to:", file.path(output_dir, filename), "\n")
  invisible(p)
}

#---------------------------------
# Plot observed vs predicted by group
#---------------------------------

#' Plot observed vs predicted catch length distribution by group
#'
#' Creates a faceted bar and line plot comparing observed and predicted catch
#' for each group individually. Works for both Frequency and Length data types.
#' Only runs when more than one group is present.
#'
#' @param results_all Output from \code{run_grouped_and_pooled()}.
#' @param lengthCompObj A LengthComp S4 object.
#' @param SizeBins A list with Linc and ToSize (from \code{make_size_bins()}).
#' @param sampleCatch The raw data frame loaded from CSV.
#' @param output_dir Character. Folder to save the plot. Created if it does not exist.
#' @param filename Character. File name for the saved plot (default = "fit_by_group.jpeg").
#' @export
#' @examples
#' \dontrun{
#' plot_fit_by_group(results_all, lengthComp, SizeBins, sampleCatch, output_dir = "plots")
#' }

plot_fit_by_group <- function(results_all, lengthCompObj, SizeBins, sampleCatch,
                              output_dir = "plots", filename = "fit_by_group.jpeg") {
  if (!dir.exists(output_dir)) dir.create(output_dir)

  grp_names <- names(results_all$grouped$group_results)

  # Only plot if more than one group
  if (length(grp_names) <= 1) {
    cat("Only one group present — by-group plot skipped.\n")
    return(invisible(NULL))
  }

  all_groups_df <- data.frame()

  for (grp in grp_names) {
    lengthComp_grp <- lengthCompObj

    if (lengthCompObj@dataType == "Frequency") {
      lengthComp_grp@dt <- sampleCatch[, c("Length", grp)]
    } else {
      lengthComp_grp@dt <- data.frame(sampleCatch[, grp, drop = FALSE])
    }

    obs_grp <- get_observed_freq(lengthComp_grp, SizeBins)

    temp_df <- data.frame(
      Length    = obs_grp$Length,
      Observed  = obs_grp$Observed,
      Predicted = results_all$grouped$group_results[[grp]]$PredLen,
      FM  = round(as.numeric(results_all$grouped$group_results[[grp]]$lbPars["F/M"]), 3),
      SPR = round(as.numeric(results_all$grouped$group_results[[grp]]$lbPars["SPR"]), 3),
      Group     = grp
    )
    all_groups_df <- rbind(all_groups_df, temp_df)
  }

  all_groups_df$label <- paste0(all_groups_df$Group,
                                "\nF/M = ", all_groups_df$FM,
                                "  SPR = ", all_groups_df$SPR)

  p <- ggplot(all_groups_df, aes(x = Length)) +
    geom_bar(aes(y = Observed), stat = "identity", fill = "grey70", alpha = 0.8) +
    geom_line(aes(y = Predicted), color = "red", linewidth = 1.0) +
    facet_wrap(~ label, ncol = 3) +
    labs(
      title = "Observed vs predicted catch length distribution — by group",
      x     = "Fish length (cm)",
      y     = "Number of fish"
    ) +
    theme_bw() +
    theme(strip.text = element_text(size = 8))

  ggsave(file.path(output_dir, filename), plot = p, width = 10, height = 6, dpi = 300)
  cat("Plot saved to:", file.path(output_dir, filename), "\n")
  invisible(p)
}

