


#---------------------------------
#Froese - proportion of catch mature
#----------------------------------

#Roxygen header
#'Proportion of the catch that is above length at maturity
#'
#'Froese, R., 2004. Keep it simple: three indicators to deal with overfishing. Fish. Fish 5, 86-91. Proportion of the length observations that are equal to or greater than the length at 50% maturity parameter from a fishSimGTG LifeHistory object.
#'
#' @param LifeHistoryObj  A LifeHistory object from fishSimGTG.
#' @param LengthCompObj A LengthComp object
#' @param byGroup A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.
#' @import fishSimGTG
#' @export
#' @examples
#' library(fishSimGTG)
#' PmatFunc(fishSimGTG::LifeHistoryExample, fishLengthAssess::LengthCompExampleFreq, byGroup = FALSE)

PmatFunc<-function(LifeHistoryObj, LengthCompObj, byGroup = FALSE) {

  if(class(LifeHistoryObj) != "LifeHistory" ||
     class(LengthCompObj) != "LengthComp" ||
     length(LifeHistoryObj@L50) != 1 ||
     length(LengthCompObj@dt) == 0 ||
     length(LengthCompObj@dataType) != 1 ||
     !(LengthCompObj@dataType %in%  c("Frequency", "Length")) ||
     LengthCompObj@L_source == "FI"
  ) {
    return(NULL)
  } else {

    #Frequency data
    #Sum of bin counts >= L50
    if(LengthCompObj@dataType == "Frequency") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      binNum<-which(dt[,1]>=LifeHistoryObj@L50)
      return(sapply(X=2:NCOL(dt), function(X){
        sum(dt[binNum,X], na.rm = TRUE) / sum(dt[,X], na.rm=TRUE)
      }))
    }

    #Length data
    #Count obs. that meet criteria / total obs
    if(LengthCompObj@dataType == "Length") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      return(sapply(X=1:NCOL(dt), function(X){
        len = dt[,X]
        len=len[!is.na(len)]
        length(len[len>=LifeHistoryObj@L50])/length(len)
      }))
    }
  }
}


#-------------------------------------------------------------------
#Froese - proportion of specimens within optimum length range in the catch
#-------------------------------------------------------------------

#Roxygen header
#'Proportion of specimens within optimum length range in the catch
#'
#'Froese, R., 2004. Keep it simple: three indicators to deal with overfishing. Fish. Fish 5, 86-91. Proportion of the length observations that are equal to or greater than 0.9*Lopt and equal to or less than 1.1*Lopt. Lopt obtained from LoptFunc
#'
#' @param LifeHistoryObj  A LifeHistory object from fishSimGTG.
#' @param LengthCompObj A LengthComp object
#' @param byGroup A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.
#' @import fishSimGTG
#' @export
#' @examples
#' library(fishSimGTG)
#' PoptFunc(fishSimGTG::LifeHistoryExample, fishLengthAssess::LengthCompExampleFreq, byGroup = FALSE)

PoptFunc<-function(LifeHistoryObj, LengthCompObj, byGroup = FALSE) {

  Lopt<-LoptFunc(LifeHistoryObj) # calculates optimum harvest length

  if(class(LifeHistoryObj) != "LifeHistory" ||
     class(LengthCompObj) != "LengthComp" ||
     is.null(Lopt) ||
     length(LengthCompObj@dt) == 0 ||
     length(LengthCompObj@dataType) != 1 ||
     !(LengthCompObj@dataType %in%  c("Frequency", "Length")) ||
     LengthCompObj@L_source == "FI"
  ) {
    return(NULL)
  } else {

    #Frequency data
    #Sum of bin counts >= Lopt*0.9 and <= Lopt*1.1
    if(LengthCompObj@dataType == "Frequency") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      binNum<-which(dt[,1]>=Lopt*0.9 & dt[,1]<=Lopt*1.1)
      return(sapply(X=2:NCOL(dt), function(X){
        sum(dt[binNum,X], na.rm = TRUE) / sum(dt[,X], na.rm=TRUE)
      }))
    }

    #Length data
    #Count obs. that meet criteria / total obs
    if(LengthCompObj@dataType == "Length") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      return(sapply(X=1:NCOL(dt), function(X){
        len = dt[,X]
        len=len[!is.na(len)]
        length(len[len>=Lopt*0.9 & len<=Lopt*1.1])/length(len)
      }))
    }
  }
}


#--------------------------------------------------------------------------------------------------
#Froese - Proportion of mega-spawners, re-cast to 'avoid' mega spawners so that 100% is the target
#--------------------------------------------------------------------------------------------------

#Roxygen header
#'Proportion of mega-spawners, re-cast to 'avoid' mega spawners so that 100% is the target
#'
#'Froese, R., 2004. Keep it simple: three indicators to deal with overfishing. Fish. Fish 5, 86-91. 1 - proportion of the length observations that are equal to or greater than 1.1*Lopt. Lopt obtained from LoptFunc
#'
#' @param LifeHistoryObj  A LifeHistory object from fishSimGTG.
#' @param LengthCompObj A LengthComp object
#' @param byGroup A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.
#' @import fishSimGTG
#' @export
#' @examples
#' library(fishSimGTG)
#' PmegaFunc(fishSimGTG::LifeHistoryExample, fishLengthAssess::LengthCompExampleFreq, byGroup = FALSE)

PmegaFunc<-function(LifeHistoryObj, LengthCompObj, byGroup = FALSE) {

  Lopt<-LoptFunc(LifeHistoryObj)

  if(class(LifeHistoryObj) != "LifeHistory" ||
     class(LengthCompObj) != "LengthComp" ||
     is.null(Lopt) ||
     length(LengthCompObj@dt) == 0 ||
     length(LengthCompObj@dataType) != 1 ||
     !(LengthCompObj@dataType %in%  c("Frequency", "Length")) ||
     LengthCompObj@L_source == "FI"
  ) {
    return(NULL)
  } else {

    #Frequency data
    #Sum of bin counts >= L50
    if(LengthCompObj@dataType == "Frequency") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      binNum<-which(dt[,1]>=Lopt*1.1)
      return(sapply(X=2:NCOL(dt), function(X){
        1 - sum(dt[binNum,X], na.rm = TRUE) / sum(dt[,X], na.rm=TRUE)
      }))
    }

    #Length data
    #Count obs. that meet criteria / total obs
    if(LengthCompObj@dataType == "Length") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      return(sapply(X=1:NCOL(dt), function(X){
        len = dt[,X]
        len=len[!is.na(len)]
        1 - length(len[len>=Lopt*1.1])/length(len)
      }))
    }
  }
}

#--------------------------------------------------------------------------------------------------
#Proportion of catch above the minumum size limit (compliance)
#--------------------------------------------------------------------------------------------------

#Roxygen header
#'Proportion of catch above the minumum size limit
#'
#'Proportion of the length observations that are equal to or greater than a minimum size limit
#'
#' @param Lc  A minimum size limit in the same units of measure (L_units) and type of measurement (L_type) as the corresponding LengthComp object.
#' @param LengthCompObj A LengthComp object
#' @param byGroup A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.
#' @export
#' @examples
#' library(fishSimGTG)
#' PLcFunc(50, fishLengthAssess::LengthCompExampleFreq, byGroup = FALSE)


PLcFunc<-function(Lc, LengthCompObj, byGroup = FALSE) {

  if(class(LengthCompObj) != "LengthComp" ||
     !is.numeric(Lc) ||
     length(LengthCompObj@dt) == 0 ||
     length(LengthCompObj@dataType) != 1 ||
     !(LengthCompObj@dataType %in%  c("Frequency", "Length")) ||
     LengthCompObj@L_source == "FI"
  ) {
    return(NULL)
  } else {

    #Frequency data
    if(LengthCompObj@dataType == "Frequency") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      binNum<-which(dt[,1]>=Lc)
      return(sapply(X=2:NCOL(dt), function(X){
        sum(dt[binNum,X], na.rm = TRUE) / sum(dt[,X], na.rm=TRUE)
      }))
    }

    #Length data
    #Count obs. that meet criteria / total obs
    if(LengthCompObj@dataType == "Length") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      return(sapply(X=1:NCOL(dt), function(X){
        len = dt[,X]
        len=len[!is.na(len)]
        length(len[len>=Lc])/length(len)
      }))
    }
  }
}


#--------------------------------------------------------------------------------------------------
#LBSPR
#--------------------------------------------------------------------------------------------------

#Roxygen header
#'Estimation of SPR using LBSPR method (Hordyk et al. 2016)
#'
#'Hordyk, A.R., Ono, K., Prince, J.D., Walters, C.J., 2016. A simple length-structured model based on life history ratios and incorporating size-dependent selectivity: application to spawning potential ratios for data-poor stocks. Can. J. Fish. Aquat. Sci. 73, 1787–1799. https://doi.org/10.1139/cjfas-2015-0422
#'
#' @param LengthCompObj A LengthComp object
#' @param LifeHistoryObj  A LifeHistory object from fishSimGTG.
#' @param Lc  A minimum size limit in the same units of measure (L_units) and type of measurement (L_type) as the corresponding LengthComp object. Only used for analysis of fishery-independent length composition data.
#' @param binWidth  Used to create length-frequencies when LengthComp dataType is raw lengths.
#' @param cvLinf  Variability in Linf, used in LBSPR fitting. Default value same as that used in LBSPR library
#' @param byGroup A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.
#' @import fishSimGTG LBSPR
#' @importFrom graphics hist
#' @export
#' @examples
#' library(fishSimGTG)
#' library(LBSPR)
#' lbsprWrapper(fishSimGTG::LifeHistoryExample, fishLengthAssess::LengthCompExampleFreq)


lbsprWrapper<-function(LifeHistoryObj, LengthCompObj, Lc = 0, binWidth=1, cvLinf = 0.1, byGroup = FALSE, modtype = "GTG") {

  if(class(LifeHistoryObj) != "LifeHistory" ||
     class(LengthCompObj) != "LengthComp" ||
     !is.numeric(binWidth) ||
     !is.numeric(cvLinf) ||
     !is.logical(byGroup) ||
     !is.numeric(Lc) ||
     length(LifeHistoryObj@Linf) == 0 ||
     length(LifeHistoryObj@L50) == 0 ||
     length(LifeHistoryObj@MK) == 0 ||
     length(LengthCompObj@dt) == 0 ||
     length(LengthCompObj@dataType) != 1 ||
     !(LengthCompObj@dataType %in%  c("Frequency", "Length")) ||
     !(LengthCompObj@L_source %in%  c("FI", "FD")) ||
     LifeHistoryObj@Linf < 0 ||
     LifeHistoryObj@L50 < 0 ||
     LifeHistoryObj@L95delta <= 0 ||
     LifeHistoryObj@MK < 0 ||
     LifeHistoryObj@L50 >= LifeHistoryObj@Linf ||
     cvLinf < 0 ||
     LengthCompObj@L_source == "FI" & Lc < 0
  ) {
    return(NULL)
  } else {

    #-------------------------------
    #Create life history pars list
    #-------------------------------
    MyPars <- new("LB_pars")
    MyPars@CVLinf <- cvLinf
    MyPars@Linf <- LifeHistoryObj@Linf
    MyPars@L50 <- LifeHistoryObj@L50
    MyPars@L95 <- LifeHistoryObj@L50 + LifeHistoryObj@L95delta
    MyPars@MK <- LifeHistoryObj@MK
    if(length(LifeHistoryObj@LW_A) > 0) MyPars@Walpha <- LifeHistoryObj@LW_A
    if(length(LifeHistoryObj@LW_B) > 0) MyPars@Wbeta <- LifeHistoryObj@LW_B
    if(length(LifeHistoryObj@Steep) > 0) MyPars@Steepness<-LifeHistoryObj@Steep
    if(length(LifeHistoryObj@L_units) >0) MyPars@L_units <- LifeHistoryObj@L_units
    if(length(LifeHistoryObj@LW_A) >0)  MyPars@Walpha <- LifeHistoryObj@LW_A
    if(length(LifeHistoryObj@Walpha_units) >0)  MyPars@Walpha_units <-LifeHistoryObj@Walpha_units
    if(length(LifeHistoryObj@LW_B) >0)   MyPars@Wbeta <-  MyPars@FecB <- LifeHistoryObj@LW_B

    #-----------------
    #Data formatting
    #-----------------
    #Frequency data
    if(LengthCompObj@dataType == "Frequency") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      if(LengthCompObj@L_source == "FI") dt<-dt[which(dt[,1] >= Lc),] # only for fishery-independent
    }

    #If length data
    if(LengthCompObj@dataType == "Length") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      if(binWidth <= 0) binWidth <- 1
      Lbins_new<-seq(from=min(dt, na.rm=TRUE),
                     to=max(max(dt + binWidth, na.rm=TRUE), (MyPars@Linf + binWidth)),
                     by=binWidth
      )
      tmp<-sapply(X = 1:NCOL(dt), function(X){
        if(LengthCompObj@L_source == "FI") dtIn<-dt[which(dt[,1] >= Lc),X]
        if(LengthCompObj@L_source == "FD") dtIn<-dt[,X]
        dtIn=dtIn[!is.na(dtIn)]
        hist(dtIn,Lbins_new, plot=F)$count
      })
      Lmids_new<-Lbins_new[1:(NROW(Lbins_new)-1)]+binWidth/2
      dt<-cbind( Lmids_new,tmp)
    }

    #-----------------
    #Fitting
    #-----------------
    show_condition <- function(code) {
      tryCatch({
        x<-code
        x
      },  error = function(c) NULL
      )
    }

    Len<-new("LB_lengths", dataType="freq", verbose=FALSE)
    Len@LMids<-dt[,1]
    Len@LData<-as.matrix(dt[,-1])
    Len@NYears<-NCOL(dt[,-1])
    if(byGroup) {
      if(LengthCompObj@header) {
        if(LengthCompObj@dataType == "Frequency")  Len@Years<-colnames(LengthCompObj@dt[,-1])
        if(LengthCompObj@dataType == "Length") Len@Years<-colnames(LengthCompObj@dt)
      } else {
        if(LengthCompObj@dataType == "Frequency")  Len@Years<-seq(1, NCOL(LengthCompObj@dt[,-1]), 1)
        if(LengthCompObj@dataType == "Length") Len@Years<-seq(1, NCOL(LengthCompObj@dt), 1)
      }
    } else {
      Len@Years<- 1
    }
    return(show_condition(LBSPRfit(LB_pars = MyPars, LB_lengths = Len, Control=list(maxFM=10, modtype = modtype), verbose=FALSE)))
  }
}


#--------------------------------------------------------------------------------------------------
#Beverton-Holt mean length mortality estimator
#--------------------------------------------------------------------------------------------------

#Roxygen header
#'Calculate the equilibrium Beverton-Holt estimator of instantaneous total mortality from length data. Wrapper function from fishmethods R package.
#'
#'  Nelson GA (2023). fishmethods: Fishery Science Methods and Models. R package version 1.12-1, <https://CRAN.R-project.org/package=fishmethods>.
#'
#' @param LengthCompObj A LengthComp object
#' @param LifeHistoryObj  A LifeHistory object from fishSimGTG.
#' @param Lsel  Length at full selectivity. User can apply LcFunc or LcFuncKernel to the length data previously.
#' @param byGroup A logical indicating whether quantity is to be calculated separately for each of multiple length comp groups (TRUE) or to length comp is to be pooled across groups prior to calculating quantity (default = FALSE). When TRUE, pooling is ignored if only a single group exists.
#' @import fishSimGTG fishmethods
#' @export
#' @examples
#' library(fishSimGTG)
#' library(fishmethods)
#' bheqWrapper(fishSimGTG::LifeHistoryExample, fishLengthAssess::LengthCompExampleFreq, type=1,Lsel = 40)

bheqWrapper <- function (LifeHistoryObj, LengthCompObj, byGroup = FALSE, type, Lsel) {

  #Lc <- LcFunc(LengthCompObj) # not sure about this

  # add all checks for function arguments
  if(class(LifeHistoryObj) != "LifeHistory" ||
     class(LengthCompObj) != "LengthComp" ||
     !is.logical(byGroup) ||
     length(LifeHistoryObj@K) != 1 ||
     length(LifeHistoryObj@Linf) != 1 ||
     !is.numeric(Lsel) ||
     length(LengthCompObj@dt) == 0 ||
     length(LengthCompObj@dataType) != 1 ||
     !(LengthCompObj@dataType %in%  c("Frequency", "Length")) ||
     LengthCompObj@L_source == "FI" ||
     LifeHistoryObj@Linf < 0 ||
     LifeHistoryObj@K < 0
  ) {
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
          bh_freq <- bheq(rep(LMids,frequencies),type=type,K=LifeHistoryObj@K,Linf=LifeHistoryObj@Linf,Lc=Lsel)
        })
      }))
    }

    #Length data
    if(LengthCompObj@dataType == "Length") {
      dt<-poolLengthComp(LengthCompObj, byGroup)
      return(sapply(X=1:NCOL(dt), function(X){
        show_condition({
          bh_length <- bheq(dt[,X],type=type,K=LifeHistoryObj@K,Linf=LifeHistoryObj@Linf,Lc=Lsel)
          bh_length
        })
      }))

    }
  }
}


