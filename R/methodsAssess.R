


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
#' PmatFunc(fishSimGTG::LifeHistoryExample, fishLengthAssess::LengthCompExample, byGroup = FALSE)

PmatFunc<-function(LifeHistoryObj, LengthCompObj, byGroup = FALSE) {

  if(length(LifeHistoryObj@L50) != 1 || length(LengthCompObj@dt) == 0 || length(LengthCompObj@dataType) != 1 ||  !(LengthCompObj@dataType %in%  c("Frequency", "Length"))) {
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
#' PoptFunc(fishSimGTG::LifeHistoryExample, fishLengthAssess::LengthCompExample, byGroup = FALSE)

PoptFunc<-function(LifeHistoryObj, LengthCompObj, byGroup = FALSE) {

  Lopt<-LoptFunc(LifeHistoryObj)

  if(is.null(Lopt) || length(LengthCompObj@dt) == 0 || length(LengthCompObj@dataType) != 1 ||  !(LengthCompObj@dataType %in%  c("Frequency", "Length"))) {
    return(NULL)
  } else {

    #Frequency data
    #Sum of bin counts >= L50
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
#' PmegaFunc(fishSimGTG::LifeHistoryExample, fishLengthAssess::LengthCompExample, byGroup = FALSE)

PmegaFunc<-function(LifeHistoryObj, LengthCompObj, byGroup = FALSE) {

  Lopt<-LoptFunc(LifeHistoryObj)

  if(is.null(Lopt) || length(LengthCompObj@dt) == 0 || length(LengthCompObj@dataType) != 1 ||  !(LengthCompObj@dataType %in%  c("Frequency", "Length"))) {
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
#' @import fishSimGTG
#' @export
#' @examples
#' library(fishSimGTG)
#' PLcFunc(50, fishLengthAssess::LengthCompExample, byGroup = FALSE)


PLcFunc<-function(Lc, LengthCompObj, byGroup = FALSE) {

  if(!is.numeric(Lc) || length(LengthCompObj@dt) == 0 || length(LengthCompObj@dataType) != 1 ||  !(LengthCompObj@dataType %in%  c("Frequency", "Length"))) {
    return(NULL)
  } else {

    #Frequency data
    #Sum of bin counts >= L50
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
