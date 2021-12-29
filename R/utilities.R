

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
#' @export
#' @examples
#' library(fishSimGTG)
#' LoptFunc(fishSimGTG::LifeHistoryExample)

LoptFunc<-function(LifeHistoryObj) {
  if(length(LifeHistoryObj@M) != 1 | length(LifeHistoryObj@Linf) != 1 | length(LifeHistoryObj@K) != 1) {
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

  if(length(LengthCompObj@dt) == 0 || length(LengthCompObj@dataType) != 1 ||  !(LengthCompObj@dataType %in%  c("Frequency", "Length"))) {
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

