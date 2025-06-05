#' @importFrom stats dnorm median dbeta optim qlogis quantile sd setNames
#' @importFrom utils write.csv
#' @importFrom grDevices rainbow
NULL

#' GTG Dome-Shaped LBSPR Simulation
#'
#' Simulates length-based spawning potential ratio (SPR) using a genetically-based
#' growth type group (GTG) approach, with support for dome-shaped selectivity patterns.
#'
#' @param lifeHistoryObj S4 life history object containing stock parameters:
#' \itemize{
#'   \item Linf: Asymptotic length (accessible via lifeHistoryObj@Linf)
#'   \item MK: Natural mortality to growth coefficient ratio (lifeHistoryObj@MK)
#'   \item L50: Length at 50% maturity (lifeHistoryObj@L50)
#'   \item L95delta: Delta from L50 to L95 (lifeHistoryObj@L95delta)
#'   \item LW_A: Weight-length relationship intercept (lifeHistoryObj@LW_A)
#'   \item LW_B: Weight-length relationship exponent (lifeHistoryObj@LW_B)
#'   \item Steep: Stock-recruitment steepness (lifeHistoryObj@Steep)
#'   \item R0: Unfished recruitment (lifeHistoryObj@R0)
#' }
#' Additional parameters can be set as attributes:
#' \itemize{
#'   \item NGTG: Number of growth type groups (default 13)
#'   \item GTGLinfBy: Increment between Linf values (alternative to NGTG)
#'   \item CVLinf: Coefficient of variation in length-at-age (default 0.1)
#'   \item MaxSD: Maximum SD for Linf variation (default 2)
#'   \item FecB: Fecundity exponent (default 3)
#'   \item Mpow: Exponent for M/K ratio (default 0)
#' }
#' @param FleetPars A list containing fleet/selectivity parameters:
#' \itemize{
#'   \item selectivityCurve: Type of selectivity ("Logistic", "Normal.loc", "Normal.sca", 
#'         "logNorm", "binorm.sca", "bilognorm", "Knife")
#'   \item SL1-SL5: Selectivity parameters (depends on curve type)
#'   \item SLmesh: Mesh sizes (for gillnet models)
#'   \item SLMin: Minimum landing limit
#'   \item use_aggregated: Whether to use aggregated selectivity (default TRUE)
#'   \item fishery_mesh: Specific mesh size for fishery
#'   \item FM: Fishing mortality rate
#'   \item MLLKnife: Minimum length for knife-edge selectivity
#' }
#' @param SizeBins Optional list containing:
#' \itemize{
#'   \item Linc: Length bin increment
#'   \item ToSize: Maximum size bin
#' }
#'
#' @return A list containing simulation results:
#' \itemize{
#'   \item SPR: Spawning potential ratio
#'   \item Yield: Equilibrium yield
#'   \item YPR: Yield per recruit
#'   \item LCatchFished: Expected length composition of catch (fished)
#'   \item LPopFished: Expected length composition of population (fished)
#'   \item LCatchUnfished: Expected length composition of catch (unfished)
#'   \item LPopUnfished: Expected length composition of population (unfished)
#'   \item LenMids: Length bin midpoints
#'   \item VulLen2: Selectivity-at-length
#'   \item RelRec: Relative recruitment
#'   \item SPRatsize: Cumulative spawning-per-recruit by size
#' }
#'
#' @details This function simulates fish population dynamics using multiple growth 
#' type groups (GTGs) to represent variability in maximum length (Linf). It supports
#' various selectivity patterns including dome-shaped curves commonly observed in
#' gillnet fisheries.
#'
#' The simulation proceeds by:
#' \enumerate{
#'   \item Setting up GTGs with different Linf values
#'   \item Distributing recruits across GTGs
#'   \item Computing selectivity-at-length based on the specified curve
#'   \item Simulating population dynamics for fished and unfished scenarios
#'   \item Calculating spawning potential ratio and yield metrics
#' }
#'
#' @seealso \code{\link{DoOptDome}}, \code{\link{OptFunDome}}, 
#'   \code{\link{processLengthCompData}}
#' @examples
#' # example code
#' # GTG LBSPR simulation with Normal.sca selectivity
#' # Create LifeHistory object
#' LifeHistoryObj <- new("LifeHistory")
#' LifeHistoryObj@Linf <- 120
#' LifeHistoryObj@K <- 0.2
#' LifeHistoryObj@L50 <- 60
#' LifeHistoryObj@L95delta <- 2
#' LifeHistoryObj@MK <- 1.5
#' LifeHistoryObj@LW_A <- 0.01
#' LifeHistoryObj@LW_B <- 3
#' LifeHistoryObj@Steep <- 0.7
#' LifeHistoryObj@R0 <- 1E6
#' 
#' # Add GTG attributes
#' attr(LifeHistoryObj, "NGTG") <- 13
#' attr(LifeHistoryObj, "CVLinf") <- 0.1
#' attr(LifeHistoryObj, "MaxSD") <- 2
#' 
#' # Define fleet parameters
#' FleetPars <- list(
#'   FM = 1,
#'   selectivityCurve = "Normal.sca",
#'   SL1 = 55,
#'   SL2 = 18,
#'   SLmesh = c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19),
#'   use_aggregated = TRUE
#' )
#' 
#' # Define size bins
#' SizeBins <- list(Linc = 1, ToSize = 144)
#' 
#' # Run simulation
#' result <- GTGDomeLBSPRSim2(LifeHistoryObj, FleetPars, SizeBins)
#' @export
GTGDomeLBSPRSim2 <- function(lifeHistoryObj, FleetPars, SizeBins=NULL)  {

# Send all warning and status messages back to the main console screen where I can see them
sink(stdout(), type="message")  #diagnostic messages printed to the console for debugging
  
# Direct S4 access to lifeHistoryObj 
# Additional attribute-based settings via attr(): this allow to attach optional, flexible metadata to an S4 object
# without modifying the object class structure (basically to avoid interfere with the main structure)   
  
  NGTG <- attr(lifeHistoryObj, "NGTG")
  if (is.null(NGTG)) NGTG <- 13  # Default
  GTGLinfBy <- attr(lifeHistoryObj, "GTGLinfBy")
  if (is.null(GTGLinfBy)) GTGLinfBy <- NA
  Linf <- lifeHistoryObj@Linf
  CVLinf <- attr(lifeHistoryObj, "CVLinf")
  if (is.null(CVLinf)) CVLinf <- 0.1
  MaxSD <- attr(lifeHistoryObj, "MaxSD")
  if (is.null(MaxSD)) MaxSD <- 2
  MK <- lifeHistoryObj@MK
  L50 <- lifeHistoryObj@L50
  L95 <- lifeHistoryObj@L50 + lifeHistoryObj@L95delta  # Convert from delta
  Walpha <- lifeHistoryObj@LW_A
  Wbeta <- lifeHistoryObj@LW_B
  FecB <- attr(lifeHistoryObj, "FecB")
  if (is.null(FecB)) FecB <- 3
  Steepness <- lifeHistoryObj@Steep
  Mpow <- attr(lifeHistoryObj, "Mpow")
  if (is.null(Mpow)) Mpow <- 0
  R0 <- lifeHistoryObj@R0
  
  
  # Fleet parameters as a list
  # This section of the code include the options of dome - shape selectivity
  
  MLLKnife <- FleetPars$MLLKnife #minimum length for knife-edge sel
  selectivityCurve <- FleetPars$selectivityCurve #specifies selectivity function type ("Logistic", "Normal.loc", "binorm.sca", etc.)
  if (is.null(FleetPars$use_aggregated)) FleetPars$use_aggregated <- TRUE # whether to use all mesh sizes aggregated or a single one
  
  # extract selectivity-at-length parameters
  if(selectivityCurve == "Logistic"){
    SL50 <- FleetPars$SL1 # Length at 50% selectivity
    SL95 <- FleetPars$SL2 # Length at 95% selectivity
  
    } else if(selectivityCurve == "Normal.loc"){  # normal selectivity with fixed spread
    SLk <- FleetPars$SL1          # Mode of the normal curve
    SLsigma <- FleetPars$SL2      # Fixed spread (SD) of the normal curve
    SLmesh <- FleetPars$SLmesh    # Mesh sizes
    SLMin <- FleetPars$SLMin      # Minimum landing limit, if necessary
    if (is.null(SLMin)) SLMin <- NA
  
    } else if(selectivityCurve == "Normal.sca"){    # normal selectivity with proportional spread
    SLk1 <- FleetPars$SL1         # Modal
    SLk2 <- FleetPars$SL2         # Proportional (SD) spread (larger mesh = wider curve)
    SLmesh <- FleetPars$SLmesh    # Mesh sizes
    SLMin <- FleetPars$SLMin      # Minimum landing limit, if necessary
    if (is.null(SLMin)) SLMin <- NA
  
    } else if(selectivityCurve == "logNorm"){   # lognormal selectivity
    SLmu <- FleetPars$SL1         # Mean - Mu parameter for lognormal (in log space)
    SLsigma <- FleetPars$SL2      # SD - Sigma parameter for lognormal (in log space)
    SLmesh <- FleetPars$SLmesh    # Mesh sizes
    SLMin <- FleetPars$SLMin      # Minimum landing limit, if necessary
    if (is.null(SLMin)) SLMin <- NA
    
  # Adding new sel fucntions
  } else if(selectivityCurve == "binorm.sca"){      #two spreads, and a mixture weight (logit)       
    Mode1 <- FleetPars$SL1       # Mode of first normal component
    SD1 <- FleetPars$SL2         # Standard deviation of first normal component
    Mode2 <- FleetPars$SL3       # Mode of second normal component
    SD2 <- FleetPars$SL4         # Standard deviation of second normal component
    P1_logit <- FleetPars$SL5    # Logit of proportion for first component
    SLmesh <- FleetPars$SLmesh   # Mesh sizes
    SLMin <- FleetPars$SLMin;    # Minimum landing limit
    if (is.null(SLMin)) SLMin <- NA
    
  } else if(selectivityCurve == "bilognorm"){      #two spreads, and a mixture weight (logit)      
    SLmu1 <- FleetPars$SL1      # Mu parameter for first lognormal component (log space)
    SLsigma1 <- FleetPars$SL2   # Sigma parameter for first lognormal component (log space)
    SLmu2 <- FleetPars$SL3      # Mu parameter for second lognormal component (log space)
    SLsigma2 <- FleetPars$SL4   # Sigma parameter for second lognormal component (log space)
    P1_logit <- FleetPars$SL5   # Logit of proportion for first component
    SLmesh <- FleetPars$SLmesh  # Mesh sizes
    SLMin <- FleetPars$SLMin;   #Minimum landing limit
    if (is.null(SLMin)) SLMin <- NA
  } else if(selectivityCurve=="Knife"){                    # Knife-edge selectivity
  # Knife-edge selectivity: either 0 (below MLLKnife) or 1 (above MLLKnife)
    MLLKnife <- FleetPars$MLLKnife                         #All fish smaller than this value are not selected. Fish above this value are fully selected 
  }  
  
  # Fishing mortality rate
  FM <- FleetPars$FM            
  
  # Set up length bins and gtgs
  SDLinf <- CVLinf * Linf # Standard Deviation of Length-at-Age # Assumed constant CV here
  
  # This section define the length bins if they are not provided
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 1   # if SizeBins is not provided, length bin increment=1
    SizeBins$ToSize <- Linf + MaxSD * SDLinf # Maximum bin size based on largest Linf + SD
  }
  
  Linc <- SizeBins$Linc     # Length bin increment
  ToSize <- SizeBins$ToSize # Maximum size bin
  
  # Error Catches #
  # Fails if neither NGTG (number of GTGs) nor GTGLinfBy (Linf step) is defined
  # Makes sure GTG configuration exists. Sets default R0 if missing.
  if (!(exists("NGTG") | exists("GTGLinfBy"))) stop("NGTG or GTGLinfBy must be specified")
  if (!exists("R0")) R0 <- 1E6 # Default recruitment if not specified
  if (is.null(R0)) R0 <- 1E6
  
  # Set up Linfs for the different GTGs 
  # the sequence DiffLinfs  represents a range of Linf values centered on the mean Linf
  # and extending out by a multiple of standard deviations.
  #(if NGT exist do this, and GTGLinfBy doesnt exist do this)
  if (exists("NGTG") & !exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=NGTG) #a vector of different Linf values based on SD range and either step (GTGLinfBy) or count (NGTG).
    GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1] #calculates step size between Linf values
  #(if GTGLinfBy exist do this, and NGTG doesnt exist do this)
  } else  if (!exists("NGTG") & exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=GTGLinfBy)
    NGTG <- length(DiffLinfs) #determine Ngtg
  # both exist
  } else if (exists("NGTG") & exists("GTGLinfBy")) {
    #is not NA
    if (!is.na(GTGLinfBy)) {
      DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=GTGLinfBy)
      NGTG <- length(DiffLinfs)
    } 
    #is  NA
    if (is.na(GTGLinfBy)) {
      DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=NGTG)
      GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1]
    }  
  } 
  # Distribute Recruits across GTGS - divides recruits proportionally in GTG
  # A normal distribution centered at Linf, so most individuals fall near the mean, fewer in the extremes.
  # It is calculating recruitment probabilities (RecProbs) for each of the different maximum length (Linf) values in DiffLinfs vector
  RecProbs <- dnorm(DiffLinfs, Linf, sd=SDLinf) / sum(dnorm(DiffLinfs, Linf, sd=SDLinf)) 
  
  # Defines bin edges and midpoints: Length Bins (LenBins: Lower bounds of length classes; LenMids: Midpoints of bins) 
  if (is.null(ToSize)) ToSize <- max(DiffLinfs, Linf + MaxSD * SDLinf)
  LenBins <- seq(from=0, by=Linc, to=ToSize)
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=(length(LenBins)-1))#midpoints of length bins
  
  # Calculates fish WAL using an allometric relationship
  Weight <- Walpha * LenMids^Wbeta
  
  # Maturity and Fecundity for each GTG (each GTG matures at a different length range)
  # Each GTG has a different Linf, so maturity and fecundity vary by GTG
  # to get the absolute length at 50% maturity for each growth type group
  L50GTG <- L50/Linf * DiffLinfs #Each GTG has a different Linf, so L50 should be scaled accordingly
  L95GTG <- L95/Linf * DiffLinfs 
  DeltaGTG <- L95GTG - L50GTG
  MatLenGTG <- sapply(seq_along(DiffLinfs), function (X) 
    1.0/(1+exp(-log(19)*(LenMids-L50GTG[X])/DeltaGTG[X]))) #maturity at each length for each GTG.
  FecLenGTG <- MatLenGTG * LenMids^FecB # Fecundity across GTGs 
  
  #----------------------------------------------------#
  # Defining selectivity options in GTGDomeLBSPRSim2() #
  #----------------------------------------------------#
  
  if(selectivityCurve=="Logistic"){
    VulLen <- 1.0/(1+exp(-log(19)*((LenBins+0.5*Linc)-SL50)/((SL95)-(SL50)))) # Selectivity-at-Length
  }  
    else if(selectivityCurve=="Normal.sca"){
      # Normal Scale
      SLk1 <- FleetPars$SL1  # Modal length parameter
      SLk2 <- FleetPars$SL2  # Spread parameter
      
      if (!FleetPars$use_aggregated) {
        mesh_used <- FleetPars$fishery_mesh
        relsize <- mesh_used / min(SLmesh)
        VulLen <- exp(-((LenBins + 0.5 * Linc) - SLk1 * relsize)^2 / 
                        (2 * SLk2 * relsize^2))
      } else {
        VulLen <- 0
        for (j in seq_along(SLmesh)) {
          relsize <- SLmesh[j] / min(SLmesh)
          VulLen <- VulLen + exp(-((LenBins + 0.5 * Linc) - SLk1 * relsize)^2 / 
                                   (2 * SLk2 * relsize^2))
        }
      }
      
      if (!is.na(SLMin)) VulLen[LenBins < SLMin] <- 0
      VulLen <- VulLen / max(VulLen)
    }
  
  else if(selectivityCurve=="Normal.loc"){
    SLk <- FleetPars$SL1  # Modal fish capture length (scaled by mesh size)
    SLsigma <- FleetPars$SL2  # Fixed spread (SD) of the normal curve
    
    if (!FleetPars$use_aggregated) {
      mesh_used <- FleetPars$fishery_mesh
      relsize <- mesh_used / min(SLmesh)
      VulLen <- exp(-(((LenBins + 0.5 * Linc) - SLk * relsize)^2) / (2 * SLsigma^2))
    } else {
      VulLen <- 0
      for (j in seq_along(SLmesh)) {
        relsize <- SLmesh[j] / min(SLmesh)
        VulLen <- VulLen + exp(-(((LenBins + 0.5 * Linc) - SLk * relsize)^2) / (2 * SLsigma^2))
      }
    }
    
    if (!is.na(SLMin)) VulLen[LenBins < SLMin] <- 0
    VulLen <- VulLen / max(VulLen)
  }
    else if(selectivityCurve=="logNorm"){
      SLmu <- FleetPars$SL1    # Mu parameter (log space)
      SLsigma <- FleetPars$SL2 # Sigma parameter (log space)
      
      if (!FleetPars$use_aggregated) {
        mesh_used <- FleetPars$fishery_mesh
        relsize <- mesh_used / min(SLmesh)

        VulLen <- (relsize / (LenBins + 0.5 * Linc)) * exp(SLmu - SLsigma^2/2)
        VulLen <- VulLen * exp(-(log(LenBins + 0.5 * Linc) - SLmu - log(relsize))^2 / 
                                 (2 * SLsigma^2))
      } else {
        VulLen <- 0
        for (j in seq_along(SLmesh)) {
          relsize <- SLmesh[j] / min(SLmesh)
          
          seln <- (relsize / (LenBins + 0.5 * Linc)) * exp(SLmu - SLsigma^2/2)
          seln <- seln * exp(-(log(LenBins + 0.5 * Linc) - SLmu - log(relsize))^2 / 
                               (2 * SLsigma^2))
          
          VulLen <- VulLen + seln
        }
      }
      
      if(!is.na(SLMin)) VulLen[LenBins < SLMin] <- 0
      VulLen <- VulLen / max(VulLen)
    }
  
  else if(selectivityCurve == "binorm.sca"){
    Mode1 <- FleetPars$SL1
    SD1 <- FleetPars$SL2
    Mode2 <- FleetPars$SL3
    SD2 <- FleetPars$SL4
    P1_logit <- FleetPars$SL5
    P1 <- exp(P1_logit) / (1 + exp(P1_logit))
    
    if (!FleetPars$use_aggregated) {
      mesh_used <- FleetPars$fishery_mesh
      relsize <- mesh_used / min(SLmesh)
      
      seln1 <- exp(-0.5 * ((LenBins + 0.5 * Linc - Mode1 * relsize) / (SD1 * relsize))^2)
      seln2 <- exp(-0.5 * ((LenBins + 0.5 * Linc - Mode2 * relsize) / (SD2 * relsize))^2)
      
      VulLen <- P1 * seln1 + (1 - P1) * seln2
    } else {
      VulLen <- 0
      for (j in seq_along(SLmesh)) {
        relsize <- SLmesh[j] / min(SLmesh)
        
        seln1 <- exp(-0.5 * ((LenBins + 0.5 * Linc - Mode1 * relsize) / (SD1 * relsize))^2)
        seln2 <- exp(-0.5 * ((LenBins + 0.5 * Linc - Mode2 * relsize) / (SD2 * relsize))^2)
        
        VulLen <- VulLen + (P1 * seln1 + (1 - P1) * seln2)
      }
    }
    
    if (!is.null(SLMin) && !is.na(SLMin)) VulLen[LenBins < SLMin] <- 0
    if (FleetPars$use_aggregated || isTRUE(FleetPars$force_normalize)) {
      VulLen <- VulLen / max(VulLen, na.rm = TRUE)
    }
  }
  
  else if(selectivityCurve == "bilognorm"){
    SLmu1 <- FleetPars$SL1
    SLsigma1 <- FleetPars$SL2
    SLmu2 <- FleetPars$SL3
    SLsigma2 <- FleetPars$SL4
    P1_logit <- FleetPars$SL5
    P1 <- exp(P1_logit) / (1 + exp(P1_logit))
    
    if (!FleetPars$use_aggregated) {
      mesh_used <- FleetPars$fishery_mesh
      relsize <- mesh_used / min(SLmesh)
      
      seln1 <- (relsize / (LenBins + 0.5 * Linc)) * exp(SLmu1 - SLsigma1^2 / 2)
      seln1 <- seln1 * exp(-0.5 * ((log(LenBins + 0.5 * Linc) - SLmu1 - log(relsize)) / SLsigma1)^2)
      
      seln2 <- (relsize / (LenBins + 0.5 * Linc)) * exp(SLmu2 - SLsigma2^2 / 2)
      seln2 <- seln2 * exp(-0.5 * ((log(LenBins + 0.5 * Linc) - SLmu2 - log(relsize)) / SLsigma2)^2)
      
      VulLen <- P1 * seln1 + (1 - P1) * seln2
    } else {
      VulLen <- 0
      for (j in seq_along(SLmesh)) {
        relsize <- SLmesh[j] / min(SLmesh)
        
        seln1 <- (relsize / (LenBins + 0.5 * Linc)) * exp(SLmu1 - SLsigma1^2 / 2)
        seln1 <- seln1 * exp(-0.5 * ((log(LenBins + 0.5 * Linc) - SLmu1 - log(relsize)) / SLsigma1)^2)
        
        seln2 <- (relsize / (LenBins + 0.5 * Linc)) * exp(SLmu2 - SLsigma2^2 / 2)
        seln2 <- seln2 * exp(-0.5 * ((log(LenBins + 0.5 * Linc) - SLmu2 - log(relsize)) / SLsigma2)^2)
        
        VulLen <- VulLen + (P1 * seln1 + (1 - P1) * seln2)
      }
    }
    
    if (!is.null(SLMin) && !is.na(SLMin)) VulLen[LenBins < SLMin] <- 0
    if (FleetPars$use_aggregated || isTRUE(FleetPars$force_normalize)) {
      VulLen <- VulLen / max(VulLen, na.rm = TRUE)
    }
  }
  
  else if(selectivityCurve=="Knife"){    # knife-edge selectivity
    # Knife-edge selectivity (binary selection at MLLKnife)
    VulLen <- 0
    VulLen[(LenBins+0.5*Linc) < MLLKnife] <- 0 # Zero below minimum legal length
    VulLen[(LenBins+0.5*Linc) > MLLKnife] <- 1 # Fully selected above minimum legal length
    SL95 <- SL50 <- NA                         # Not used for knife-edge
  }

  # Add F-mortality below MLL
  SelLen <- VulLen # Selectivity is equal to vulnerability currently
  
  # Life-History Ratios 
  MKL <- MK * (Linf/(LenBins+0.5*Linc))^Mpow # M/K ratio for each length class
  MKMat <- matrix(rep(MKL, NGTG), nrow=length(MKL), byrow=FALSE)#Replicates MKL across each GTG

  FK <- FM * MK # F/K ratio
  FKL <- FK * SelLen # F/K ratio for each length class
  ZKLMat <- MKMat + FKL # Z/K ratio (total mortality) for each GTG

  # Initialize empty matrices 
  # number-per-recruit at length - simulates both unfished and fished populations using GTG
  NPRFished <- NPRUnfished <- matrix(0, nrow=length(LenBins), ncol=NGTG) 
  NatLUnFishedPop <- NatLFishedPop <- NatLUnFishedCatch <- 
  NatLFishedCatch <- FecGTGUnfished <- matrix(0, nrow=length(LenMids), 
                                                ncol=NGTG) # number per GTG in each length class 
  # Distribute Recruits into first length class
  NPRFished[1, ] <- NPRUnfished[1, ] <- RecProbs * R0 #initialize the first row of fished and unfished NPR
  for (L in 2:length(LenBins)) { # Calc number at each size class
    NPRUnfished[L, ] <- NPRUnfished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^MKMat[L-1, ]
    NPRFished[L, ] <- NPRFished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^ZKLMat[L-1, ]
    
  # finds the GTGs  cannot reach length bin L
    ind <- DiffLinfs  < LenBins[L]  #finds which GTGs can’t grow big enough to reach this length bin
  # for those GTGs, set the fished population to 0 in this bin.
    NPRFished[L, ind] <- 0
  # same
    NPRUnfished[L, ind] <- 0
  } 
  # clean up invalid or negative values
  NPRUnfished[is.nan(NPRUnfished)] <- 0
  NPRFished[is.nan(NPRFished)] <- 0
  NPRUnfished[NPRUnfished < 0] <- 0
  NPRFished[NPRFished < 0] <- 0
  
  # now, calc. the number of fish at each size
  for (L in 1:length(LenMids)) { 
    NatLUnFishedPop[L, ] <- (NPRUnfished[L,] - NPRUnfished[L+1,])/MKMat[L, ]
    NatLFishedPop[L, ] <- (NPRFished[L,] - NPRFished[L+1,])/ZKLMat[L, ]  
    FecGTGUnfished[L, ] <- NatLUnFishedPop[L, ] * FecLenGTG[L, ]
  }
  
  if(selectivityCurve=="Logistic"){
    VulLen2 <- 1.0/(1+exp(-log(19)*(LenMids-(SL50))/((SL95)-(SL50))))# Selectivity-at-Length
  }
  else if(selectivityCurve=="Normal.sca"){
    SLk1 <- FleetPars$SL1  # Modal length parameter
    SLk2 <- FleetPars$SL2  # Spread parameter
    
    if (!FleetPars$use_aggregated) {
      mesh_used <- FleetPars$fishery_mesh
      relsize <- mesh_used / min(SLmesh)
      VulLen2 <- exp(-(LenMids - SLk1 * relsize)^2 / 
                       (2 * SLk2 * relsize^2))
    } else {
      VulLen2 <- 0
      for (j in seq_along(SLmesh)) {
        relsize <- SLmesh[j] / min(SLmesh)
        VulLen2 <- VulLen2 + exp(-(LenMids - SLk1 * relsize)^2 / 
                                   (2 * SLk2 * relsize^2))
      }
    }
    
    if (!is.na(SLMin)) VulLen2[LenMids < SLMin] <- 0
    VulLen2 <- VulLen2 / max(VulLen2)
  }
  
  else if(selectivityCurve=="Normal.loc"){
    SLk <- FleetPars$SL1    # Modal fish capture length (scaled by mesh size)
    SLsigma <- FleetPars$SL2  # Fixed spread (SD) of the normal curve
    
    if (!FleetPars$use_aggregated) {
      mesh_used <- FleetPars$fishery_mesh
      relsize <- mesh_used / min(SLmesh)
      VulLen2 <- exp(-((LenMids - SLk * relsize)^2) / (2 * SLsigma^2))
    } else {
      VulLen2 <- 0
      for (j in seq_along(SLmesh)) {
        relsize <- SLmesh[j] / min(SLmesh)
        VulLen2 <- VulLen2 + exp(-((LenMids - SLk * relsize)^2) / (2 * SLsigma^2))
      }
    }
    
    if (!is.na(SLMin)) VulLen2[LenMids < SLMin] <- 0
    VulLen2 <- VulLen2 / max(VulLen2)
  }
  
  else if(selectivityCurve=="logNorm"){
    SLmu <- FleetPars$SL1    # Mu parameter (log space)
    SLsigma <- FleetPars$SL2 # Sigma parameter (log space)
    
    if (!FleetPars$use_aggregated) {
      mesh_used <- FleetPars$fishery_mesh
      relsize <- mesh_used / min(SLmesh)
      
      VulLen2 <- (relsize / LenMids) * exp(SLmu - SLsigma^2/2)
      VulLen2 <- VulLen2 * exp(-(log(LenMids) - SLmu - log(relsize))^2 / 
                                 (2 * SLsigma^2))
    } else {
      VulLen2 <- 0
      for (j in seq_along(SLmesh)) {
        relsize <- SLmesh[j] / min(SLmesh)
        
        seln <- (relsize / LenMids) * exp(SLmu - SLsigma^2/2)
        seln <- seln * exp(-(log(LenMids) - SLmu - log(relsize))^2 / 
                             (2 * SLsigma^2))
        
        VulLen2 <- VulLen2 + seln
      }
    }
    
    if(!is.na(SLMin)) VulLen2[LenMids < SLMin] <- 0
    VulLen2 <- VulLen2 / max(VulLen2)
  }
  
  else if(selectivityCurve == "binorm.sca"){
    Mode1 <- FleetPars$SL1
    SD1 <- FleetPars$SL2
    Mode2 <- FleetPars$SL3
    SD2 <- FleetPars$SL4
    P1_logit <- FleetPars$SL5
    P1 <- exp(P1_logit) / (1 + exp(P1_logit))
    
    if (!FleetPars$use_aggregated) {
      mesh_used <- FleetPars$fishery_mesh
      relsize <- mesh_used / min(SLmesh)
      
      seln1 <- exp(-0.5 * ((LenMids - Mode1 * relsize) / (SD1 * relsize))^2)
      seln2 <- exp(-0.5 * ((LenMids - Mode2 * relsize) / (SD2 * relsize))^2)
      
      VulLen2 <- P1 * seln1 + (1 - P1) * seln2
    } else {
      VulLen2 <- 0
      for (j in seq_along(SLmesh)) {
        relsize <- SLmesh[j] / min(SLmesh)
        
        seln1 <- exp(-0.5 * ((LenMids - Mode1 * relsize) / (SD1 * relsize))^2)
        seln2 <- exp(-0.5 * ((LenMids - Mode2 * relsize) / (SD2 * relsize))^2)
        
        VulLen2 <- VulLen2 + (P1 * seln1 + (1 - P1) * seln2)
      }
    }
    
    if (!is.null(SLMin) && !is.na(SLMin)) VulLen2[LenMids < SLMin] <- 0
    if (FleetPars$use_aggregated || isTRUE(FleetPars$force_normalize)) {
      VulLen2 <- VulLen2 / max(VulLen2, na.rm = TRUE)
    }
  }
  
  else if(selectivityCurve == "bilognorm"){
    SLmu1 <- FleetPars$SL1
    SLsigma1 <- FleetPars$SL2
    SLmu2 <- FleetPars$SL3
    SLsigma2 <- FleetPars$SL4
    P1_logit <- FleetPars$SL5
    P1 <- exp(P1_logit) / (1 + exp(P1_logit))
    
    if (!FleetPars$use_aggregated) {
      mesh_used <- FleetPars$fishery_mesh
      relsize <- mesh_used / min(SLmesh)
      
      seln1 <- (relsize / LenMids) * exp(SLmu1 - SLsigma1^2 / 2)
      seln1 <- seln1 * exp(-0.5 * ((log(LenMids) - SLmu1 - log(relsize)) / SLsigma1)^2)
      
      seln2 <- (relsize / LenMids) * exp(SLmu2 - SLsigma2^2 / 2)
      seln2 <- seln2 * exp(-0.5 * ((log(LenMids) - SLmu2 - log(relsize)) / SLsigma2)^2)
      
      VulLen2 <- P1 * seln1 + (1 - P1) * seln2
    } else {
      VulLen2 <- 0
      for (j in seq_along(SLmesh)) {
        relsize <- SLmesh[j] / min(SLmesh)
        
        seln1 <- (relsize / LenMids) * exp(SLmu1 - SLsigma1^2 / 2)
        seln1 <- seln1 * exp(-0.5 * ((log(LenMids) - SLmu1 - log(relsize)) / SLsigma1)^2)
        
        seln2 <- (relsize / LenMids) * exp(SLmu2 - SLsigma2^2 / 2)
        seln2 <- seln2 * exp(-0.5 * ((log(LenMids) - SLmu2 - log(relsize)) / SLsigma2)^2)
        
        VulLen2 <- VulLen2 + (P1 * seln1 + (1 - P1) * seln2)
      }
    }
    
    if (!is.null(SLMin) && !is.na(SLMin)) VulLen2[LenMids < SLMin] <- 0
    if (FleetPars$use_aggregated || isTRUE(FleetPars$force_normalize)) {
      VulLen2 <- VulLen2 / max(VulLen2, na.rm = TRUE)
    }
  }
  
  else if(selectivityCurve=="Knife"){   # knife-edge selectivity
    VulLen2 <- 0
    VulLen2[LenMids < MLLKnife] <- 0
    VulLen2[LenMids > MLLKnife] <- 1
    SL95 <- SL50 <- NA
  }
  
  NatLUnFishedCatch <- NatLUnFishedPop * VulLen2 # Unfished Vul Pop
  NatLFishedCatch <- NatLFishedPop * VulLen2 # Catch Vul Pop
  
  # Expected Length Structure - standardised 
  ExpectedLenCatchFished <- apply(NatLFishedCatch, 1, sum)/sum(apply(NatLFishedCatch, 1, sum)) # sums across GTGs for each length bin (row-wise).
  ExpectedLenPopFished <- apply(NatLFishedPop, 1, sum)/sum(apply(NatLFishedPop, 1, sum))       # same but using the pop not catch
  ExpectedLenCatchUnfished <- apply(NatLUnFishedCatch, 1, sum)/sum(apply(NatLUnFishedCatch, 1, sum)) #expected length composition of catch assuming no fishing mortality
  ExpectedLenPopUnfished <- apply(NatLUnFishedPop, 1, sum)/sum(apply(NatLUnFishedPop, 1, sum))      #length composition of the unfished population (natural mortality only)
  
  # Calc SPR
  EPR0 <- sum(NatLUnFishedPop * FecLenGTG) # Eggs-per-recruit Unfished
  EPRf <- sum(NatLFishedPop * FecLenGTG) # Eggs-per-recruit Fished
  SPR <- EPRf/EPR0 
  
  # Equilibrium Relative Recruitment
  recK <- (4*Steepness)/(1-Steepness) # Goodyear compensation ratio 
  reca <- recK/EPR0
  recb <- (reca * EPR0 - 1)/(R0*EPR0)
  RelRec <- max(0, (reca * EPRf-1)/(recb*EPRf))
  # RelRec/R0 - relative recruitment 
  YPR <- sum(NatLFishedPop  * Weight * VulLen2) * FM 
  Yield <- YPR * RelRec
  
  # Calc Unfished Fitness (total reproductive success (eggs) per recruit, for each GTG)
  Fit <- apply(FecGTGUnfished, 2, sum, na.rm=TRUE) # Total Fecundity per GTG Group
  FitPR <- Fit/RecProbs # Fitness per-recruit for each gtg
  FitPR <- FitPR/median(FitPR) #normalize for comparison
  ## Debugging
  # plot(FitPR, ylim=c(0,2)) # Should be relatively flat for equal fitness across GTG
  
  # how unequal is fitness?
  # sum of squared deviations from the median fitness value
  ObjFun <- sum((FitPR - median(FitPR, na.rm=TRUE))^2, na.rm=TRUE) # how much each GTG’s fitness deviates from the mediam
  #MKMat=the matrix of M/K ratios for each length bin and GTG
  #every value in this matrix should be positive 
  Pen <- 0; if (min(MKMat) <= 0 ) Pen <- (1/abs(min(MKMat)))^2 * 1E12
  
  ObjFun <- ObjFun + Pen
  
  # Calculate cumulative spawning-per-recruit (SPR) across size classes
  # N fish at each length across gtg x fec at each length for each gtgt
  SPRatsize <- cumsum(rowSums(NatLUnFishedPop * FecLenGTG)) #summ across gtg
  SPRatsize <- SPRatsize/max(SPRatsize) # normalize to max 1 [0 (no spawning yet) to 1 (100% of spawning done)]
  
  Output <- NULL 
  Output$SPR <- SPR
  Output$Yield <- Yield 
  Output$YPR <- YPR
  Output$LCatchFished <- ExpectedLenCatchFished
  Output$LPopFished <- ExpectedLenPopFished
  Output$LCatchUnfished <- ExpectedLenCatchUnfished
  Output$LPopUnfished <- ExpectedLenPopUnfished
  Output$NatLPopFished <- NatLFishedPop
  Output$NatLPopUnFish <- NatLUnFishedPop
  Output$NatLCatchUnFish <- NatLUnFishedCatch
  Output$NatLCatchFish <- NatLFishedCatch
  Output$LenBins <- LenBins
  Output$LenMids <- LenMids
  Output$NGTG <- NGTG
  Output$GTGdL <- DiffLinfs[2] - DiffLinfs[1]
  Output$DiffLinfs <- DiffLinfs
  Output$RecProbs <- RecProbs
  Output$Weight <- Weight
  Output$Winf <- Walpha * Linf^Wbeta
  Output$FecLen <- FecLenGTG 
  Output$MatLen <- MatLenGTG 
  Output$SelLen <- SelLen
  Output$MKL <- MKL
  Output$MKMat <- MKMat 
  Output$FKL <- FKL 
  Output$ZKLMat <- ZKLMat 
  Output$ObjFun <- ObjFun 
  Output$Pen <- Pen
  Output$FitPR <- FitPR
  Output$Diff <- range(FitPR)[2] - range(FitPR)[1]
  Output$L50GTG <- L50GTG 
  Output$L95GTG <- L95GTG
  Output$SPRatsize <- SPRatsize
  Output$RelRec <- RelRec
  Output$VulLen2 <- VulLen2 # added for comparison
  sink(type = "message")
  return(Output)
}

#' Optimization Objective Function for Dome-Shaped Selectivity
#'
#' Objective function for optimizing fishing mortality and optionally logistic selectivity parameters by comparing predicted vs observed length compositions.
#'
#' @param tryFleetPars Numeric vector of parameters to optimize (typically log(F/M)
#'   and optionally selectivity parameters)
#' @param fixedFleetPars List of fixed fleet parameters including selectivity type
#'   and other non-optimized selectivity parameters (i.e., dome-shaped selectivity)
#' @param LenDat Numeric vector of observed length composition data
#' @param lifeHistoryObj S4 life history object (same format as GTGDomeLBSPRSim2)
#' @param SizeBins Optional list of size bin specifications
#' @param mod Character string specifying model type ("GTG" or "LBSPR")
#'
#' @return Negative log-likelihood value (including penalties if applicable)
#'
#' @details This function serves as the objective function for optimization routines.
#' It calls the GTG simulation with trial parameters, compares predicted length
#' composition to observed data using multinomial likelihood, and returns the
#' negative log-likelihood.
#'
#' For logistic selectivity estimation, it includes penalty functions to prevent
#' biologically unreasonable parameter estimates (e.g., SL50 too close to Linf).
#' @seealso \code{\link{DoOptDome}}, \code{\link{GTGDomeLBSPRSim2}}
#' @keywords internal
#' @examples
#' # example code
#'  # Create LifeHistory object
#' LifeHistoryObj <- new("LifeHistory")
#' LifeHistoryObj@Linf <- 120
#' LifeHistoryObj@K <- 0.2
#' LifeHistoryObj@L50 <- 60
#' LifeHistoryObj@L95delta <- 2
#' LifeHistoryObj@MK <- 1.5
#' LifeHistoryObj@LW_A <- 0.01
#' LifeHistoryObj@LW_B <- 3
#' LifeHistoryObj@Steep <- 0.7
#' LifeHistoryObj@R0 <- 1E6
#' attr(LifeHistoryObj, "NGTG") <- 13
#' attr(LifeHistoryObj, "CVLinf") <- 0.1
#' attr(LifeHistoryObj, "MaxSD") <- 2
#'
#' # Fixed dome-shaped selectivity parameters
#' fixedFleetPars <- list(
#' selectivityCurve = "Normal.sca",
#' SL1 = 55,
#' SL2 = 18,
#' SLmesh = c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19),
#' use_aggregated = TRUE
#')
#'
#' # Trial parameter (log F/M)
#' tryFleetPars <- c(0.5)  # corresponds to F/M = exp(0.5) = 1.65
#' # Simulated observed data
#' data(gtg_catch_frequency)
#' LenDat <- gtg_catch_frequency$Catch_1

#' # Calculate negative log-likelihood
#' nll <- OptFunDome(tryFleetPars, fixedFleetPars, LenDat, LifeHistoryObj, 
#'                  SizeBins = list(Linc = 1, ToSize = 144), mod = "GTG")
#' @export
OptFunDome <- function(tryFleetPars, fixedFleetPars, LenDat, lifeHistoryObj, SizeBins=NULL, 
                       mod=c("GTG", "LBSPR")) {
  
  
  Fleet <- NULL
  Fleet$selectivityCurve <- fixedFleetPars$selectivityCurve
  if(Fleet$selectivityCurve=="Logistic"){
    if(length(tryFleetPars) == 3 & length(fixedFleetPars) == 1 & 
       c("selectivityCurve") %in% names(fixedFleetPars)){
      Fleet$SL1 <- exp(tryFleetPars[2]) * lifeHistoryObj@Linf 
      Fleet$SL2 <- Fleet$SL1  + (exp(tryFleetPars[3]) * lifeHistoryObj@Linf )
    } else {
      Fleet$SL1 <- fixedFleetPars$SL1
      Fleet$SL2 <- fixedFleetPars$SL2
    }
  }else if(Fleet$selectivityCurve=="Knife"){
    Fleet$MLLKnife <- fixedFleetPars$MLLKnife
    
  }else if(Fleet$selectivityCurve %in% c("Normal.sca", "Normal.loc", "logNorm",
                                         "binorm.sca","bilognorm")){
    Fleet$SL1 <- fixedFleetPars$SL1 
    Fleet$SL2 <- fixedFleetPars$SL2
    Fleet$SL3 <- fixedFleetPars$SL3
    Fleet$SL4 <- fixedFleetPars$SL4
    Fleet$SL5 <- fixedFleetPars$SL5
    Fleet$SLmesh <- fixedFleetPars$SLmesh
    
    #ensure that important optional settings from fixedFleetPars get correctly passed into the Fleet object used in the simulation
    if(!is.null(fixedFleetPars$SLMin)) Fleet$SLMin <- fixedFleetPars$SLMin
    if(!is.null(fixedFleetPars$use_aggregated)) Fleet$use_aggregated <- fixedFleetPars$use_aggregated
    if(!is.null(fixedFleetPars$fishery_mesh)) Fleet$fishery_mesh <- fixedFleetPars$fishery_mesh
  }
  
  #The only free parameter when dome selectivty
  Fleet$FM <- exp(tryFleetPars[1]) # changed to 1 from 3, as other parameters are fixed
  
  if (mod == "GTG") {
    runMod <- GTGDomeLBSPRSim2(lifeHistoryObj, Fleet, SizeBins)
  } else if (mod == "LBSPR") {
    stop("LBSPR mode is not supported with the current modifications")
  }
  
  
  LenDat <- LenDat + 1E-15 # is the  the observed length data - add tiny constant for zero catches
  LenProb <- LenDat/sum(LenDat) # proportions
  predProb <- runMod$LCatchFished # predicted proportionof catch at each length bin, not the number
  predProb <- predProb + 1E-15 # add tiny constant for zero catches
  NLL <- -sum(LenDat * log(predProb/LenProb)) ## Multinomial negative log-likelihood
  
  # an observaton: the real multinomial is NLL <- -sum(LenDat * log(predProb))
  # this form penalizes deviations between predicted and observed proportions
   
  #Penalty functions are used when model estimates the selectivity  
  # add penalty for SL50 when estimated and it becomes biologically unreasonable (too close to linf)
  if(length(tryFleetPars) == 3 & is.null(fixedFleetPars)){
    trySL50 <- exp(tryFleetPars[2])
    PenVal <- NLL
  
  # dbeta has a sharp peak near 1 
  # If trySL50 is close to 1 (i.e., near Linf : SL50 / Linf), the penalty is large
  # If trySL50 is smaller (more realistic selectivity), the penalty is small
  # The penalty ensures SL50 wont be unrealistically close to Linf
    Pen <- dbeta(trySL50, shape1=5, shape2=0.01) * PenVal  
    #if(!is.finite(NLL)) return(1E9 + runif(1, 1E4, 1E5))
    if (Pen == 0) {Pen <- PenVal * trySL50}
    # plot(xx, dbeta(xx, shape1=5, shape2=0.01) )
    NLL <- NLL+Pen
  }
  
  return(NLL)
}

#' Optimize Dome-Shaped Selectivity Parameters
#'
#' Main optimization function that estimates fishing mortality and optionally
#' selectivity parameters (only for logistic selectivity) by fitting to observed length composition data.
#'
#' @param lifeHistoryObj S4 life history object (same format as GTGDomeLBSPRSim2)
#' @param fixedFleetPars List of fleet parameters including selectivity type.
#'   If only selectivity type is provided for logistic selectivity, SL50 and SL95 will be estimated
#' @param LenDat Numeric vector of observed length composition (catch numbers by length)
#' @param SizeBins Optional list of size bin specifications  
#' @param mod Character string specifying model type ("GTG" or "LBSPR").
#'   Currently only "GTG" is supported
#'
#' @return A list containing:
#' \itemize{
#'   \item lbPars: Named vector of estimated parameters including F/M, SPR, 
#'         and selectivity parameters (if estimated)
#'   \item lbStdErrs: Standard errors for estimated parameters
#'   \item fixedFleetPars: Original fixed fleet parameters
#'   \item PredLen: Predicted length composition from fitted model
#'   \item NLL: Negative log-likelihood at optimum
#'   \item optimOut: Full output from optim() function
#'   \item MLE: Data frame with parameter estimates, initial values, and standard errors
#' }
#'
#' @details This function performs maximum likelihood estimation by:
#' \enumerate{
#'   \item Setting up optimization based on selectivity type
#'   \item Using optim() with either BFGS (multiple parameters) or Brent (F/M only)
#'   \item Calculating standard errors using the inverse Hessian
#'   \item Re-running the simulation with optimal parameters
#'   \item Returning comprehensive results
#' }
#'
#' For logistic selectivity, if only the selectivity type is specified, both
#' F/M and selectivity parameters (SL50, SL95) are estimated. 
#' For dome-shaped selectivity types, only F/M is estimated with selectivity parameters fixed.
#' 
#' @seealso \code{\link{GTGDomeLBSPRSim2}}, \code{\link{OptFunDome}}
#'
#' #' @examples
#' \donttest{
#' 
#' 
#' # Optimize F/M 
#' # Create LifeHistory object
#' LifeHistoryObj <- new("LifeHistory")
#' LifeHistoryObj@Linf <- 120
#' LifeHistoryObj@K <- 0.2
#' LifeHistoryObj@L50 <- 60
#' LifeHistoryObj@L95delta <- 2
#' LifeHistoryObj@MK <- 1.5
#' LifeHistoryObj@LW_A <- 0.01
#' LifeHistoryObj@LW_B <- 3
#' LifeHistoryObj@Steep <- 0.7
#' LifeHistoryObj@R0 <- 1E6
#' attr(LifeHistoryObj, "NGTG") <- 13
#' attr(LifeHistoryObj, "CVLinf") <- 0.1
#' attr(LifeHistoryObj, "MaxSD") <- 2
#'  
#' # Fixed dome-shaped selectivity
#' fixedFleetPars <- list(
#' selectivityCurve = "Normal.sca",
#' SL1 = 55,
#' SL2 = 18,
#' SLmesh = c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19),
#'    use_aggregated = TRUE
#'  )
#'  
#'  # Observed length composition data
#'  data(gtg_catch_frequency)
#'  LenDat <- gtg_catch_frequency$Catch_1
#'  
#'  # Run optimization
#'  result <- DoOptDome(LifeHistoryObj, fixedFleetPars, LenDat, 
#'                      SizeBins = list(Linc = 1, ToSize = 144), mod = "GTG")
#'  
#'  # Check results
#'  print(result$lbPars)
#'  }
#' @export
DoOptDome <- function(lifeHistoryObj, fixedFleetPars, LenDat, SizeBins=NULL, mod=c("GTG", "LBSPR")) {
  
 # setting up bins  
  CVLinf <- attr(lifeHistoryObj, "CVLinf")
  if (is.null(CVLinf)) CVLinf <- 0.1
  
  MaxSD <- attr(lifeHistoryObj, "MaxSD")
  if (is.null(MaxSD)) MaxSD <- 2
  
  SDLinf <- CVLinf * lifeHistoryObj@Linf
  
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 1
    SizeBins$ToSize <- lifeHistoryObj@Linf + MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) 
    SizeBins$ToSize <- lifeHistoryObj@Linf + MaxSD * SDLinf
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
  
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  
  # control parameters (for optimization)
  control_opt <- list(maxit = 500, reltol = 1e-8, REPORT = 10, trace = 1)
  
  selectivityCurve <- fixedFleetPars$selectivityCurve
  sFM <- 0.5  #starting value for f/M
  
  # detect whether the user wants to estimate SL50 and SL95 
  # or whether they are already provided and should be held fixed
  # If fixedFleetPars contains only $selectivityCurve = "Logistic", so length(fixedFleetPars) == 1 (TRUE)
  # so SL50 and SL95 are not provided, and should be estimated togehter with F/M
  if(fixedFleetPars$selectivityCurve=="Logistic" && length(fixedFleetPars) == 1){ # general logistic
    
    # Starting guesses
    sSL50 <- LenMids[which.max(LenDat)]/lifeHistoryObj@Linf  # Direct S4 access
    sDel <- 0.2 * LenMids[which.max(LenDat)]/lifeHistoryObj@Linf  # Direct S4 access
    Start <- log(c(sFM, sSL50, sDel))  #tryFleetPars
    # lowerBound <- c(-Inf, log(0.01), 0.0 ) # not used in BFGS optime
    # upperBound <- c(Inf, log(1+StockPars$CVLinf*StockPars$MaxSD), 1.0) # not used in BFGS optim
    
    methodOpt <- "BFGS"
    opt <- optim(par = Start, fn = OptFunDome, gr = NULL, 
                 fixedFleetPars=fixedFleetPars, LenDat=LenDat, lifeHistoryObj=lifeHistoryObj, SizeBins=SizeBins, mod=mod, 
                 method = methodOpt, control= list(maxit=500, abstol=1E-20),
                 hessian = TRUE)
    mleNames <-  c("log(F/M)", "SL50/Linf", "Sdelta/Linf")
    
    #uses the Brent optimization method when estimating F/M only
    
  } else{ # dome-shaped or fixed selectivity logistic
    Start <- log(c(sFM))  #tryFleetPars - it only optimizes F/M
    lowerBound <- -20
    upperBound <- 20
    methodOpt <- "Brent"
    opt <- optim(par = Start, fn = OptFunDome, gr = NULL, 
                 fixedFleetPars=fixedFleetPars, LenDat=LenDat, lifeHistoryObj=lifeHistoryObj, SizeBins=SizeBins, mod=mod, 
                 method = methodOpt, lower = lowerBound, upper = upperBound, 
                 control= list(maxit=500, abstol=1E-20),
                 hessian = TRUE)
    mleNames <- c("log(F/M)")
  }
  
  # negative log-likelihood
  newNLL<-opt$value # replaces objective in nlminb
  
  # variance-covariance matrix for std error calculation
  varcov <- solve(opt$hessian) # inverse of hessian matrix
  
  # maximum likelihood estimators and fishing parameters
  MLE <- data.frame(Parameter = mleNames, "Initial" = Start, "Estimate" = opt$par, "Std. Error" = diag(varcov),
                    check.names = FALSE)
  
  # back-transform MLE to obtain fishing parameters
  newFleet <- NULL 
  newFleet$selectivityCurve <- selectivityCurve
  newFleet$FM <- exp(opt$par[1])
  lbPars <- c("F/M"  = exp(opt$par[1]))
  if(fixedFleetPars$selectivityCurve=="Logistic"){
    if(length(fixedFleetPars) == 1){
      newFleet$SL1 <- exp(opt$par[2]) * lifeHistoryObj@Linf  # Direct S4 access
      newFleet$SL2 <- newFleet$SL1 + exp(opt$par[3]) * lifeHistoryObj@Linf  # Direct S4 access
      lbPars <- c(lbPars, 
                  "SL50" = exp(opt$par[2])*lifeHistoryObj@Linf,  # Direct S4 access
                  "SL95" = (exp(opt$par[2]) + exp(opt$par[3]))*lifeHistoryObj@Linf)  # Direct S4 access
    } else{
      newFleet$SL1 <- fixedFleetPars$SL1
      newFleet$SL2 <- fixedFleetPars$SL2
    }
  }else if(selectivityCurve=="Knife"){
    newFleet$MLLKnife <- fixedFleetPars$MLLKnife
  } else if(selectivityCurve %in% c("Normal.sca", "Normal.loc", "logNorm","binorm.sca", "bilognorm")){ # prescribed values, not optimised
    newFleet$SL1 <- fixedFleetPars$SL1
    newFleet$SL2 <- fixedFleetPars$SL2
    newFleet$SL3 <- fixedFleetPars$SL3
    newFleet$SL4 <- fixedFleetPars$SL4
    newFleet$SL5 <- fixedFleetPars$SL5
    newFleet$SLmesh <- fixedFleetPars$SLmesh
    if(!is.null(fixedFleetPars$SLMin)) newFleet$SLMin <- fixedFleetPars$SLMin
    if(!is.null(fixedFleetPars$use_aggregated)) newFleet$use_aggregated <- fixedFleetPars$use_aggregated
    if(!is.null(fixedFleetPars$fishery_mesh)) newFleet$fishery_mesh <- fixedFleetPars$fishery_mesh
  }
  
  # delta method to approximate standard error, CIs of estimates
  
  # ML estimators are log(F/M)
  sderr <- c(sqrt(exp(opt$par[1])*varcov[1,1]))
  names(sderr) = "F/M"
  if(fixedFleetPars$selectivityCurve=="Logistic" && length(fixedFleetPars) == 1){
    # log(SL50/Linf), log((SL95-SL50)/Linf) in log-space
    sderrSL50 <- sqrt((lifeHistoryObj@Linf*exp(opt$par[2]))^2*varcov[2,2])  # Direct S4 access
    sderrSL95 <- sqrt((lifeHistoryObj@Linf^2)*exp(opt$par[2])^2*varcov[2,2] +   # Direct S4 access
                        exp(opt$par[3])^2*varcov[3,3] + 2*exp(opt$par[3])*exp(opt$par[2])*varcov[2,3])
    sderr <-  c(sderr, SL50 = sderrSL50, SL95 = sderrSL95)
  } 
  
  # re-run the model using the new estimated parameters from optimization
  
  if (mod == "GTG") {
    runMod <- GTGDomeLBSPRSim2(lifeHistoryObj, newFleet, SizeBins)
  } else if (mod == "LBSPR") {
    stop("LBSPR mode is not supported with the current modifications")
  }
  
  lbPars <- c(lbPars, "SPR" = runMod$SPR)
  
  Out <- NULL 
  Out$lbPars <- lbPars      # fishing mortality, selectivity (where applicable), SPR
  Out$lbStdErrs <- sderr    # standard error for fishing mortality,... 
  Out$fixedFleetPars <- fixedFleetPars
  Out$PredLen <- runMod$LCatchFished * sum(LenDat)
  Out$NLL <- newNLL
  Out$optimOut <- opt
  Out$MLE <- MLE
  return(Out)
}

#' Optimize Dome Selectivity with LengthComp Objects
#'
#' Wrapper function for DoOptDome that handles S4 LengthComp objects with support
#' for both grouped and pooled analysis approaches.
#'
#' @param lifeHistoryObj S4 life history object (same format as GTGDomeLBSPRSim2)
#' @param fixedFleetPars List of fixed fleet parameters including selectivity type
#' @param LengthCompObj S4 LengthComp object containing length composition data
#' @param SizeBins Optional list of size bin specifications with Linc and ToSize elements
#' @param byGroup Logical, whether to analyze each group separately (TRUE) or pool all data (FALSE).
#'   Default is FALSE
#' @param Lc Length at first capture for removing small fish from analysis (default 0).
#'   Only applies to fishery-independent data
#' @param mod Character string specifying model type ("GTG" or "LBSPR"). 
#'   Currently only "GTG" is supported
#'
#' @return When byGroup = FALSE (pooled): Same as DoOptDome() plus processed_data and original_n_groups.
#'   When byGroup = TRUE (grouped): List containing:
#' \itemize{
#'   \item group_results: List of DoOptDome() results for each group
#'   \item processed_data: Processed length composition data
#'   \item n_groups: Number of groups analyzed
#' }
#'
#' @details This function provides a bridge between S4 LengthComp objects and the 
#' optimization routines. It automatically processes the S4 object using 
#' processLengthCompData() and then applies the appropriate optimization approach:
#' 
#' - If byGroup = TRUE and multiple groups exist: Fits separate models to each group
#' - If byGroup = FALSE or single group: Pools data and fits one model
#'
#' @seealso \code{\link{DoOptDome}}, \code{\link{processLengthCompData}}, 
#'   \code{\link{DoOptDome.aggregated}}, \code{\link{run_grouped_and_pooled}}
#' @examples
#' \donttest{
#' # Optimize using LengthComp object (pooled analysis)
#' # Create LifeHistory object
#' LifeHistoryObj <- new("LifeHistory")
#' LifeHistoryObj@Linf <- 120
#' LifeHistoryObj@K <- 0.2
#' LifeHistoryObj@L50 <- 60
#' LifeHistoryObj@L95delta <- 2
#' LifeHistoryObj@MK <- 1.5
#' LifeHistoryObj@LW_A <- 0.01
#' LifeHistoryObj@LW_B <- 3
#' LifeHistoryObj@Steep <- 0.7
#' LifeHistoryObj@R0 <- 1E6
#' attr(LifeHistoryObj, "NGTG") <- 13
#' attr(LifeHistoryObj, "CVLinf") <- 0.1
#' attr(LifeHistoryObj, "MaxSD") <- 2
#' 
#' # Create LengthComp object
#' data(gtg_catch_frequency)
#' LengthCompObj <- new("LengthComp",
#'                      dt = gtg_catch_frequency,
#'                      dataType = "Frequency",
#'                      L_source = "FD",
#'                      header = TRUE)
#' 
#' # Fixed fleet parameters
#' fixedFleetPars <- list(
#'   selectivityCurve = "Normal.sca",
#'   SL1 = 55,
#'   SL2 = 18,
#'   SLmesh = c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19),
#'   use_aggregated = TRUE
#' )
#' 
#' # Pooled analysis
#' result_pooled <- DoOptDome.LengthComp(LifeHistoryObj, fixedFleetPars, 
#'                                       LengthCompObj, byGroup = FALSE)
#' print(result_pooled$lbPars)
#' 
#' # By-group analysis
#' result_groups <- DoOptDome.LengthComp(LifeHistoryObj, fixedFleetPars, 
#'                                       LengthCompObj, byGroup = TRUE)
#' print(names(result_groups$group_results))
#' }
#' @examples
#' \donttest{
#' # Optimize using LengthComp object (pooled analysis)
#' # Create LifeHistory object
#' LifeHistoryObj <- new("LifeHistory")
#' LifeHistoryObj@Linf <- 120
#' LifeHistoryObj@K <- 0.2
#' LifeHistoryObj@L50 <- 60
#' LifeHistoryObj@L95delta <- 2
#' LifeHistoryObj@MK <- 1.5
#' LifeHistoryObj@LW_A <- 0.01
#' LifeHistoryObj@LW_B <- 3
#' LifeHistoryObj@Steep <- 0.7
#' LifeHistoryObj@R0 <- 1E6
#' attr(LifeHistoryObj, "NGTG") <- 13
#' attr(LifeHistoryObj, "CVLinf") <- 0.1
#' attr(LifeHistoryObj, "MaxSD") <- 2
#' 
#' # Create LengthComp object
#' data(gtg_catch_frequency)
#' LengthCompObj <- new("LengthComp",
#'                      dt = gtg_catch_frequency,
#'                      dataType = "Frequency",
#'                      L_source = "FD",
#'                      header = TRUE)
#' 
#' # Fixed fleet parameters
#' fixedFleetPars <- list(
#'   selectivityCurve = "Normal.sca",
#'   SL1 = 55,
#'   SL2 = 18,
#'   SLmesh = c(13.5, 14.0, 14.8, 15.4, 15.9, 16.6, 17.8, 19),
#'   use_aggregated = TRUE
#' )
#' 
#' # Pooled analysis
#' result_pooled <- DoOptDome.LengthComp(LifeHistoryObj, fixedFleetPars, 
#'                                       LengthCompObj, byGroup = FALSE)
#' print(result_pooled$lbPars)
#' 
#' # By-group analysis
#' result_groups <- DoOptDome.LengthComp(LifeHistoryObj, fixedFleetPars, 
#'                                       LengthCompObj, byGroup = TRUE)
#' print(names(result_groups$group_results))
#' }
#' @export
DoOptDome.LengthComp <- function(lifeHistoryObj, fixedFleetPars, LengthCompObj, SizeBins=NULL, 
                                 byGroup = FALSE, Lc = 0, mod=c("GTG", "LBSPR")) {
  

  
  # Process the S4 LengthComp object
  processed_data <- processLengthCompData(LengthCompObj, byGroup, SizeBins, Lc)
  
  # If multiple groups and not pooled - analyze each group separately
  if(processed_data$n_groups > 1 && !processed_data$was_pooled) {
    # Multiple groups - optimize each separately
    results <- list()
    
    for(i in 1:processed_data$n_groups) {
      cat(paste("Processing group:", processed_data$group_names[i], "\n"))
      
      # Extract data for this group (by group data)
      if(is.matrix(processed_data$LenDat)) {
        LenDat_group <- processed_data$LenDat[,i]
      } else {
      # Extract data for this group (pooled data)
        LenDat_group <- processed_data$LenDat
      }
      
      # Call your original DoOptDome function
      group_result <- DoOptDome(lifeHistoryObj, fixedFleetPars, LenDat_group, SizeBins, mod)
      group_result$group_name <- processed_data$group_names[i]
      results[[i]] <- group_result
    }
    
    names(results) <- processed_data$group_names
    return(list(
      group_results = results,
      processed_data = processed_data,
      n_groups = processed_data$n_groups
    ))
    
  } else {
    # Single group or pooled data
    result <- DoOptDome(lifeHistoryObj, fixedFleetPars, processed_data$LenDat, SizeBins, mod)
    
    # Add metadata about the original data
    result$processed_data <- processed_data
    result$original_n_groups <- if(is.matrix(LengthCompObj@dt)) {
      if(LengthCompObj@dataType == "Frequency") NCOL(LengthCompObj@dt) - 1 else NCOL(LengthCompObj@dt)
    } else {
      1       
    }
    
    return(result)
  }
}

#' Optimize with Aggregated Length Data
#'
#' Convenience wrapper that always pools (aggregates) multiple groups before 
#' optimization, regardless of the original data structure.
#'
#' @param lifeHistoryObj S4 life history object (same format as GTGDomeLBSPRSim2)
#' @param fixedFleetPars List of fixed fleet parameters including selectivity type
#' @param LengthCompObj S4 LengthComp object containing length composition data
#' @param SizeBins Optional list of size bin specifications with Linc and ToSize elements
#' @param Lc Length at first capture for removing small fish from analysis (default 0).
#'   Only applies to fishery-independent data
#' @param mod Character string specifying model type ("GTG" or "LBSPR"). 
#'   Currently only "GTG" is supported
#'
#' @return Same as DoOptDome() plus:
#' \itemize{
#'   \item processed_data: Processed and pooled length composition data
#'   \item original_n_groups: Number of groups in the original data before pooling
#' }
#'
#' @details This function ensures data is always pooled before optimization by
#' calling processLengthCompData() with byGroup = FALSE. It's useful when the user
#' want to analyze combined data from multiple sources/years without having to
#' remember to set byGroup = FALSE.
#'
#' @seealso \code{\link{DoOptDome}}, \code{\link{DoOptDome.LengthComp}}, 
#'   \code{\link{processLengthCompData}}
#'
#' @export
DoOptDome.aggregated <- function(lifeHistoryObj, fixedFleetPars, LengthCompObj, SizeBins=NULL, 
                                 Lc = 0, mod=c("GTG", "LBSPR")) {
  
  # Always aggregate (pool) multiple groups
  processed_data <- processLengthCompData(LengthCompObj, byGroup = FALSE, SizeBins, Lc)
  
  # Call your original DoOptDome function
  result <- DoOptDome(lifeHistoryObj, fixedFleetPars, processed_data$LenDat, SizeBins, mod)
  
  # Add metadata about the original data
  result$processed_data <- processed_data
  result$original_n_groups <- if(is.matrix(LengthCompObj@dt)) {
    if(LengthCompObj@dataType == "Frequency") NCOL(LengthCompObj@dt) - 1 else NCOL(LengthCompObj@dt)
  } else {
    1
  }
  
  return(result)
}

#' Run Both Grouped and Pooled Optimization
#'
#' Convenience function that runs both grouped and pooled optimization analyses
#' on the same dataset, allowing comparison between approaches.
#'
#' @param lifeHistoryObj S4 life history object (same format as GTGDomeLBSPRSim2)
#' @param fixedFleetPars List of fixed fleet parameters including selectivity type
#' @param LengthCompObj S4 LengthComp object containing length composition data
#' @param SizeBins Optional list of size bin specifications with Linc and ToSize elements
#' @param Lc Length at first capture for removing small fish from analysis (default 0).
#'   Only applies to fishery-independent data
#' @param mod Character string specifying model type (default "GTG"). 
#'   Currently only "GTG" is supported
#'
#' @return List containing:
#' \itemize{
#'   \item grouped: Results from DoOptDome.LengthComp() with byGroup = TRUE
#'   \item pooled: Results from DoOptDome.aggregated() 
#' }
#'
#' @details This function simplifies comparative analysis by automatically running:
#' \enumerate{
#'   \item Grouped analysis (each group fitted separately)
#'   \item Pooled analysis (all groups combined)
#' }
#' 
#' This allows users to compare whether pooling data affects parameter estimates
#' and model fit.
#'
#' @seealso \code{\link{DoOptDome.LengthComp}}, \code{\link{DoOptDome.aggregated}}
#'
#' @export
run_grouped_and_pooled <- function(lifeHistoryObj, fixedFleetPars, LengthCompObj, SizeBins = NULL, Lc = 0, mod = "GTG") {
  
  message("Running optimization by group...")
  grouped_results <- DoOptDome.LengthComp(lifeHistoryObj, fixedFleetPars, LengthCompObj, SizeBins, byGroup = TRUE, Lc = Lc, mod = mod)
  
  message("Running pooled optimization...")
  pooled_result <- DoOptDome.aggregated(lifeHistoryObj, fixedFleetPars, LengthCompObj, SizeBins, Lc = Lc, mod = mod)
  
  return(list(
    grouped = grouped_results,
    pooled = pooled_result
  ))
}
