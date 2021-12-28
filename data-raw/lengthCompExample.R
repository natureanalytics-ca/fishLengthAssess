
devtools::load_all()
library(LBSPR)
library(fishSimGTG)

#Create a life history - same life history as in fishSimGTG:LifeHistoryExample
LifeHistoryExample<-new("LifeHistory")
LifeHistoryExample@Linf<-100
LifeHistoryExample@L50<-66
LifeHistoryExample@L95<-67
LifeHistoryExample@MK<-1.5
LifeHistoryExample@LW_A<-0.01
LifeHistoryExample@LW_B<-3
LifeHistoryExample@title<-"Example life history"
LifeHistoryExample@shortDescription<-"Simulated life history of a fish based on B-H invariants"
LifeHistoryExample@speciesName<-"fish"
LifeHistoryExample@L_units<-"cm"
LifeHistoryExample@Walpha_units<-"g"
LifeHistoryExample@K<-0.2
LifeHistoryExample@M<-0.3
LifeHistoryExample@Tmax<- floor(-log(0.01)/0.3)
LifeHistoryExample@Steep<-0.8

#Use LBSPR sim to create a length comp
MyPars <- new("LB_pars")
MyPars@Linf <- LifeHistoryExample@Linf
MyPars@L50 <- LifeHistoryExample@L50
MyPars@L95 <- LifeHistoryExample@L95
MyPars@MK <- LifeHistoryExample@MK
MyPars@BinWidth <- 2
MyPars@Steepness <- LifeHistoryExample@Steep
MyPars@L_units <- LifeHistoryExample@L_units
MyPars@Walpha <- LifeHistoryExample@LW_A
MyPars@Walpha_units <- LifeHistoryExample@Walpha_units
MyPars@Wbeta <- LifeHistoryExample@LW_B
MyPars@FecB <- LifeHistoryExample@LW_B
MyPars@BinMin <- 30

#Year 1
MyPars@FM <- 1
MyPars@SL50 <- LifeHistoryExample@L50
MyPars@SL95 <- LifeHistoryExample@L50*1.8
sim<-LBSPRsim(MyPars)
dt<-data.frame(
  LMids = sim@LMids
)
dt$`2015` = sim@pLCatch[,1]*1000


#Year 2
MyPars@FM <- 1.4
MyPars@SL50 <- LifeHistoryExample@L50
MyPars@SL95 <- LifeHistoryExample@L50*1.6
sim<-LBSPRsim(MyPars)
dt$`2016` = sim@pLCatch[,1]*1000

#Year 3
MyPars@FM <- 1.6
MyPars@SL50 <- LifeHistoryExample@L50
MyPars@SL95 <- LifeHistoryExample@L50*1.15
sim<-LBSPRsim(MyPars)
dt$`2017` = sim@pLCatch[,1]*1000

#Year 4
MyPars@FM <- 1.8
MyPars@SL50 <- LifeHistoryExample@L50*0.8
MyPars@SL95 <- LifeHistoryExample@L50*1.15
sim<-LBSPRsim(MyPars)
dt$`2018` = sim@pLCatch[,1]*1000

#Year 5
MyPars@FM <- 1.8
MyPars@SL50 <- LifeHistoryExample@L50*0.6
MyPars@SL95 <- LifeHistoryExample@L50*1.15
sim<-LBSPRsim(MyPars)
dt$`2019` = sim@pLCatch[,1]*1000


#Frequency
LengthCompExample<-new("LengthComp")
LengthCompExample@title<-"Example length composition"
LengthCompExample@description<-"Simulated data using LBSPRsim"
LengthCompExample@L_units<-"cm"
LengthCompExample@L_type<-"TL"
LengthCompExample@L_source<-"FD"
LengthCompExample@dataType<-"Frequency"
LengthCompExample@observationGroup<-"Year"
LengthCompExample@header<-TRUE
LengthCompExample@dt<-dt
usethis::use_data(LengthCompExample, overwrite=TRUE)
