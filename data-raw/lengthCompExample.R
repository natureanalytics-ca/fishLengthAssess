
devtools::load_all()
library(LBSPR)
library(fishSimGTG)

#Create a life history - same life history as in fishSimGTG:LifeHistoryExample
LifeHistoryExample<-new("LifeHistory")
LifeHistoryExample@Linf<-100
LifeHistoryExample@L50<-66
LifeHistoryExample@L95delta<-1
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
MyPars@L95 <- LifeHistoryExample@L50 + LifeHistoryExample@L95delta
MyPars@MK <- LifeHistoryExample@MK
MyPars@BinWidth <- 2
MyPars@Steepness <- LifeHistoryExample@Steep
MyPars@L_units <- LifeHistoryExample@L_units
MyPars@Walpha <- LifeHistoryExample@LW_A
MyPars@Walpha_units <- LifeHistoryExample@Walpha_units
MyPars@Wbeta <- LifeHistoryExample@LW_B
MyPars@FecB <- LifeHistoryExample@LW_B
MyPars@BinMin <- 10

#Year 1
MyPars@FM <- 1
MyPars@SL50 <- LifeHistoryExample@L50
MyPars@SL95 <- LifeHistoryExample@L50*1.15
sim<-LBSPRsim(MyPars)
dt<-data.frame(
  LMids = sim@LMids
)
dt$`2015` <- rmultinom(1, 158, sim@pLCatch[,1])
y1<-c(rep(sim@LMids, dt$`2015`), rep(NA, (483-158)))

#Year 2
MyPars@FM <- 1.4
MyPars@SL50 <- LifeHistoryExample@L50
MyPars@SL95 <- LifeHistoryExample@L50*1.15
sim<-LBSPRsim(MyPars)
dt$`2016` <- rmultinom(1, 200, sim@pLCatch[,1])
y2<-c(rep(sim@LMids, dt$`2016`), rep(NA, (483-200)))

#Year 3
MyPars@FM <- 1.6
MyPars@SL50 <- LifeHistoryExample@L50
MyPars@SL95 <- LifeHistoryExample@L50*1.15
sim<-LBSPRsim(MyPars)
dt$`2017` <- rmultinom(1, 235, sim@pLCatch[,1])
y3<-c(rep(sim@LMids, dt$`2017`), rep(NA, (483-235)))

#Year 4
MyPars@FM <- 1.8
MyPars@SL50 <- LifeHistoryExample@L50*0.8
MyPars@SL95 <- LifeHistoryExample@L50*0.8*1.15
sim<-LBSPRsim(MyPars)
dt$`2018` <- rmultinom(1, 378, sim@pLCatch[,1])
y4<-c(rep(sim@LMids, dt$`2018`), rep(NA, (483-378)))

#Year 5
MyPars@FM <- 2
MyPars@SL50 <- LifeHistoryExample@L50*0.6
MyPars@SL95 <- LifeHistoryExample@L50*0.6*1.15
sim<-LBSPRsim(MyPars)
dt$`2019` <- rmultinom(1, 483, sim@pLCatch[,1])
y5<-rep(sim@LMids, dt$`2019`)

#Frequency
LengthCompExampleFreq<-new("LengthComp")
LengthCompExampleFreq@title<-"Example length composition"
LengthCompExampleFreq@description<-"Simulated data using LBSPRsim"
LengthCompExampleFreq@L_units<-"cm"
LengthCompExampleFreq@L_type<-"TL"
LengthCompExampleFreq@L_source<-"FD"
LengthCompExampleFreq@dataType<-"Frequency"
LengthCompExampleFreq@observationGroup<-"Year"
LengthCompExampleFreq@header<-TRUE
LengthCompExampleFreq@dt<-dt
usethis::use_data(LengthCompExampleFreq, overwrite=TRUE)

#Length
LengthCompExampleLength<-new("LengthComp")
LengthCompExampleLength@title<-"Example length composition"
LengthCompExampleLength@description<-"Simulated data using LBSPRsim"
LengthCompExampleLength@L_units<-"cm"
LengthCompExampleLength@L_type<-"TL"
LengthCompExampleLength@L_source<-"FD"
LengthCompExampleLength@dataType<-"Length"
LengthCompExampleLength@observationGroup<-"Year"
LengthCompExampleLength@header<-TRUE

dt<-data.frame(cbind(y1, y2, y3, y4, y5))
colnames(dt)<-c("2015", "2016", "2017", "2018", "2019")
LengthCompExampleLength@dt<-dt
usethis::use_data(LengthCompExampleLength, overwrite=TRUE)
