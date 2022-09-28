
library(fishSimGTG)
library(tidyverse)
LifeHistoryObj<-new("LifeHistory")
LifeHistoryObj@MK<-1.5
LifeHistoryObj@Linf<-87
LifeHistoryObj@L50<-30
LifeHistoryObj@L95delta<-1


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

LifeHistoryObj<-LifeHistoryExample

LengthCompObj<-fishLengthAssess::LengthCompExampleLength

dt<-read_csv("D:/Cloud/TNC/Web app/poseidonApp/www/LengthComp/LengthComp_Length_noHeader_Year.csv", col_names = FALSE)
LengthCompObj@L_type<-"TL"
LengthCompObj@L_source<-"FD"
LengthCompObj@observationGroup<-"Year"
LengthCompObj@dataType = "Length"
LengthCompObj@header = FALSE
LengthCompObj@dt<-dt

lbsprWrapper(LifeHistoryObj, LengthCompObj, Lc = 0, binWidth=1, cvLinf = 0.1, byGroup = FALSE)
