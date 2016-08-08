###  COMPUTATION OF SPI TREND IN THE MEDITERRANEAN BASIN  ########
#
#  Edmondo Di Giuseppe
#  5 August 2016
#
#  this file is self consistent
##################################################################
# Dataset: CRU TS 3.23;
# Compute SPI following library(SPEI):
#Santiago Beguería and Sergio M. Vicente-Serrano (2013). SPEI: Calculation of the Standardised
#Precipitation-Evapotranspiration Index. R package version 1.6.
#https://CRAN.R-project.org/package=SPEI
###############################################################
library(raster)
library(SPEI)
library(gdata)
source("TrendFunctions.R")

### Load data from local directory:
#---------------------------------
### Prec
precCRU<-brick("/Users/edmondo/DATA/HOMER/Fda_shp_Locale/dati/CRU_TS_3.23/cru_ts3.23.1901.2014.pre.dat.nc")
MedArea<-extent(c(-9.05, 37, 29.95, 48.5)) 

#Crop over Mediterranean Area:
precCRUmed<-crop(precCRU,MedArea)
plot(precCRUmed,c(1,dim(precCRUmed)[3]))
rm(precCRU)
#---------------------------------

#Set times::
Xtimes<-getZ(precCRUmed)
TS.start<-c(as.numeric(getYear(Xtimes[1])),as.numeric(getMonth(Xtimes[1])))
LongTermBegin<-which(Xtimes=="1971-01-16")
LongTermEnd<-which(Xtimes=="2000-12-16")
#---------------
####   Analysis of 1971-2000 monthly precipitation  ####
precCRUmed.7100<-subset(precCRUmed,LongTermBegin:LongTermEnd)
plot(mean(precCRUmed.7100))

AveragedMonthlyPrec.7100<-ClinoMensili(precCRUmed.7100,median)
#---------------

###  SPI COMPUTATION ####
#------------------------

###  Cycling over scale times 3, 6, 12:

for (i in c(3,6,12)) {
  tutti1<-SPI_Raster(x=precCRUmed,filename=paste("SPI_",i,"_CRU.nc",sep=""),
                     na.rm=TRUE,scale=i,ts.start = TS.start,
                     ref.start = c(1971,1),ref.end = c(2000,12))
  }

plot(min(tutti1,na.rm=T))

#------------------------

