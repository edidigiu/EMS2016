###  ANALISYS OF SPI TREND IN THE MEDITERRANEAN BASIN  ########
#
#  Edmondo Di Giuseppe
#  20 Luglio 2016
###############################################################
# 1) read SPI3(6)(12).nc data in gdrive  (period gen1981-may2016: 425 months)
# 2) select seasons:
#     SPI3  ->  Feb, May, Aug, Nov
#     SPI6  ->  Feb, May
#     SPI12 ->  Aug
# and four negative SPI classes:
#   - moderately dry (-1>SPI>-1.5)
#   - severely dry (-1.5>SPI>-2) 
#   - extremely dry (SPI< -2)
#   - extremely&severely dry (-1000>SPI>-1.5) 
#
# 3) counting of 10-years ahead moving window
# 4) trend analysis for 7*3 time series listed in 2)
# 5) test of trend significance:::
#     *1  for positive trend
#     *-1 for negative trend
#     *0  for no trend OR no drought events  !!!!!!
#---------------------------------------------------------------

#Load library and Functions:
#--------------------------
library(raster)
library(gdata)    # for getYears()
source("TrendFunctions.R")
#---------------


#Set data and stuff:
#------------------------------------------------
Tbegin<-as.Date(strptime("1981-01-15",format="%Y-%m-%d"))
Tend<-as.Date(strptime("2016-05-15",format="%Y-%m-%d"))
Tseries<-seq.Date(Tbegin,Tend,by="month")

# set Twindow ONLY if you want subdivide events into "Twindow periods":
Twindow<-5 #numero di anni per la moving window

#SPI time scale:
scaleSPI<-c()

#Test acceptance I type error:
alpha=0.10

#drought class selection:
moderateT<-c(-1.5,-1)
severeT<-c(-2,-1.5)
extremeT<-c(-1000,-2)
ExtremeSevereT<-c(-1000,-1.5)
DroughtT<-c(-1000,-1)

drought.classes<-list(moderateT,severeT,extremeT,ExtremeSevereT,DroughtT)
names(drought.classes)<-c("ModerateT","SevereT","ExtremeT","Severe+ExtremeT","DroughtT")
#------------------------------------------------


#------------------------------------------------
#Read data  (period gen1981-may2016: 425 months):
#------------------------------------------------

DataList<-list.files(paste(getwd(),"/Dati",sep=""))
for(i in 1:length(DataList)){
  
  if(grepl("12",DataList[i])==TRUE){scaleSPI=12; N.ts=1; monthSPI<-"August"}
  if(grepl("3",DataList[i])==TRUE){scaleSPI=3; N.ts=4; monthSPI<-c("February", "May", "August", "November")}
  if(grepl("6",DataList[i])==TRUE){scaleSPI=6; N.ts=2; monthSPI<-c("February", "May")}
  
  datax<-brick(paste(getwd(),"/Dati/",DataList[i],sep=""))
  datax<-setZ(datax,Tseries)
  names(datax)<-paste("X",Tseries,sep="")

  for (m in 1:N.ts) {
    #months selection in time series:
    #if(scaleSPI==12 | scaleSPI==6){begin<-(match(monthSPI[m],month.name)+12)}
    #if(scaleSPI==3){
    #  if(monthSPI[m]=="February"){
        begin<-(match(monthSPI[m],month.name)+12)
    #  }else{
    #    begin<-match(monthSPI[m],month.name)
    #  }
    #}
    single_seq<-seq(begin,dim(datax)[3],by=12)
    #-----
    
    #season selection:
    datax1<-datax[[single_seq]]
    datax1<-setZ(datax1,as.Date(substr(names(datax1),2,10),format = "%Y.%m.%d"))
    # (Nyears=dim(datax1)[3])
    # plot(datax1,1:4)

    
    #Cycling on drought classes:
    for (cl in 1:length(drought.classes)) {
      ClassName<-substr(names(drought.classes[cl]),1,(nchar(names(drought.classes[cl]))-1))
      
      # Counting events by means of plugged in function:
      EventsNum<-Trend_Raster(x=datax1,filename=paste("Outputs/EventsNumber_",monthSPI[m],"_SPI-",scaleSPI,"_DroughtClass_",ClassName,".nc",sep = ""),
                              threshold=drought.classes[[cl]],xSum=TRUE, monthSPI=monthSPI[m])
      #Maps (number of events)::
      pdf(paste("Plots/EventsNumber/SPI-",scaleSPI,"/EventsNumber_",monthSPI[m],"_SPI-",scaleSPI,"_DroughtClass_",ClassName,".pdf",sep = ""))
      plot(EventsNum,interpolate=FALSE, col=terrain.colors(11),breaks=c(0,1,2,3,4,5,6,7,8,9,10),
           main=paste(monthSPI[m]," SPI-",scaleSPI," events in ",ClassName," class",sep="")
           )
      mtext(paste(getYear(getZ(datax1)[1]),"-",getYear(getZ(datax1)[length(getZ(datax1))]),sep=""),line = -1.5)
      dev.off()
      #------
      
      # Testing trend by means of plugged in function:
      TrendTest<-Trend_Raster(x=datax1,filename=paste("Outputs/TrendTest_alpha",alpha,"_",monthSPI[m],"_SPI-",scaleSPI,"_DroughtClass_",ClassName,".nc",sep = ""),
                              alpha=alpha,threshold=drought.classes[[cl]],xSum=FALSE, monthSPI=monthSPI[m])
      
      #Maps (time-trend)::
      pdf(paste("Plots/TrendTest/SPI-",scaleSPI,"/TrendTest_alpha",alpha,"_",monthSPI[m],"_SPI-",scaleSPI,"_DroughtClass_",ClassName,".pdf",sep = ""))
      plot(TrendTest,interpolate=FALSE, col=c("green4","lightsteelblue","orangered" ),breaks=c(-1.1,-0.4,0.4,1.1),
           lab.breaks=c("decr","notrend","notrend","incr"),
           main=paste(monthSPI[m]," SPI-",scaleSPI," events in ",ClassName," class",sep="")
      )
      mtext(paste("trend ",getYear(getZ(datax1)[1]),"-",getYear(getZ(datax1)[length(getZ(datax1))])," (alpha=",alpha,")",sep=""),line = -1.5)
      dev.off()
      #------
      
      
    } #closing "cl" in drought.classes
  } #closing "m" in monthsSPI
} #closing "i" in DataList







  