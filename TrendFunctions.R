####          FUNCTIONS HUB FOR TREND ANALYSIS           #############
######################################################################
# 1) Put1 to cases with drought severity and 0 otherwise (single cell)
# 1bis) if selected,count of 10-years ahead moving window
# 2) verifying if it's a Poisson distribution 
# 3) computing frequency of cases during the entire period
# 4) identifying trend sign and its significance
# 5) extending trend calculation to raster
# 6) maps
#---------------------------------------------------------------------

### Function for single cell
# Put1 according to classes threshold-:: 
#Per le prove:
# cl=1
# ClassName<-substr(names(drought.classes[cl]),1,(nchar(names(drought.classes[cl]))-1))
# x<-getValues(datax1)[2,]  #solo NA
# x<-getValues(datax1)[318160,]  
# threshold=drought.classes[[cl]]
#------------

Put1_SingleCell<-function(x,threshold,xSum=FALSE){
  #   #threshold:levels for identifying drought class
  #   #xSum: if TRUE return vector sum (default =FALSE)
  x0<-numeric(length = length(x))
  x00<-na.exclude(x)
  if(any(x00 <= (-1))){
  x0[x<=threshold[2] & x>threshold[1]]<-1
  x<-x0
  }else{
  x<-rep(NA,length.out = length(x))
  }
  ifelse(xSum==FALSE,return(x),return(sum(x)))
  }
#(SingleCell<-Put1_SingleCell(x,threshold = threshold,xSum = FALSE))
#----------------------------------------------



##Test for Poisson distribution:
#(mu=sum(SingleCell))      #expected number of events in 35 years
#(lambda=mean(SingleCell))
##(expec=dpois(0:1,lambda=mean(SingleCell))*length(SingleCell))

# (Var=var(SingleCell))
# (R=Var/lambda)   #ratio
#(poisTest<-poisson.test(x=sum(SingleCell),T=length(SingleCell)))
#--------------------------------


### Power Law and Chisquare test
#plot.ecdf(SingleCell,do.points=FALSE,xval=1:length(SingleCell))
#plot(cumsum(SingleCell),type="l",xlab="Years",ylab=paste("cumulated number of ",ClassName,"events"))



### FUNCTION TO EVALUATE SPI TREND ###
##   --from volcanoes notation--
# x<-getValues(datax1)[318160,]  
# x<-getValues(datax1)[2,]  #solo NA
#   # x is raster
#   # alpha: test acceptance level-- 0.05, 010, etc.  -- (default=0.05)
#   # threshold: levels for identifying drought class
#   # xSum: if TRUE return vector sum otherwise perform time trend test (default =FALSE)
trendSPI<-function(x,alpha=0.05,threshold,xSum = FALSE){
  if(xSum==TRUE){
    x<-Put1_SingleCell(x,threshold = threshold,xSum = TRUE)
    return(x)
  }else{
    x<-Put1_SingleCell(x,threshold = threshold,xSum = FALSE)
    if(length(x[is.na(x)])!=0){#==length(x)){
      TestResult<-NA
    }else{
      #no drought events, i.e. sequence of zero's:
      if(length(x[x==0])==length(x)){
        TestResult<-0
      }else{
        n=sum(x)
        t=length(x)
        S=0
        t_i<-which(x>0)
        for (i in 1:n) {S=S+log(t/t_i[i])}
        #parameter estimates:
        Beta_hat=n/S
        lambda_hat=n*Beta_hat/t   #(no. of droughts/year)
        #---
        #Chisquare test:
        Statistic=2*S
        LowerTail=qchisq(alpha/2, 2*n, ncp = 0, lower.tail = TRUE, log.p = FALSE)
        UpperTail=qchisq(alpha/2, 2*n, ncp = 0, lower.tail = FALSE, log.p = FALSE)
        if(Statistic > LowerTail & Statistic <= UpperTail) Test<-"Accept H0: no trend";TestResult<-0
        if(Statistic <= LowerTail | Statistic > UpperTail) {
          Test<-"Reject H0: significative trend"
          if(Beta_hat>1) TestResult=1
          if(Beta_hat<1) TestResult=(-1)
        }  
      }
    }
    return(TestResult)
  }
}

#(TrendTest<-trendSPI(x,alpha=0.05,threshold=threshold,xSum = FALSE))

####   FUNCTION FOR RASTER   #####
# "x" is a RasterStack or RasterBrick  
# "filename" is output filename
# "threshold" vector of 2 elements to delimit a bin for the events
# "xSum" TRUE to count the events in each bin, FALSE to perform trend test (default=F)
# "monthSPI" character, month name (not abbreviation) of SPI calculation
#x<-crop(datax1,extent(c(11.3,12.6,41.5,43.8)))
Trend_Raster<-function(x,filename,alpha=0.05,threshold,xSum=FALSE,monthSPI){  
  #out<-brick(x,values=F) #contenitore rasterBrick vuoto con le stesse caratteristiche di x
  out<-raster(x) #contenitore raster vuoto con le stesse caratteristiche di x
  
  out <- writeStart(out, filename, format="CDF",overwrite=TRUE)#,zname="time",zunit="months")    
  for (i in 1:nrow(out)) {   
    v <- getValues(x, i)
    v <- t(apply(v,MARGIN=1,FUN=trendSPI,alpha=alpha,threshold=threshold,xSum=xSum))  
    out <- writeValues(out,v , i) 
  }  #chiudo ciclo "i" sulle rows del raster x
  
  out <- writeStop(out)
  names(out)<-monthSPI
  return(out)  
}

# tutti1<-Trend_Raster(x,filename="prova.nc",alpha=0.05,threshold=drought.classes[[cl]],xSum = FALSE,
#                      monthSPI = monthSPI[m])
#----------------------------------------------------------------------


#####################################################
######   FUNCTION HUB for CRU  ######################
#::::::::::::::::::::::::::::::::::::::::::::::::::::
#### Function for 30-years monthly average:
#------------------------------------------
ClinoMensili<-function(db.raster,fun){  #db.raster=db su cui calcolare le medie; #funzione=c(media,mediana)   
  FUN<-match.fun(fun)
  clino.12mesi<-raster(db.raster)
  for (m in 1:12){
    mese.tutti<-subset(db.raster,seq(m,dim(db.raster)[3],by=12))
    clino.mese<-round(calc(mese.tutti,FUN),1)
    clino.12mesi<-addLayer(clino.12mesi,clino.mese)
  }
  names(clino.12mesi)<-month.abb
  return(clino.12mesi)
}
#-----------------------------------------

##  SPI COMPUTATION  #####
#x<-crop(precCRUmed,extent(c(11.0,14.5,41.5,44.5)))
#x<-crop(precCRUmed,extent(c(11.5,12.0,41.5,42.0)))   ### solo NA
#x<-getValues(x)[1,]
## Internal function for extracting fitted values from spi() function in SPEI package::
#   -- set reference period for estimation phase to 1971-2000 ---
SPI.fitted<-function(x,ts.start,scale,na.rm,ref.start,ref.end){
  #transform into ts object:
  x<-ts(x,start = ts.start,frequency = 12)
  y<-spi(x,scale=scale,na.rm=na.rm,ref.start=ref.start,ref.end=ref.end)$fitted
  return(y)
}

#SPI.fitted(x,ts.start = ts.start,scale = 3,na.rm=T,ref.start = c(1971,1),ref.end = c(2000,12))

####  SPI APPLIED OVER RASTER  ######
#  #ts.start=c(year,month)   initialize time series
#  #scale=1,3,6,12,48 for spi() 
#  #ref.start=c(year,month)  subsetting 30-years period for estimating Gamma parameters:
SPI_Raster<-function(x,filename,ts.start,scale,na.rm,ref.start,ref.end){  
  out<-brick(x,values=F) #contenitore rasterBrick vuoto con le stesse caratteristiche di x
  #out<-raster(x) #contenitore raster vuoto con le stesse caratteristiche di x
  
  out <- writeStart(out, filename, format="CDF",overwrite=TRUE)#,zname="time",zunit="months")    
  for (i in 1:nrow(out)) {   
    v <- getValues(x, i)
    v <- t(apply(v,MARGIN=1,FUN=SPI.fitted,scale=scale,na.rm=na.rm,
                 ts.start=ts.start,ref.start=ref.start,ref.end=ref.end))  
    out <- writeValues(out,v , i) 
  }  #chiudo ciclo "i" sulle rows del raster x
  
  out <- writeStop(out)
  #names(out)<-monthSPI
  return(out)  
}


#::::::::::::::::::::::::::::::::::::::::::::::::::




######################################################################
####   spi() and spei()  functions from  SPEI package   ##############  
######################################################################
##spei() is needed to apply spi()::

#################################################################

### OTHER ATTEMPT FUNCTIONS  ########
# ## Same function with the option of setting a moving window::
# Put1_SingleCell_MovWindow<-function(x,threshold,MovWindow){     
#   #threshold:levels for identifying drought class
#   #MovWindow: define number of years of moving window
#   x0<-numeric(length = length(x))
#   x00<-numeric(length = length(x))
#   names(x00)<-names(x)
#   x000<-na.exclude(x)
#   if(any(x000 <= (-1))){
#     x0[x<=threshold[2] & x>threshold[1]]<-1
#     for (i in MovWindow:length(x0)) {x00[i]<-sum(x0[i:(i-MovWindow)])}
#     x<-x00[MovWindow:length(x00)]
#   }else{
#     x<-rep(NA,length.out = (length(x)-MovWindow-1))
#   }
#   return(x)
# } 
# prova<-Put1_SingleCell_MovWindow(x,threshold = moderateT,MovWindow=Twindow)
# prova
#----------------------------------------


