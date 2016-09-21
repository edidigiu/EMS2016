# EMS2016
SPI calculation from CRU TS 3.23 gridded dataset over Mediterranean basin and trend analysis

Standardized Precipitation Index (SPI) is widely used for drought quantification and monitoring. 
It was introduced in literature by @McKee1993 and is based on **precipitation** variable, exclusively.
For trend analysis, we consider the use of a special case of a nonhomogeneous Poisson process (NHPP)’s: the power law process defined in
@Crow1974a. The power law approach suggests the use of the $\chi^2$ test for time-trend analysis.

## Methodological steps
The analysis is done according to the following steps: 

1) calculate SPI index from precipitation values and generate SPI(3)(6)(12)_CRU.nc data (period **gen1901-dec2014**: 1328 months)

2) select seasons:

    - SPI3  ->  Feb, May, Aug, Nov
    
    - SPI6  ->  Feb, May
    
    - SPI12 ->  Aug

3)  and five negative SPI classes:

    - moderately dry (-1>SPI>-1.5)
    
    - severely dry (-1.5>SPI>-2)
    
    - extremely dry (SPI< -2)
    
    - extremely&severely dry (-1.5 >SPI> -1000)
    
    - drought risk (-1 >SPI> -1000)
    
 4) trend analysis for 7*5 time series listed in 2-3)
 
 5) test of trend significance


## Few technical details
Some details are given in each script files.


###Required libraries
```{r Library, include=FALSE}
library(RColorBrewer)
library(raster)
library(rasterVis)
library(gdata)
library(SPEI)
```

###From single cell to raster
We compute SPI following `library(SPEI)` developed by Santiago Beguería and Sergio M. Vicente-Serrano (2013). SPEI: Calculation of the Standardised
Precipitation-Evapotranspiration Index. R package version 1.6. https://CRAN.R-project.org/package=SPEI.
We develop some functions for applying SPI() from 'library(SPEI)' to a raster which are included into

```{r Functions, include=FALSE}
source("TrendFunctions.R")
```
###Output
The output is composed of as many .nc files as drought classes per each time scale, and each file is described as follows:

     - code 1  for positive trend
     - code -1 for negative trend
     - code 0  for no trend OR no drought events  
---------------------------------------------------------------
