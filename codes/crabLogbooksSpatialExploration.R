##  crabLogbooksSpatialExploration.R
# Primary Author: Isaac Kaplan (May 1 2018)
# isaac.kaplan@noaa.gov
# Secondary Author: Emily Norton
# emilyln@uw.edu

# This script helps explore spatial autocorrelation through 1) variograms
# (https://gis.stackexchange.com/questions/109824/variogram-model-fit-compatability-among-geor-gstat-and-nlme-packages-in-r),
# 2) Moran's I statistic (https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/),
# and 3) concurvity test in mgcv (https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/concurvity.html)


#-------------------------------------------------------------------------------

#This program requires the following packages 
library('ncdf4')
library('tidyr')
library('geosphere')
library('tidyverse')
library('rgdal')
library('maptools')
library(gstat) # variogram model fitting
library(sp) # classes and methods for spatial data
library(geoR) # variogram model fitting
library(nlme) # linear mixed models with autocorrelation
library(ggplot2) #graphics
library(ape)
library(mgcv)

setwd('~/Documents/DungenessCPUEModel/Rcode/')


#--------------------------------------
# Load Observed fishing behaviors + ocn conditions
load("../Output/WA_OR_LogbookDataWithOceanConditionsHistoricalFuture_reallyreadyforGAM.RData")   # load df with all of the observed fishing behaviors from 2007-2018 

# Load GAM
save(gamDiS_SDL_less, file="../Output/gamDiS_SDL_less.RData")   # load GAM

gam.model <- gamDiS_SDL_less 
#----------------------------------

# TEST WITH SMALL SUBSET OF DATA
#LogbookDataWithOceanConditions <-LogbookDataWithOceanConditionsHistorical[1:5000,]

# Alternatively, use the full historical data set
LogbookDataWithOceanConditions <- LogbookDataWithOceanConditionsHistorical

#----------
# XY conversion from Lat Lon to meters (to ease interpreting of variograms)
#Note this leaves columns named as "Lons" and "Lats", which will be remedied farther below
zone <-10  
LonsLatsResids <-as.data.frame(cbind(gam.model$model$Lons,gam.model$model$Lats,gam.model$residuals)) #,colnames=c("Lats","Lons","residuals"))

colnames(LonsLatsResids) = c("Lons","Lats","residuals")

 coordinates(LonsLatsResids) <- c("Lons", "Lats")
  proj4string(LonsLatsResids) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(LonsLatsResids, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  XYResids <- as.data.frame(res)
 
  # Give columns the proper names (X and Y in meters, since they are no longer actually lons and lats)
  XYResids$Ymeters<- XYResids$Lats
  XYResids$Xmeters <- XYResids$Lons
   
   # Add jitter,  a small amount of noise so that no two points are directly overlapping. important for some spatial correlation functions ! 
   set.seed(456)
   XYResids$YmetersJitter <- jitter(XYResids$Ymeters)
   XYResids$XmetersJitter <- jitter(XYResids$Xmeters)
   
   # end of XY conversion from Lat Long
   #---------------


   #----------------------------------
   # Try some empirical variograms: https://gis.stackexchange.com/questions/109824/variogram-model-fit-compatability-among-geor-gstat-and-nlme-packages-in-r
   # Gstat empirical variogram

   # ADD A DUMMY GROUPING so we can use LME below: 
   XYResids$dGroup <- 1
  
   # Fitting a linear model to residuals of the GAM, because ultimately we want to see if there is 
   #   spatial autocorrelation in the GAM's residuals
   m1 <- lme(residuals ~ 1 , random = ~ 1|dGroup, data=XYResids) # MINI)
   
 #-------------- 
   
   # Use the simple residuals from the GAM
  dataframeXY <- data.frame(x=XYResids$XmetersJitter, y=XYResids$YmetersJitter, resid = XYResids$residuals)
  
  
  #define coordinates
   coordinates(dataframeXY) = ~x+y
   
   #--------------------------
   # VARIOGRAMS  
   #-------------------
   #CALCULATE THE EMPIRICAL VARIOGRAM using Gstat package
   v <- variogram(resid~1, dataframeXY)    # this takes some time
   dev.new()
   plot(v)
   # FIT THE VARIOGRAM SHAPE TO THE EMPIRICAL VARIOGRAM
     # note below that the call to vgm can look like this: vgm(1, "Exp", 300, 1) which specifies sill, model, range, nugget
  gstExp <- fit.variogram(v, vgm(model="Exp"), fit.method=1)
   # INSPECT THE FITTED VARIOGRAM -- SILL (ASYMPTOTE) AND RANGE (DISTANCE IN METERS Where variogram hits 95 PERCENT OF SILL)
   gstExp   
  # tiff("variogram.tiff")
   plot(v, model = gstExp, main = 'gstat')
#   dev.off()
   
   #--------------------------
   # VARIOGRAM, limited to 2000m or 2km only
   #-------------------
   # at scales of 2000 meters only, CALCULATE THE EMPIRICAL VARIOGRAM with one bin per 100m of distance
   v2 <- variogram(resid~1, dataframeXY, cutoff=2000,width=100)
   dev.new()
   plot(v2)
   # FIT THE VARIOGRAM SHAPE TO THE EMPIRICAL VARIOGRAM
   # note below that the call to vgm can look like this: vgm(1, "Exp", 300, 1) which specifies sill,model , range, nugget
   gstExp2 <- fit.variogram(v2, vgm(model="Exp"), fit.method=1)
   # INSPECT THE FITTED VARIOGRAM -- SILL (ASYMPTOTE) AND RANGE (DISTANCE IN METERS Where variogram hits 95 PERCENT OF SILL)
   gstExp2   
#   tiff("variogram_2000mCutoff_110521.tiff")
   plot(v2, model = gstExp2, main = 'gstat 2000m cutoff')
#dev.off()
  
   
 #----------------------------
   # Moran's I: 
#    https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/
#-----------------------------

   LonsLatsResids <-as.data.frame(cbind(gam.model$model$Lons,gam.model$model$Lats,gam.model$residuals)) 
   colnames(LonsLatsResids) = c("Lons","Lats","residuals")
   
   # Need to jitter so no 0 distances between points
   LonsLatsResids$Lons<-jitter(LonsLatsResids$Lons)
   LonsLatsResids$Lats<-jitter(LonsLatsResids$Lats)
   
   
   # SETTING UP SOME INFO BECAUSE WE will NEED TO RESAMPLE FROM THINNED DATA SET
   numPoints<-length(LonsLatsResids$Lats)  # counting rows of data
   p <-numeric(0)
   numReplicates <- 10   #  CHOOSE THE NUMBER OF REPLICATE TEST DATA SETS WE WANT
   fractionForTestSet <- 0.05  # CHOOSE THE FRACTION OF THE DATA TO RESERVE AS TEST DATA
   
   
   # NOW CALCULATE MORAN'S I , BUT FROM THINNED VERSIONS OF THE DATASET
   for(i in 1:numReplicates) {
   
   # Subset because can only handle a small thinned amount of data
   reshuffledPointNumbers <- sample(1:numPoints, numPoints, replace=F)  # reshuffle the points
   numPointsForTestSet  <- round(fractionForTestSet*numPoints) 
   LonsLatsResidsSubset <-LonsLatsResids[reshuffledPointNumbers[1:numPointsForTestSet],]
   
    # calculate distances between pairs of points
    resid.dists <- as.matrix(dist(cbind(LonsLatsResidsSubset$Lons, LonsLatsResidsSubset$Lats)))
   
    # invert the distance matrix. 
    resid.dists.inv <- 1/resid.dists
    diag(resid.dists.inv) <- 0   # set elements on diagonal to 0. 

    # calculate Moran's I statistic
    moranOut<-Moran.I(LonsLatsResidsSubset$residuals, resid.dists.inv,na.rm=TRUE)
    # pvalue less than 0.05 means  reject the null hypothesis that there is zero spatial autocorrelation present in the variable 
   
    print(moranOut$p.value)
    # save p.value into p
    p<-rbind(p,moranOut$p.value)

   }
   print(p)
   
   # end of Moran's I calculation. 
   #---------------------
   
   
   # Explore concurvity with mgcv function: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/concurvity.html
   # 1. Concurvity of each predictor compared to the whole model:
   concurvity(gam.model, full=TRUE)
   
   # 2. Concurvity for each pairwise comparison of predictor variables included in the model:
   concurvity(gam.model, full=FALSE)
   
   