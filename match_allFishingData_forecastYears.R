# match_allfishingdata_forecastYears.R
# Primary Author: Emily Norton (Dec 15 2020)
# emilyln@uw.edu
# Secondary Author: Isaac Kaplan 
# isaac.kaplan@noaa.gov

# This script is required for forecasting crab CPUE. This script matching ocean 
# conditions from the forecast year to ALL of the previous fishing behaviors
# observed. This must be run before "forecastGAM_fishingBehavSelect.R" can be 
# used to forecast the crab CPUE. 

#-------------------------------------------------------------------------------

# The following packages are required to run this script:
library(ncdf4)
library('tidyr')
library('geosphere')
library(tidyverse)
library('rgdal')
library('maptools')
 
setwd('~/Documents/DungenessCPUEModel/Rcode/')

forecastYear <- 2017    # Set forecastYear

#----------------
## Load the existing matrices with all J-SCOPE forecast + hindcast ocn data:
load(file="../Output/oceanMatrixAllMonthsYears.Rdata")
load(file="../Output/LagoceanMatrixAllYears.Rdata")  # LagoceanMatrixAllYears

#----------------
# Load prepped logbook data, with outliers removed, all available years (2007-2018), broken into "Historical" (2007-2015) and "Future" (2016-2018)
save(LogbookDataWithOceanConditionsHistorical, LogbookDataWithOceanConditionsFuture, file="../Output/WA_OR_LogbooksAndOceanConditions_HistoricalFuture_readyforGAM.RData")
AllLogbookDataWithOceanConditions <- rbind(LogbookDataWithOceanConditionsHistorical,LogbookDataWithOceanConditionsFuture)   

#----------------
# Winnow down the dataframe to the fishing obs and other important metadata
FishingObsColNames <- c("LnKgsPerPot","LbsPerPot","TotalLbs","Lons","Lats","SoakTime","julianDay","dayinseason","CrabYearS","state","distToNearestGridpointM","bathymetry","calendarMonthNum")
AllFishingObs <- AllLogbookDataWithOceanConditions[,FishingObsColNames]
names(AllFishingObs)[names(AllFishingObs) == 'LnLbsPerPot'] <- 'LnLbsPerPot_Obs'

#----------------
## Match all obs to J-SCOPE ocn forecasts for our forecast year:
source('MatchJSCOPEnetcdfToLogbookForecastYear.R')

forecastLogbookDataWithOceanConditions<- MatchJSCOPEnetcdfToLogbookForecastYear(oceanMatrixAllMonthsYears,AllFishingObs,forecastYear)

#save(forecastLogbookDataWithOceanConditions,file=paste("../Output/AllFishingObs_SubsamplingApproach/forecastLogbookDataWithOceanConditions_forecastYear", forecastYear,".Rdata",sep="")) 


#----------------
# Now match lagged ocean conditions from hindcasts:
source("MatchLaggedJSCOPEnetcdfToLogbook.R")

LagYrs <- c(3:4)    # consider 3- and 4-yr lags

for(y in 1:length(LagYrs))
{
  forecastLogbookDataWithOceanConditions <-  MatchLaggedJSCOPEnetcdfToLogbook(LagoceanMatrixAllYears,forecastLogbookDataWithOceanConditions, forecastYear, LagYrs[y])

  }

forecastLaggedLogbookDataWithOceanConditions <- forecastLogbookDataWithOceanConditions
#save(forecastLaggedLogbookDataWithOceanConditions,file=paste("../Output/AllFishingObs_SubsamplingApproach/forecastLaggedLogbookDataWithOceanConditions_forecastYear", forecastYear, ".Rdata",sep="")) 

#----------------
# Now match lagged (summed) PDO:
fPDO <- "../Data/PDO/pdo_JantoJuly_monthlySum_fromJISAOwebsite.csv"
PDOAllYears <- as.matrix(read.csv(fPDO))  

LagYrsPDO <- c(4)

source('MatchLaggedPDOToLogbook.R')  

for(y in 1:length(LagYrsPDO)) {
  forecastLaggedLogbookDataWithOceanConditions <-  MatchLaggedPDOToLogbook(PDOAllYears,forecastLaggedLogbookDataWithOceanConditions, forecastYear, LagYrsPDO[y])
}

forecastLaggedPDOLogbookDataWithOceanConditions <- forecastLaggedLogbookDataWithOceanConditions

#save(forecastLaggedPDOLogbookDataWithOceanConditions,file=paste("../Output/AllFishingObs_SubsamplingApproach/forecastLaggedPDOLogbookDataWithOceanConditions_forecastYear", forecastYear, ".Rdata",sep="")) 


# Done building our dataframe. Now start subsampling and making predictions using "forecastGAM_fishingBehavSelect.R"

