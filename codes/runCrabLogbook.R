##  runCrabLogbook.R
# Primary Author: Isaac Kaplan (May 1 2018)
# isaac.kaplan@noaa.gov
# Secondary Author: Emily Norton
# emilyln@uw.edu

# This is the parent function used to generate GAMs based on logbook data
# matched to J-SCOPE oceanography for Norton et al. (in prep)


#------------------

rm(list=ls())

#This program requires the following packages 
library(ncdf4)
library('tidyr')
library('geosphere')
library(tidyverse)
library('rgdal')
library('maptools')
library(stringr)
library(mgcv)
library(dplyr)
#library(PBSmapping)   # required by PlotCoastForSardine.R et al.
#library(lubridate)    # required by AssignDaysSinceStartDatePerArea.R
#library(stringr)      # required by AssignDaysSinceStartDatePerArea.R

# Users need to set to your working directory
setwd('~/Documents/DungenessCPUEModel/Rcode/')

SepForecasts <- TRUE   #Sep-initialized forecasts being used, which span Sep-May

# ------------Part I: Read in ocean conditions from J-SCOPE---------------------
#--------------
#   Get and count the forecast netcdf files
#--------------
# retrieve a list of [Sep-initialized forecast, monthly avg'd] nc files 
fdir <- "../JSCOPE_ocn_data_masked/SepForecasts/"
flist <- list.files(fdir, pattern = "^.*\\.(nc|NC|Nc|Nc)$")

print("length of flist is")
print(length(flist))

#--------------
#   Go through the NetCDF files, extract them,  and add them to a very large matrix (oceanMatrix) that has all ocean lat lon points in it, with corresponding ocean variables
#--------------
source("ExtractOceanMatrixFromOneMonthNCFile.R")
for (thisMonthYear in 1:length(flist)) {
  thisNCFileName <- paste0(fdir,flist[thisMonthYear])
  print(thisNCFileName)

  # Extract the NC data
  oceanMatrix<-ExtractOceanMatrixFromOneMonthNCFile(thisNCFileName)

  if (thisMonthYear ==1) {
     oceanMatrixAllMonthsYears <- oceanMatrix }
  else {
     oceanMatrixAllMonthsYears <-rbind(oceanMatrixAllMonthsYears, oceanMatrix)}
}

#--------------
# CLEAN UP LAND/MASKED POINTS and remove them from ocean Jscope matrix 
#--------------
rowsNotLandFore <- which(oceanMatrixAllMonthsYears[,"bathymetryVector"]>4)
oceanMatrixAllMonthsYears <- oceanMatrixAllMonthsYears[rowsNotLandFore,]

# Omit any rows with missing data, since this will make matching easier in the future: 
oceanMatrixAllMonthsYears <- na.omit(oceanMatrixAllMonthsYears)

save(oceanMatrixAllMonthsYears,file="../Output/oceanMatrixAllMonthsYears.Rdata")


#--------------
#   Now do the same thing, but for the lagged ocean conditions: go through the hindcast NetCDF files, extract them,  and add them all to a very large matrix (LaggedoceanMatrix) that has all ocean lat lon point in it, with corresponding ocean variables
#--------------
source("ExtractOceanMatrixFromOneYearNCFile.R")
require(stringr)
annavgFP <- "../JSCOPE_ocn_data_masked/Hindcasts/"
fLaglist <- list.files(annavgFP, pattern="^.*\\.(nc|NC|Nc|Nc)$")
for (thisYear in 1:length(fLaglist)) {
  thisLagNCFileName <- paste0(annavgFP,fLaglist[thisYear])    
  print(thisLagNCFileName)
  
  # Extracts the NC data
  LagoceanMatrix<-ExtractOceanMatrixFromOneYearNCFile(thisLagNCFileName)
  
  if (thisYear ==1) {
    LagoceanMatrixAllYears <- LagoceanMatrix
  } else {
    LagoceanMatrixAllYears <-rbind(LagoceanMatrixAllYears, LagoceanMatrix)
  }
}

# Omit any rows with missing data, since this will make matching easier in the future: 
LagoceanMatrixAllYears <- na.omit(LagoceanMatrixAllYears)

save(LagoceanMatrixAllYears,file="../Output/LagoceanMatrixAllYears.Rdata")


#--------------
#   Now do the same thing, but for the PDO values (summed Jan-July, from Shanks et al., 2010; Shanks, 2013): 
#--------------
source("ExtractOceanMatrixFromOneYearNCFile.R")
fPDO <- "../Data/PDO/pdo_JantoJuly_monthlySum_fromJISAOwebsite.csv"
PDOAllYears <- as.matrix(read.csv(fPDO))   


# ----------------Part II: Read in logbook data---------------------------------
#--------------
# Part IIa: READ IN WASHINGTON LOGBOOK CRAB DATA
#--------------
WAyrs <- c(2009:2018)  # loop through all the years for which we have Washington logbooks

for (wy in 1:length(WAyrs)) {
  WAyr <- WAyrs[wy]
  
fWAdata <- "../Data/Logbooks/WA/"

WDFWLogbookData <- read.csv(paste(fWAdata,"WDFW_Crab_Logbooks_Data_",WAyr,"Only.csv",sep=""),header=TRUE)

WDFWLogbookData$calendarMonthNum <- as.POSIXlt(as.Date(WDFWLogbookData$Set.Date,format='%m/%d/%Y'))$mo+1


#--------------
# Just use autumn logbook data since right now roms is just Autumn nc files
#--------------
#use dplyr to filter data
if (SepForecasts) {                             
  WDFWLogbookData <- WDFWLogbookData %>%    
    filter(calendarMonthNum > 10 | calendarMonthNum < 6)   # since the fishery doesn't open till November at the earliest, filter accordingly
} 


#--------------
# Match logbooks to J-SCOPE forecast oceanogrphy
#--------------
# i.e., match Logbook locations to corresponding Ocean Matrix grid locations and their ocean variables (note will add columns with  soak time, days since season, etc., below)
# This creates columns to logbook data with counts of crabs per pot [for each string] and matching ocean jscope conditions
# The third argument here should match the autumn that fishing is expected to start.

source("MatchJSCOPEnetcdfToLogbook.R")

WDFWCrabLbsAndOceanConditionsTemp <-  MatchJSCOPEnetcdfToLogbook(oceanMatrixAllMonthsYears,WDFWLogbookData, WAyr)


#--------------
# Match logbooks to lagged conditions (i.e., J-SCOPE hindcasts)
#--------------
# i.e., match Logbook locations to *Lagged* Ocean Matrix grid locations and their variables 
# This creates columns to logbook with counts of crabs per pot [for each string] and matching ocean jscope conditions
# The third argument here should match the autumn that fishing is expected to start. 
# Loop through the number of lagged years we're interested in (3-4)

LagYrs <- c(3,4)    # consider 3- and 4-yr lags

source("MatchLaggedJSCOPEnetcdfToLogbook.R")

for(y in 1:length(LagYrs)) {
WDFWCrabLbsAndOceanConditionsTemp <-  MatchLaggedJSCOPEnetcdfToLogbook(LagoceanMatrixAllYears,WDFWCrabLbsAndOceanConditionsTemp, WAyr, LagYrs[y])
}    


#--------------
## Now match lagged PDO to WA logbook
#--------------
LagYrsPDO <- c(4)   # consider only 4-yr lagged PDO (Shanks et al., 2010; Shanks, 2013)

source('MatchLaggedPDOToLogbook.R')
for(y in 1:length(LagYrsPDO)) {
  WDFWCrabLbsAndOceanConditionsTemp <-  MatchLaggedPDOToLogbook(PDOAllYears,WDFWCrabLbsAndOceanConditionsTemp, WAyr, LagYrsPDO[y])
}    


#--------------
# Add Opening Dates 
#--------------
# Assign Opening Dates per fishing area and add that column to the dataframe.
# i.e., allocate logbook entries to fishing zone, allocate start date per zone, 
# and then calculate DaysSinceStartDatePerArea
# Note: AssignDaysSinceStartDatePerArea(OceanConditionsMatrix, crabYear,forecastYear)
# where crabYear can be some previous year to inherit fishing characteristics such as soak time.But forecastYear is the actual ROMS oceanography year. 
# Forecast year is simply appended to the dataframe for later use or subsetting. 

source('point in spatial polygons function.R')   
source("AssignDaysSinceStartDatePerArea.R")
WDFWCrabLbsAndOceanConditions <-  AssignDaysSinceStartDatePerArea(WDFWCrabLbsAndOceanConditionsTemp, WAyr,WAyr)

# concatenate logbook/ocn condition data for all WA yrs
if (wy==1){
  WashingtonLogbooksAndOceanConditions <- WDFWCrabLbsAndOceanConditions
  } else {
    WashingtonLogbooksAndOceanConditions <- rbind(WashingtonLogbooksAndOceanConditions,WDFWCrabLbsAndOceanConditions) }

}    # end for loop through wy WAyrs


save(WashingtonLogbooksAndOceanConditions,file="../Output/WashingtonLogbooksAndOceanConditions.RData")

#----------------------------
# Part IIb: READ IN OREGON LOGBOOK CRAB DATA
#----------------------------

fORdata <- "../Data/Logbooks/OR/"

OregonLogbookData20072008to20152016 <- read.csv(paste(fORdata,"Kaplan_crablog_041318_051818.csv",sep=""),header=TRUE)
OregonLogbookData20162017to20172018 <- read.csv(paste(fORdata,"OR_CrabLogbooks_fromKellyCorbett_16to17_elnmod.csv",sep=""),header=TRUE)

OregonLogbookData20072008to20152016$calendarMonthNum<- as.POSIXlt(as.Date(OregonLogbookData20072008to20152016$DetailDate,format='%d-%b-%Y'))$mo+1
OregonLogbookData20162017to20172018$calendarMonthNum<- as.POSIXlt(as.Date(OregonLogbookData20162017to20172018$DetailDate,format='%d-%b-%Y'))$mo+1

#--------------
# Just use autumn logbook data since right now roms is just Autumn nc files
#--------------
#use dplyr to filter data
require(dplyr)   

if (SepForecasts) {

  OregonLogbookData20072008to20152016 <- OregonLogbookData20072008to20152016 %>%     # since the fishery doesn't open till November at the earliest, filter accordingly
    filter(calendarMonthNum > 10 | calendarMonthNum < 6)
  
  OregonLogbookData20162017to20172018 <- OregonLogbookData20162017to20172018 %>% 
    filter(calendarMonthNum > 10 | calendarMonthNum < 6)
} 

#--------------
# Prepare to match logbook data to ocn
#--------------

ORyrs <- c(2007:2017)  # loop through all the years for which we have Oregon logbooks

for (oy in 1:length(ORyrs)) {
  ORyr <- ORyrs[oy]
  
#--------------
# Pull out individual years of Oregon data (to make matching to ocn easier)
#--------------

  if (ORyr < 2016) {
    rowsThisYr <- which ( OregonLogbookData20072008to20152016$CrabYear== paste(ORyr,"-",ORyr+1,sep=""))
    ORLogbookData<-OregonLogbookData20072008to20152016[rowsThisYr , ]
  } else {
    rowsThisYr <- which ( OregonLogbookData20162017to20172018$CrabYear== paste(ORyr,"-",ORyr+1,sep=""))
    ORLogbookData<-OregonLogbookData20162017to20172018[rowsThisYr , ]
  }

#--------------
# Match logbooks to J-SCOPE forecast oceanogrphy
#--------------
# i.e., match Logbook locations to corresponding Ocean Matrix grid locations and their ocean variables (note will add columns with  soak time, days since season, etc., below)
# This creates columns to logbook data with counts of crabs per pot [for each string] and matching ocean jscope conditions
# The third argument here should match the autumn that fishing is expected to start.
# Note different Matching function used here because Oregon logbook formats differ from Washington
source("MatchJSCOPEnetcdfToOregonLogbook.R")
# note: MatchJSCOPEnetcdfToOregonLogbook <- function(oceanMatrix,logbookData,startYearForFishingSeason,startYearForOceanography,forecastYear
# This lets you for instance use fishingSeason fishing characteristics from 2015, ocean from 2012, but really intend to forecast 2020 (if you wanted to)
#     forecastYear is really just appended to the dataframe for later use.  

ORLogbooksAndOceanConditions <-  MatchJSCOPEnetcdfToOregonLogbook(oceanMatrixAllMonthsYears,ORLogbookData, ORyr,ORyr,ORyr)


#--------------
# Match logbooks to lagged conditions (i.e., J-SCOPE hindcasts)
#--------------
# i.e., match Logbook locations to *Lagged* Ocean Matrix grid locations and their variables 
# This creates columns to logbook with counts of crabs per pot [for each string] and matching ocean jscope conditions
# The third argument here should match the autumn that fishing is expected to start. 
# Loop through the number of lagged years we're interested in (3-4)

LagYrs <- c(3:4)        # consider 3- and 4-yr lagged conditions

source("MatchLaggedJSCOPEnetcdfToOregonLogbook.R")

for(y in 1:length(LagYrs)) {
  ORLogbooksAndOceanConditions <-  MatchLaggedJSCOPEnetcdfToOregonLogbook(LagoceanMatrixAllYears,ORLogbooksAndOceanConditions, ORyr, LagYrs[y])
}

#--------------
## Now match lagged PDO to OR logbook
#--------------
LagYrsPDO <- c(4)   # consider only 4-yr lagged PDO (Shanks et al., 2010; Shanks, 2013)

source('MatchLaggedPDOToOregonLogbook.R')

for(y in 1:length(LagYrsPDO)) {
  ORLogbooksAndOceanConditions <-  MatchLaggedPDOToOregonLogbook(PDOAllYears,ORLogbooksAndOceanConditions, ORyr, LagYrsPDO[y])
  graphics.off()
}

# concatenate logbook/ocn condition data for all OR yrs
if (oy==1){
  OregonLogbooksAndOceanConditions <- ORLogbooksAndOceanConditions
} else {
  OregonLogbooksAndOceanConditions <- rbind(OregonLogbooksAndOceanConditions,ORLogbooksAndOceanConditions)}
}    # end for loop through oy ORyrs


save(OregonLogbooksAndOceanConditions,file="../Output/OregonLogbooksAndOceanConditions.RData")


# ----------Part III: Add add'l fields, QA/QC, and combine OR + WA data---------
#--------------
# Convert to day of year
#--------------
WashingtonLogbooksAndOceanConditions$julianDay <- as.POSIXlt(WashingtonLogbooksAndOceanConditions$Set.Date, format = "%m/%d/%y")$yday
OregonLogbooksAndOceanConditions$julianDay <- as.POSIXlt(OregonLogbooksAndOceanConditions$DetailDate, format = "%d-%b-%y")$yday

#--------------
# Match points to grain size
#--------------
source('MatchEnvironmentLayerToLatLon.R')

WashingtonLogbooksAndOceanConditions <-  MatchEnvironmentLayerToLatLon(WashingtonLogbooksAndOceanConditions,"../Data/grainsize_pts_geo.csv","GRAINSIZE")
OregonLogbooksAndOceanConditions <-  MatchEnvironmentLayerToLatLon(OregonLogbooksAndOceanConditions,"../Data/grainsize_pts_geo.csv","GRAINSIZE")

# & rename some columns, including those that we just added:
names(WashingtonLogbooksAndOceanConditions)[names(WashingtonLogbooksAndOceanConditions) == 'Soak.Time..days.'] <- 'SoakTime'
names(WashingtonLogbooksAndOceanConditions)[names(WashingtonLogbooksAndOceanConditions) == 'NewHabLayer'] <- 'GrainSize'       
names(WashingtonLogbooksAndOceanConditions)[names(WashingtonLogbooksAndOceanConditions) == 'crabYear'] <- 'CrabYearS'
names(WashingtonLogbooksAndOceanConditions)[names(WashingtonLogbooksAndOceanConditions) == 'Crab.Retained..lbs.'] <- 'TotalLbs'
names(WashingtonLogbooksAndOceanConditions)[names(WashingtonLogbooksAndOceanConditions) == 'Pots.Fished'] <- 'NumPotsFished'

names(OregonLogbooksAndOceanConditions)[names(OregonLogbooksAndOceanConditions) == 'NewHabLayer'] <- 'GrainSize'      
names(OregonLogbooksAndOceanConditions)[names(OregonLogbooksAndOceanConditions) == 'AdjPounds'] <- 'TotalLbs'
names(OregonLogbooksAndOceanConditions)[names(OregonLogbooksAndOceanConditions) == 'NumPots'] <- 'NumPotsFished'

#--------------
# Further standardize OR + WA data so they can be combined
#--------------

# For OR data only, screen out points that have been flagged for a spatial error 

OregonLogbooksAndOceanConditions <- OregonLogbooksAndOceanConditions[which(OregonLogbooksAndOceanConditions$Spatial.Error==FALSE), ]

#Make CrabYearStandardized that can be used for both States
# WASHINGTON:
WashingtonLogbooksAndOceanConditions <- WashingtonLogbooksAndOceanConditions %>% 
  mutate(CrabYearStandardized = CrabYearS)

# OREGON: (note - CrabYear was for instance 2007-2008, this needs to be turned into characters, and split at the hyphen into two columsn (we only wan the first column))
OregonLogbooksAndOceanConditions$CrabYearStandardized <- data.frame(do.call('rbind', strsplit(as.character(OregonLogbooksAndOceanConditions$CrabYear),'-',fixed=TRUE)))[,1]
# Make CrabYearStandardized a character string, as it is for Washington: 
OregonLogbooksAndOceanConditions$CrabYearStandardized <- as.character(OregonLogbooksAndOceanConditions$CrabYearStandardized)


require(dplyr)
library(tidyverse)
OregonLogbooksAndOceanConditions <- OregonLogbooksAndOceanConditions %>%
  mutate(CrabYearS = as.numeric(CrabYearStandardized)) %>% 
  add_column(state = 'OR')

WashingtonLogbooksAndOceanConditions <- WashingtonLogbooksAndOceanConditions %>% 
  mutate(CrabYearS = as.numeric(CrabYearStandardized)) %>% 
  mutate(LagPDO = as.numeric(PDO_4yrLag)) %>%
  add_column(state = 'WA')



#------------
# Combine OR + WA data
#------------
commonColNames <- c("NumPotsFished","LbsPerPot","TotalLbs","Lons","Lats","bathymetry","bottomTemp","bottomOxy","SSHVector","Chla2m","bottomArag","bottompH","bottomSalt","GrainSize","SoakTime","julianDay","dayinseason","distToNearestGridpointM","CrabYearS","calendarMonthNum","PDO_4yrLag","bottomOxy_4yrLag","bottomTemp_4yrLag","Chla2m_4yrLag","Chla2m_3yrLag","bottomOxy_3yrLag","bottomTemp_3yrLag","forecastYear","state")  #Include chlorophyll lag variables too!

LogbookDataWithOceanConditions<-rbind(OregonLogbooksAndOceanConditions[,commonColNames], WashingtonLogbooksAndOceanConditions[,commonColNames])
save(LogbookDataWithOceanConditions, file = "../Output/WA_OR_LogbooksAndOceanConditions.RData")

#------------
# Deal with lbs vs kg and log-transform
#------------
LogbookDataWithOceanConditions$LnLbsPerPot <- log(LogbookDataWithOceanConditions$LbsPerPot+0.01)
LogbookDataWithOceanConditions$LnTotalLbs <- log(LogbookDataWithOceanConditions$TotalLbs+0.01)

#Convert Lbs to Kgs, calculate log(KgsPerPot), and combine with existing matrix
KgsPerPot <-LogbookDataWithOceanConditions$LbsPerPot*0.45359   #Convert lbs to kgs
LnKgsPerPot <-log(KgsPerPot+0.01)
LogbookDataWithOceanConditions<-cbind(LogbookDataWithOceanConditions,KgsPerPot,LnKgsPerPot)
TotalKgs <-LogbookDataWithOceanConditions$TotalLbs*0.45359   #Convert lbs to kgs
LnTotalKgs <-log(TotalKgs+0.01)
LogbookDataWithOceanConditions<-cbind(LogbookDataWithOceanConditions,TotalKgs,LnTotalKgs)

#------------
# Deal with NaN/0/inf. problems: 
#------------
# in logbook data: 
allNonNANonInfRows <- which(!is.na(LogbookDataWithOceanConditions$LnLbsPerPot) & !is.infinite(LogbookDataWithOceanConditions$LnLbsPerPot) ) 
LogbookDataWithOceanConditions <- LogbookDataWithOceanConditions [allNonNANonInfRows ,  ]   

# in ocean data: 
NaNOceanRows <- which(LogbookDataWithOceanConditions$bottomTemp == "NaN"  ) 
NaNLogbookDataWithOceanConditions <- LogbookDataWithOceanConditions [NaNOceanRows ,  ]

#------------
# Other data QA/QC
#------------
# Cutting out all bathmetry =4 (land) logbook sets
LogbookDataWithOceanConditions <-LogbookDataWithOceanConditions[LogbookDataWithOceanConditions$bathymetry>4,]

# CLEAN UP LAND POINTS THAT ARE POTS TOO FAR FROM OCEAN ROMS POINTS
LogbookDataWithOceanConditions <-LogbookDataWithOceanConditions[LogbookDataWithOceanConditions$distToNearestGridpointM<2000,]   

#------------
#Check for lingering outliers and remove them before fitting gams
#------------
LogbookDataWithOceanConditionsWithOutliers <- LogbookDataWithOceanConditions

#Get rid of outliers in catches (>50 kgs) 
goodi_kgs <- which(LogbookDataWithOceanConditionsWithOutliers$KgsPerPot<=50,)  

#Get rid of outliers in bathymetry (>220 m)
goodi_bathy <- which(LogbookDataWithOceanConditionsWithOutliers$bathymetry<=220,)  

#Get rid of outliers in soaktime (>30 days)
goodi_soak <- which(LogbookDataWithOceanConditionsWithOutliers$SoakTime<=30)  

#cut out day in season < 0 or > 365
goodi_dis1 <- which(LogbookDataWithOceanConditionsWithOutliers$dayinseason>=0)
goodi_dis2 <- which(LogbookDataWithOceanConditionsWithOutliers$dayinseason<367)
goodi_dis <- intersect(goodi_dis1,goodi_dis2)

#Get rid of bad aragonite values (due to algorithm)
goodi_arag <- which(LogbookDataWithOceanConditionsWithOutliers$bottomArag>0)

# find which inds have high Julian Day and high day in season (we'll remove these)
goodi_JD1 <- which(LogbookDataWithOceanConditionsWithOutliers$julianDay<=200)   # these are good, they fall before the June cut-off for the Sep forecast
questi_JD2 <- which(LogbookDataWithOceanConditionsWithOutliers$julianDay>200)
goodi_JD2 <- which(LogbookDataWithOceanConditionsWithOutliers$dayinseason[questi_JD2]<=100)
goodi_JD <- union(goodi_JD1, questi_JD2[goodi_JD2])

goodinds <- Reduce(intersect, list(goodi_kgs,goodi_bathy,goodi_soak,goodi_dis,goodi_arag,goodi_JD))

LogbookDataWithOceanConditionsNoOUT <- LogbookDataWithOceanConditionsWithOutliers[goodinds,]
LogbookDataWithOceanConditions <- LogbookDataWithOceanConditionsNoOUT

# Omit any rows with missing data, since this will make prediction hard: 
LogbookDataWithOceanConditions <- na.omit(LogbookDataWithOceanConditions)


#------------
# Subset "historical" (i.e. model training) vs "future" (i.e. model testing) years and save
#------------

# Just note the historical data that will be used to fit the GAM: 
LogbookDataWithOceanConditionsHistorical <- LogbookDataWithOceanConditions[which(LogbookDataWithOceanConditions$forecastYear<2016) ,]

# This is not used below, but just keeping track of the true forecast years. 
LogbookDataWithOceanConditionsFuture <- LogbookDataWithOceanConditions[which(LogbookDataWithOceanConditions$forecastYear>2015) ,]

save(LogbookDataWithOceanConditionsHistorical, LogbookDataWithOceanConditionsFuture, file="../Output/WA_OR_LogbooksAndOceanConditions_HistoricalFuture_readyforGAM.RData")


# -----------------------Part IV: Fit GAMs--------------------------------------
#--------------
## Fit univariate gams
#--------------
require(mgcv)

gamU <- gam(formula= LnKgsPerPot ~  s(dayinseason,k=3), family=gaussian, data = LogbookDataWithOceanConditionsHistorical, scale=-1,select=TRUE)   #univariate gams


#--------------
## Fit multivariate gams
#--------------
gamDiS <- gam(formula= LnKgsPerPot ~  s(dayinseason,k=3), family=gaussian, data = LogbookDataWithOceanConditionsHistorical, scale=-1,select=TRUE)   # day in season only
gamDiS_D <- gam(formula= LnKgsPerPot ~  s(dayinseason,k=3)+s(SSHVector,k=3)+s(bottomArag,k=3)+s(bottompH,k=3)+s(bottomSalt,k=3)+s(bottomTemp,k=3)+s(bottomOxy,k=3)+s(Chla2m,k=3), family=gaussian, data = LogbookDataWithOceanConditionsHistorical, scale=-1,select=TRUE)    # day in season + dynamic vars
gamDiS_L <- gam(formula= LnKgsPerPot ~  s(dayinseason,k=3)+s(bottomTemp_4yrLag,k=3)+s(bottomOxy_4yrLag,k=3)+s(Chla2m_4yrLag,k=3)+s(bottomTemp_3yrLag,k=3)+s(bottomOxy_3yrLag,k=3)+s(Chla2m_3yrLag,k=3)+s(PDO_4yrLag,k=3), family=gaussian, data = LogbookDataWithOceanConditionsHistorical, scale=-1,select=TRUE)   # day in season + dynamic + lagged vars
gamDiS_DL <- gam(formula= LnKgsPerPot ~  s(dayinseason,k=3)+s(SSHVector,k=3)+s(bottomArag,k=3)+s(bottompH,k=3)+s(bottomSalt,k=3)+s(bottomTemp,k=3)+s(bottomOxy,k=3)+s(Chla2m,k=3)+s(bottomTemp_4yrLag,k=3)+s(bottomOxy_4yrLag,k=3)+s(Chla2m_4yrLag,k=3)+s(bottomTemp_3yrLag,k=3)+s(bottomOxy_3yrLag,k=3)+s(Chla2m_3yrLag,k=3)+s(PDO_4yrLag,k=3), family=gaussian, data = LogbookDataWithOceanConditionsHistorical, scale=-1,select=TRUE)   # day in season + dynamic + lagged vars
gamDiS_SDL <- gam(formula= LnKgsPerPot ~  s(dayinseason,k=3)+s(bathymetry,k=3)+s(SoakTime,k=3)+s(CrabYearS,k=1)+te(Lats,Lons,k=3)+s(SSHVector,k=3)+s(bottomArag,k=3)+s(bottompH,k=3)+s(bottomSalt,k=3)+s(bottomTemp,k=3)+s(bottomOxy,k=3)+s(Chla2m,k=3)+s(bottomTemp_4yrLag,k=3)+s(bottomOxy_4yrLag,k=3)+s(Chla2m_4yrLag,k=3)+s(bottomTemp_3yrLag,k=3)+s(bottomOxy_3yrLag,k=3)+s(Chla2m_3yrLag,k=3)+s(PDO_4yrLag,k=3), family=gaussian, data = LogbookDataWithOceanConditionsHistorical, scale=-1,select=TRUE)   # day in season + static + dynamic + lagged vars
gamDiS_SDL_less <- gam(formula= LnKgsPerPot ~  s(dayinseason,k=3)+s(bathymetry,k=3)+s(SoakTime,k=3)+s(CrabYearS,k=9)+te(Lats,Lons,k=3)+s(bottomArag,k=3)+s(bottomSalt,k=3)+s(bottomTemp,k=3)+s(bottomOxy,k=3)+s(Chla2m,k=3)+s(bottomTemp_4yrLag,k=3)+s(bottomOxy_4yrLag,k=3)+s(Chla2m_4yrLag,k=3)+s(bottomTemp_3yrLag,k=3)+s(bottomOxy_3yrLag,k=3), family=gaussian, data = LogbookDataWithOceanConditionsHistorical, scale=-1,select=TRUE)   # day in season + static + dynamic + lagged vars (simplified, but with same AIC)


bestgam <- gamDiS_SDL_less

save(gamDiS_SDL_less, file="../Output/gamDiS_SDL_less.RData")


###--------------- Part V: Explore GAM & Make Figures --------------------------
summary(bestgam)
AIC(bestgam)
bestgam$gcv
anova(bestgam)

dev.new()
gam.check(bestgam)


# Plot subplots for each explanatory variable from the GAM 
dev.new()
#tiff("../Manuscript/Tables_Figures/GAM_SDL_BEST_101620.tiff", width=11, height=8.5, units='in', res=300)
par(mfrow=c(3,6),
    oma = c(5,1,1,0) + 0.1,
    mar = c(3,2,1,1) + 0.1)
plot.gam(bestgam, residuals=F, pch=16, all.terms=T, rugplot=T, se = T, rug=T, scale=0, shade=T)   #option 1: scaled ylim to fit var + shaded SE
#plot.gam(bestgam, residuals=F, pch=16, all.terms=T, rugplot=T, se = T, rug=T)        #option 2: consistent ylim; dashed SE lines
dev.off()
