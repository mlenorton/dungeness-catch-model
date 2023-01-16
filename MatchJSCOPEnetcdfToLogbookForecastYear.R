MatchJSCOPEnetcdfToLogbookForecastYear <- function(oceanMatrixAllMonthsYears,LogbookDataWithOceanConditionsHistorical,startYearForROMSForecast)
{
#------
  # Primary Author: Isaac Kaplan (Feb 20 2019) 
  # isaac.kaplan at noaa.gov
  # Secondary Author: Emily Norton (Dec 2020)
  # emilyln@uw.edu
  
  # Matches existing logbook fishing behaviors [i.e., 'LogbookDataWithOceanConditionsHistorical'] 
  # (which has columns for oceanography) to a new single year of forecast oceanography. 
  # This matching (i.e., to a single year of oceanography), will retain historical logbook's 
  # fishing characteristics, but replace that hiostircal year's ocean characteristics per pot. 
  # Ocean matching based on calendar day (but with the intended mismatch of years.)
  # Note AFTER this function is run, we can resample from the output data frame to draw many replicates of 
  # ~1 years worth of sample points. 
  
  # Arguments:
  #       oceanMatrixAllMonthsYears: very large matrix of data for every grid point of ROMS Jscope NC files
  #       LogbookDataWithOceanConditionsHistorical: ALL YEARS of existing logbook data
  #       startYearForROMSForecast:  would be for instance 2017 for ROMS autum 2107--spring 2018.
  #                                   Important because this is used to find NC file we want to match to this logbook record.
 
  
  #   Returns: LogbookDataWithOceanConditionsForecast:  each row has OBSERVED logbook record (including catch rate per trap) 
                                         # but paired based on calendar day with associated ocean conditions for a future forecast month/year. 

  
  # Note (12/15/20) - grab the following columns from the Historical dataframe for the new forecast dataframe:
  # ForeColNames <- c("LbsPerPot","TotalLbs","Lons","Lats","SoakTime","julianDay","dayinseason","CrabYearS","state","distToNearestGridpointM","bathymetry")
#----------

  print("STARTING TO MATCH JSCOPE netcdf FILE TO CRAB LOGBOOK FOR 1 YEAR")

source("PlotCoastForSardine.r")

rowsForThisYear <- which(oceanMatrixAllMonthsYears[,"yearROMSstart"]==startYearForROMSForecast)

if(length(rowsForThisYear)==0)
{
  print("Error: rows of ocean Jscope data for this fishing year are 0. Check if year (and months) match up in Jscope versus fishery data")
  LogbookDataWithOceanConditionsForecast$bottomTemp[i] <- NaN   # Modified so the program can keep running even if it encounters an row that it doesn't have data for...pull these out at the end before making a GAM prediction 
  LogbookDataWithOceanConditionsForecast$bottomOxy[i] <- NaN
  LogbookDataWithOceanConditionsForecast$SSHVector[i] <- NaN
  LogbookDataWithOceanConditionsForecast$Chla2m[i] <- NaN
  LogbookDataWithOceanConditionsForecast$bottomArag[i] <- NaN
  LogbookDataWithOceanConditionsForecast$bottompH[i] <- NaN
  LogbookDataWithOceanConditionsForecast$bottomSalt[i] <- NaN
  # Note that "forecastYear" is the year when a forecast starts i.e. the fall when ROMS run starts, even if it is a spring logbook observation
  LogbookDataWithOceanConditionsForecast$forecastYear[i] <- startYearForROMSForecast
}
else {

oceanMatrixThisYear <- oceanMatrixAllMonthsYears[rowsForThisYear, ]


require(tidyr)

#-----------
# Note calendarMonthNum (along with the other fishing behavior characteristics) 
# already exists in "the"LogbookDataWithOceanConditionsHistorical" so we'll rely on that 
LogbookDataWithOceanConditionsForecast <- LogbookDataWithOceanConditionsHistorical

#---------------

#AIM NOW IS TO ADD COLUMNS OF BEST MATCHING ROMS JSCOPE OCEAN POINTS AND OCEAN VARIABLES (LAT, LON, AND COLUMNS FOR EACH VARIABLE)
# THIS IS ADDED to the crab logbook data as new columns, sort of like adding ocean observations to what the fisher saw in the trap.

require(geosphere)  # LOAD PACKAGE

for(i in 1:nrow(LogbookDataWithOceanConditionsForecast))
  {
 
  #-------------
   # For this logbook row's fishing month, subset oceanMatrix to only this month.
  calendarMonthNum <- LogbookDataWithOceanConditionsForecast$calendarMonthNum[i]
  rowsForThisMonth <- which(oceanMatrixThisYear[,"monthROMSOceanConditions"]==calendarMonthNum)

  # A quick check to make sure there's ROMS ocean data for this logbook point
   if(length(rowsForThisMonth)==0)
  {
    print("row")
    print(i)
    print("month")
    print(calendarMonthNum)
    print("rows of ocean Jscope data for this fishing month are 0. Check if year (and months) match up in Jscope versus fishery data")
    LogbookDataWithOceanConditionsForecast$bottomTemp[i] <- NaN    # Modified so the program can keep running even if it encounters an row that it doesn't have data for. Pull these out before making a GAM prediction.
    LogbookDataWithOceanConditionsForecast$bottomOxy[i] <- NaN
    LogbookDataWithOceanConditionsForecast$SSHVector[i] <- NaN
    LogbookDataWithOceanConditionsForecast$Chla2m[i] <- NaN
    LogbookDataWithOceanConditionsForecast$bottomArag[i] <- NaN
    LogbookDataWithOceanConditionsForecast$bottompH[i] <- NaN
    LogbookDataWithOceanConditionsForecast$bottomSalt[i] <- NaN
    # Note that "forecastYear" is the year when a forecast starts i.e. the fall when ROMS it starts, even if it is a spring logbook observation
    LogbookDataWithOceanConditionsForecast$forecastYear[i] <- startYearForROMSForecast
  }
  else {
  oceanMatrixThisYearThisMonth <- oceanMatrixThisYear[rowsForThisMonth, ]

  #----------


  ## reformat ocean lat,lon vectors into a matrix
oceanLonLatMatrix <- matrix(c(oceanMatrixThisYearThisMonth[,"nc_lonVector"], oceanMatrixThisYearThisMonth[ ,"nc_latVector"]), ncol = 2)

#  NOW FORMAT LON AND LAT
## subset lat,lon for single hake point
logbookLonLat <- c(LogbookDataWithOceanConditionsForecast$Lons[i],LogbookDataWithOceanConditionsForecast$Lats[i])

# THIS IS THE MAIN FUNCTION FROM MIKE MALICK TO GET DISTANCES BETWEEN LOGBOOK AND ROM POINTS
## calc dist b/w CRAB point and all ocean obs
dist.i <- geosphere::distHaversine(logbookLonLat, oceanLonLatMatrix)


## find index of ocean obs w/ min distance
ocean.ind <- which(dist.i == min(dist.i))[1]

## assign ocean obs w/ min distance to crab data

#-----------
# NOW ADD THE OCEAN VARIABLES AS NEW COLUMNS TO THE CRAB LOBGOOK DATA
#-----------
LogbookDataWithOceanConditionsForecast$bottomTemp[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomTempVector"]
LogbookDataWithOceanConditionsForecast$bottomOxy[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomOxyVector"]
LogbookDataWithOceanConditionsForecast$SSHVector[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "SSHVector"]
LogbookDataWithOceanConditionsForecast$Chla2m[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "Chla2mVector"]
LogbookDataWithOceanConditionsForecast$bottomArag[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomAragVector"]
LogbookDataWithOceanConditionsForecast$bottompH[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottompHVector"]
LogbookDataWithOceanConditionsForecast$bottomSalt[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomSaltVector"]

# Note that "forecastYear" is the year when a forecast starts i.e. the fall when ROMS it starts, even if it is a spring logbook observation
LogbookDataWithOceanConditionsForecast$forecastYear[i] <- startYearForROMSForecast

#-----------
#   JUST PRINT SOME INFO TO SCREEN PERIODICALLY SO YOU KNOW IT IS RUNNING (SLOWLY)
if(i%%1000==0)
{
print("row is")
print(i)
}
#---------------

}
}   # if months are included in ROMS
}   # if years are included in ROMS

return(LogbookDataWithOceanConditionsForecast)

}
