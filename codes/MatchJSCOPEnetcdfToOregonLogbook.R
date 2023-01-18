MatchJSCOPEnetcdfToOregonLogbook <- function(oceanMatrixAllMonthsYears,Oregon_Crab_Logbooks_Data_ThisYearOnly,startYearForFishingSeason,startYearForOceanography,forecastYear)
{
  #----------
  # Author: Isaac Kaplan (Feb 20 2019) 
  # isaac.kaplan at noaa.gov
  
  #  Matches the ocean points up to logbook locations, and assign ocean variables as new columns of the logbook rows.
  #  Arguments:
  #       oceanMatrixAllMonthsYears: very large matrix of data for every grid point of ROMS Jscope NC files
  #       Oregon_Crab_Logbooks_Data_ThisYearOnly: single year of Oregon logbook data
  #       startYearForFishingSeason:  would be for instance 2017 for fishing year autumn 2017--spring 2018.
  #                                   Important because this is used to find NC file we want to match to this logbook record.
  #       startYearForOceanography: the (autumn) year of modeled ocean conditions that you'd like to match the fishing behavior to 
  #       forecastYear: adds a new column indicating which year you'd like to use these conditions when forecasting (usually same as startYearForOceanography, but could be different)
  #  Returns: WDFWThisYearLbsAndCountsConvertedToLbs,  each row has logbook record (including cath rate per trap) and associated ocean conditions predicted by ROMS nearest point for that month.
  
  # NOTE THIS IS A SLIGHLY MODIFIED VERSION OF MatchJSCOPEnetcdfToLogbook(  ), which was for Washington.
  # SEE THAT FILE FOR BETTER COMMENTS.
  # THIS TWEAKED VERSION IS NECESSARY BECAUSE OREGON LOGBOOK DATA LOOKS DIFFERENT FROM WASHINGTON LOGBOOKS.
  #-------------

print("STARTING TO MATCH JSCOPE netcdf FILE TO CRAB LOGBOOK FOR 1 YEAR")

  #--------------
source("PlotCoastForSardine.r")
#setwd("C:/Users/Isaac.Kaplan/Documents/FATE/CrabModelingJSCOPE/Rcode")

rowsForThisYear <- which(oceanMatrixAllMonthsYears[,"yearROMSstart"]==startYearForOceanography)
if(length(rowsForThisYear)==0)
{
  stop("rows of ocean Jscope data for this fishing year are 0. Check if year (and months) match up in Jscope versus fishery data")
}

oceanMatrixThisYear <- oceanMatrixAllMonthsYears[rowsForThisYear, ]

require(tidyr)

OregonThisYearLbsAndOceanConditions <- Oregon_Crab_Logbooks_Data_ThisYearOnly
OregonThisYearLbsAndOceanConditions <- OregonThisYearLbsAndOceanConditions[which(!is.na(OregonThisYearLbsAndOceanConditions$AdjPounds)) ,]

# Standardize Lon and Lat to names and negative longitude(west) as used for Washington
OregonThisYearLbsAndOceanConditions$Lons<- -1*(OregonThisYearLbsAndOceanConditions$SetLongDecimal )
OregonThisYearLbsAndOceanConditions$Lats<-OregonThisYearLbsAndOceanConditions$SetLatDecimal

OregonThisYearLbsAndOceanConditions$LbsPerPot<-OregonThisYearLbsAndOceanConditions$AdjPounds/OregonThisYearLbsAndOceanConditions$NumPots

#----------------
# REMOVE ROWS THAT LACK LAT OR LON:  these NAs cause geosphere::distHaversine to fail.
rowsWithLonLat <- which(!is.na(OregonThisYearLbsAndOceanConditions$Lons) & !is.na(OregonThisYearLbsAndOceanConditions$Lats) )
OregonThisYearLbsAndOceanConditions <- OregonThisYearLbsAndOceanConditions[rowsWithLonLat, ]
#----------------

# Now need to add column of days since start of fishing season for that port
# OPENERS
# ---------------
SeasonOpeningsByPort <- read.csv(file="../Data/SeasonOpenings/OR_SeasonOpeningByPort5copies2.csv",header=TRUE,as.is=FALSE )
SeasonOpeningsByPort$OpenerDate <- as.POSIXct(as.Date(SeasonOpeningsByPort$StartDateMod,format='%m/%d/%Y'))

#----------------
logbookPlotName<-paste(startYearForFishingSeason,"Oregon Logbook Lbs")

PlotCoastForSardine(thisPlotName = logbookPlotName,myXlim=c(-126.0, -123.0), myYlim=c(42.0, 46.5))  # c(-126.0, -123.0),myYlim=c(42.5, 50.5)
title(main= logbookPlotName,cex.main=0.8)
points(OregonThisYearLbsAndOceanConditions$Lons,OregonThisYearLbsAndOceanConditions$Lats, pch=15,cex=0.2)
dev.copy(png,logbookPlotName)
dev.off()

#----------------
require(geosphere)
#AIM NOW IS TO ADD COLUMNS OF BEST MATCHING ROMS JSCOPE POINTS (LAT, LON, AND COLUMNS FOR EACH VARIABLE)

for(i in 1:nrow(OregonThisYearLbsAndOceanConditions)) {
  # For this logbook row's fishing month, subset oceanMatrix to only this month.
    #         for WA: calendarMonthNum <- as.POSIXlt(as.Date(WDFWLogbookData2009$Set.Date[i],format='%m/%d/%Y'))$mo+1
   calendarMonthNum<- OregonThisYearLbsAndOceanConditions$calendarMonthNum[i]

   rowsForThisMonth <- which(oceanMatrixThisYear[,"monthROMSOceanConditions"]==calendarMonthNum)
  if(length(rowsForThisMonth)==0)
  {
    print("OregonThisYearLbsAndOceanConditions$DetailDate[i]")
    print(OregonThisYearLbsAndOceanConditions$DetailDate[i])
    print("calendarMonthNum is:")
    print(calendarMonthNum)
    stop("rows of ocean Jscope data for this fishing month are 0. Check if year (and months) match up in Jscope versus fishery data")
  }
  oceanMatrixThisYearThisMonth <- oceanMatrixThisYear[rowsForThisMonth, ]
  #----------------

## reformat ocean lat,lon vectors into a matrix
oceanLonLatMatrix <- matrix(c(oceanMatrixThisYearThisMonth[,"nc_lonVector"], oceanMatrixThisYearThisMonth[ ,"nc_latVector"]), ncol = 2)

## subset lat,lon for single crab point
logbookLonLat <- c(OregonThisYearLbsAndOceanConditions$Lons[i],OregonThisYearLbsAndOceanConditions$Lats[i])

## calc dist b/w crab point and all ocean obs
dist.i <- geosphere::distHaversine(logbookLonLat, oceanLonLatMatrix)

## find index of ocean obs w/ min distance
ocean.ind <- which(dist.i == min(dist.i))[1]

## assign ocean obs w/ min distance to hake data
OregonThisYearLbsAndOceanConditions$distToNearestGridpointM[i] <- dist.i[ocean.ind]
OregonThisYearLbsAndOceanConditions$nc_lat[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "nc_latVector"]
OregonThisYearLbsAndOceanConditions$nc_lon[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "nc_lonVector"]
OregonThisYearLbsAndOceanConditions$bathymetry[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bathymetryVector"]
OregonThisYearLbsAndOceanConditions$bottomTemp[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomTempVector"]
OregonThisYearLbsAndOceanConditions$bottomOxy[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomOxyVector"]
OregonThisYearLbsAndOceanConditions$SSHVector[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "SSHVector"]
OregonThisYearLbsAndOceanConditions$Chla2m[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "Chla2mVector"]
OregonThisYearLbsAndOceanConditions$bottomArag[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomAragVector"]
OregonThisYearLbsAndOceanConditions$bottompH[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottompHVector"]
OregonThisYearLbsAndOceanConditions$bottomSalt[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomSaltVector"]
OregonThisYearLbsAndOceanConditions$forecastYear[i]<-forecastYear
  
# OPENERS
dateIndex <- which((SeasonOpeningsByPort$Port==OregonThisYearLbsAndOceanConditions$PortCode[i]) & (SeasonOpeningsByPort$CrabYear== startYearForFishingSeason))
OregonThisYearLbsAndOceanConditions$SeasonOpeningsByPort[i]  <-SeasonOpeningsByPort$OpenerDate[dateIndex]


if(i%%100==0) {
print("row is")
print(i)
print("dist is")
print(OregonThisYearLbsAndOceanConditions$distToNearestGridpointM[i])
}


}

OregonThisYearLbsAndOceanConditions$DetailDateNum <- as.numeric(as.POSIXct(as.Date(OregonThisYearLbsAndOceanConditions$DetailDate,format='%d-%b-%y')))
OregonThisYearLbsAndOceanConditions$dayinseason <-(OregonThisYearLbsAndOceanConditions$DetailDateNum-OregonThisYearLbsAndOceanConditions$SeasonOpeningsByPort)/86400


return(OregonThisYearLbsAndOceanConditions)

}
