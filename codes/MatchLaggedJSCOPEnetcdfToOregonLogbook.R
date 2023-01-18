MatchLaggedJSCOPEnetcdfToOregonLogbook <- function(LagoceanMatrixAllYears,Oregon_Crab_Logbooks_Data_ThisYearOnly,startYearForFishingSeason,LagYrs)

{
  #------
  # Primary Author: Isaac Kaplan (Feb 20 2019) 
  # isaac.kaplan at noaa.gov
  # Secondary Author: Emily Norton
  # emilyln@uw.edu
  
  #  Matches the ocean points up to logbook locations, and assigns annually-averaged, lagged ocean variables as new columns of the logbook rows.
  #  Arguments:
  #       LagoceanMatrixAllYears: very large matrix of data for every grid point of ROMS Jscope NC files (annually-averaged, lagged conditions)
  #       Oregon_Crab_Logbooks_Data_ThisYearOnly: single year of Oregon logbook data
  #       startYearForFishingSeason:  would be for instance 2017 for fishing year autumn 2107--spring 2018.
  #                                   Important because this is used to find NC file we want to match to this logbook record.
  #       LagYrs: number of years by which to lag the annual J-SCOPE conditions compared to the startYearForFishingSeason
  #
  #   Returns: OregonThisYearLbsAndOceanConditions, each row has logbook record (including cath rate per trap) and associated ocean conditiosn predicted by ROMS nearest point for that year
  
  #-------------
  # NOTE THIS IS  A SLIGHLY MODIFIED VERSION OF MatchLaggedJSCOPEnetcdfToLogbook(  ), which was for Washington.
  # SEE THAT FILE FOR BETTER COMMENTS.
  # THIS TWEAKED VERSION IS NECESSARY BECAUSE OREGON LOGBOOK DATA LOOKS DIFFERENT FROM WASHINGTON LOGBOOKS.
  #-------------

print("STARTING TO MATCH JSCOPE netcdf FILE TO CRAB LOGBOOK FOR 1 YEAR")

  #--------------
source("PlotCoastForSardine.r")

rowsForThisYear <- which(LagoceanMatrixAllYears[,"yearROMSstart"]==startYearForFishingSeason-LagYrs)
if(length(rowsForThisYear)==0)
{
  stop("rows of ocean Jscope data for this fishing year are 0. Check if year (and months) match up in Jscope versus fishery data")
}

LagoceanMatrixThisYear <- LagoceanMatrixAllYears[rowsForThisYear, ]

require(tidyr)

OregonThisYearLbsAndOceanConditions <- Oregon_Crab_Logbooks_Data_ThisYearOnly
OregonThisYearLbsAndOceanConditions <- OregonThisYearLbsAndOceanConditions[which(!is.na(OregonThisYearLbsAndOceanConditions$AdjPounds)) ,]

# Standardize Lon and Lat to names and negative longitude(west) as used for Washington
OregonThisYearLbsAndOceanConditions$Lons<- -1*(OregonThisYearLbsAndOceanConditions$SetLongDecimal )
OregonThisYearLbsAndOceanConditions$Lats<-OregonThisYearLbsAndOceanConditions$SetLatDecimal

OregonThisYearLbsAndOceanConditions$LbsPerPot<-OregonThisYearLbsAndOceanConditions$AdjPounds/OregonThisYearLbsAndOceanConditions$NumPots

#-----------
#----------------------
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

for(i in 1:nrow(OregonThisYearLbsAndOceanConditions))  {
## reformat ocean lat,lon vectors into a matrix
LagoceanLonLatMatrix <- matrix(c(LagoceanMatrixThisYear[,"nc_lonVector"], LagoceanMatrixThisYear[ ,"nc_latVector"]), ncol = 2)

## subset lat,lon for single crab point
logbookLonLat <- c(OregonThisYearLbsAndOceanConditions$Lons[i],OregonThisYearLbsAndOceanConditions$Lats[i])

## calc dist b/w crab point and all ocean obs
dist.i <- geosphere::distHaversine(logbookLonLat, LagoceanLonLatMatrix)

## find index of ocean obs w/ min distance
ocean.ind <- which(dist.i == min(dist.i))[1]

## assign ocean obs w/ min distance to crab data
OregonThisYearLbsAndOceanConditions$LagbottomTemp[i] <- LagoceanMatrixThisYear[ ocean.ind, "bottomTempVector"]
OregonThisYearLbsAndOceanConditions$LagbottomOxy[i] <- LagoceanMatrixThisYear[ ocean.ind, "bottomOxyVector"]
OregonThisYearLbsAndOceanConditions$LagChla2m[i] <- LagoceanMatrixThisYear[ ocean.ind, "Chla2mVector"]


if(i%%100==0) {
print("row is")
print(i)
print("dist is")
print(OregonThisYearLbsAndOceanConditions$distToNearestGridpointM[i])
}


}

# Rename these columns based on how many year lag they have...will need to run this script multiple times for different year lags
names(OregonThisYearLbsAndOceanConditions)[names(OregonThisYearLbsAndOceanConditions) == 'LagbottomTemp'] <- paste0('bottomTemp_', LagYrs,'yrLag')
names(OregonThisYearLbsAndOceanConditions)[names(OregonThisYearLbsAndOceanConditions) == 'LagbottomOxy'] <- paste0('bottomOxy_', LagYrs,'yrLag')
names(OregonThisYearLbsAndOceanConditions)[names(OregonThisYearLbsAndOceanConditions) == 'LagChla2m'] <- paste0('Chla2m_', LagYrs,'yrLag')


return(OregonThisYearLbsAndOceanConditions)

}
