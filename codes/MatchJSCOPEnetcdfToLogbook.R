MatchJSCOPEnetcdfToLogbook <- function(oceanMatrixAllMonthsYears,WDFW_Crab_Logbooks_Data_ThisYearOnly,startYearForFishingSeason)
{
#----------
  # Author: Isaac Kaplan (Feb 20 2019) 
  # isaac.kaplan at noaa.gov
  
  #  Matches the ocean points up to logbook locations, and assign ocean variables as new columns of the logbook rows.
  #  Arguments:
  #       oceanMatrixAllMonthsYears: very large matrix of data for every grid point of ROMS Jscope NC files
  #       WDFW_Crab_Logbooks_Data_ThisYearOnly: single year of Washington logbook data
  #       startYearForFishingSeason:  would be for instance 2017 for fishing year autumn 2017--spring 2018.
  #                                   Important because this is used to find NC file we want to match to this logbook record.
  #
  #  Returns: WDFWThisYearLbsAndCountsConvertedToLbs,  each row has logbook record (including cath rate per trap) and associated ocean conditions predicted by ROMS nearest point for that month.

#----------
print("STARTING TO MATCH JSCOPE NETCDF FILE TO CRAB LOGBOOK FOR 1 YEAR")

source("PlotCoastForSardine.r")

rowsForThisYear <- which(oceanMatrixAllMonthsYears[,"yearROMSstart"]==startYearForFishingSeason)

if(length(rowsForThisYear)==0)
{
  stop("Error: rows of ocean Jscope data for this fishing year are 0. Check if year (and months) match up in Jscope versus fishery data")
}

oceanMatrixThisYear <- oceanMatrixAllMonthsYears[rowsForThisYear, ]

require(tidyr)

#----------
# Note that Washington Data has Lbs and counts.
# Ultimately we want lbs, but we will look at each of these separately,
# then combine them below.
#-----------------
WDFWThisYearCounts <- WDFW_Crab_Logbooks_Data_ThisYearOnly
WDFWThisYearLbs <- WDFW_Crab_Logbooks_Data_ThisYearOnly

#----------------
# To avoid duplication of data rows if there were cases with both counts and Lbs

indicesHasLbs<-which(!is.na(WDFWThisYearLbs$Crab.Retained..lbs.))
indicesLacksLbs<-which(is.na(WDFWThisYearLbs$Crab.Retained..lbs.))

WDFWThisYearCounts<- WDFWThisYearCounts[indicesLacksLbs ,]
WDFWThisYearLbs <- WDFWThisYearLbs[ indicesHasLbs,]
#------------------

WDFWThisYearCounts[,"CountsPerPot"] <- NA
WDFWThisYearLbs[,"CountsPerPot"] <- NA
WDFWThisYearCounts[,"LbsPerPot"] <- NA
WDFWThisYearLbs[,"LbsPerPot"] <- NA


#-------------
# CONSIDER THE COUNTS PER POT LOCATIONS
# But also convert to Pounds since that is what is used in Oregon and (kg) will be used for GAM.

LonsCounts<- -1*(WDFWThisYearCounts$Longitude.Begin.Degrees + WDFWThisYearCounts$Longitude.Begin.Minutes/60 )
LatsCounts<-WDFWThisYearCounts$Latitude.Begin.Degrees+ WDFWThisYearCounts$Latitude.Begin.Minutes/60

# CONVERT FROM COUNTS TO LBS FOR HAULS WHERE THIS IS NECESSARY. 
avgWtPerCrabLbs <- 1.8  # 1.8 from Carol Henry WDFW per comm May 4 2018

WDFWThisYearCounts$Crab.Retained..lbs.<-WDFWThisYearCounts$Crab.Retained..count.*avgWtPerCrabLbs

#--------------
# PLOT THE LBS PER POT LOCATIONS

LonsLbs<- -1*(WDFWThisYearLbs$Longitude.Begin.Degrees + WDFWThisYearLbs$Longitude.Begin.Minutes/60 )
LatsLbs<-WDFWThisYearLbs$Latitude.Begin.Degrees+ WDFWThisYearLbs$Latitude.Begin.Minutes/60

#PlotCoastForSardine()
#title(main= "WDFWThisYear Lbs PerPot",cex.main=1.5)
#points(LonsLbs,LatsLbs, pch=15,cex=0.2)


#--------------------
# NOW COMBINE THE LOGBOOK RECORDS that recorded catch in lbs with the records that recorded catch in counts.
#-----------------

WDFWThisYearLbsAndCountsConvertedToLbs <-rbind(WDFWThisYearLbs, WDFWThisYearCounts)

WDFWThisYearLbsAndCountsConvertedToLbs$CountsPerPot<-WDFWThisYearLbsAndCountsConvertedToLbs$Crab.Retained..count./WDFWThisYearLbsAndCountsConvertedToLbs$Pots.Fished

WDFWThisYearLbsAndCountsConvertedToLbs$LbsPerPot<-WDFWThisYearLbsAndCountsConvertedToLbs$Crab.Retained..lbs./WDFWThisYearLbsAndCountsConvertedToLbs$Pots.Fished


#----------------------

Lons<- -1*(WDFWThisYearLbsAndCountsConvertedToLbs$Longitude.Begin.Degrees + WDFWThisYearLbsAndCountsConvertedToLbs$Longitude.Begin.Minutes/60 )
Lats<-WDFWThisYearLbsAndCountsConvertedToLbs$Latitude.Begin.Degrees+ WDFWThisYearLbsAndCountsConvertedToLbs$Latitude.Begin.Minutes/60

#--------------------
# PLOT THE POUNDS DATA THAT ULTIMATELY WE WILL USE FOR WASHINGTON
logbookPlotName<-paste(startYearForFishingSeason,"Washington Logbook Lbs")
PlotCoastForSardine(thisPlotName = logbookPlotName,myXlim=c(-126.0, -123.0), myYlim=c(45.0, 49.0))  # c(-126.0, -123.0),myYlim=c(42.5, 50.5)
title(main= logbookPlotName,cex.main=0.8)
points(Lons,Lats, pch=15,cex=0.2)
dev.copy(png,logbookPlotName)
dev.off()

#-----------------------

WDFWThisYearLbsAndCountsConvertedToLbs$Lats <- Lats
WDFWThisYearLbsAndCountsConvertedToLbs$Lons <- Lons

#----------------------
# REMOVE LOGBOOK ROWS THAT LACK LAT OR LON:  these NAs cause geosphere::distHaversine to fail.
rowsWithLonLat <- which(!is.na(WDFWThisYearLbsAndCountsConvertedToLbs$Lons) & !is.na(WDFWThisYearLbsAndCountsConvertedToLbs$Lats) )

WDFWThisYearLbsAndCountsConvertedToLbs <- WDFWThisYearLbsAndCountsConvertedToLbs[rowsWithLonLat, ]

#----------------

require(geosphere)  # LOAD PACKAGE

#---------------

#AIM NOW IS TO ADD COLUMNS OF BEST MATCHING ROMS JSCOPE OCEAN POINTS AND OCEAN VARIABLES (LAT, LON, AND COLUMNS FOR EACH VARIABLE)
# THIS IS ADDED to the crab logbook data as new columns, sort of like adding ocean observations to what the fisher saw in the trap.

#----------------

for(i in 1:nrow(WDFWThisYearLbsAndCountsConvertedToLbs))
  {
 
  #-------------
   # For this logbook row's fishing month, subset oceanMatrix to only this month.
  calendarMonthNum <- WDFWThisYearLbsAndCountsConvertedToLbs$calendarMonthNum[i]
  rowsForThisMonth <- which(oceanMatrixThisYear[,"monthROMSOceanConditions"]==calendarMonthNum)

  # A quick check to make sure there's ROMS ocean data for this logbook point
   if(length(rowsForThisMonth)==0)
  {
    print(WDFWThisYearLbsAndCountsConvertedToLbs[i,])
    print("row")
    print(i)
    print("month")
     print(calendarMonthNum)
    stop("rows of ocean Jscope data for this fishing month are 0. Check if year (and months) match up in Jscope versus fishery data")
  }
  oceanMatrixThisYearThisMonth <- oceanMatrixThisYear[rowsForThisMonth, ]

  #----------

  ## reformat ocean lat,lon vectors into a matrix
oceanLonLatMatrix <- matrix(c(oceanMatrixThisYearThisMonth[,"nc_lonVector"], oceanMatrixThisYearThisMonth[ ,"nc_latVector"]), ncol = 2)

#  NOW FORMAT LON AND LAT
## subset lat,lon for single crab point
logbookLonLat <- c(WDFWThisYearLbsAndCountsConvertedToLbs$Lons[i],WDFWThisYearLbsAndCountsConvertedToLbs$Lats[i])


# THIS IS THE MAIN FUNCTION FROM MIKE MALICK TO GET DISTANCES BETWEEN LOGBOOK AND ROM POINTS
## calc dist b/w CRAB point and all ocean obs
dist.i <- geosphere::distHaversine(logbookLonLat, oceanLonLatMatrix)

## find index of ocean obs w/ min distance
ocean.ind <- which(dist.i == min(dist.i))[1]

## assign ocean obs w/ min distance to crab data

#-----------
# NOW ADD THE OCEAN VARIABLES AS NEW COLUMNS TO THE CRAB LOBGOOK DATA
#-----------
WDFWThisYearLbsAndCountsConvertedToLbs$distToNearestGridpointM[i] <- dist.i[ocean.ind]
WDFWThisYearLbsAndCountsConvertedToLbs$nc_lat[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "nc_latVector"]
WDFWThisYearLbsAndCountsConvertedToLbs$nc_lon[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "nc_lonVector"]
WDFWThisYearLbsAndCountsConvertedToLbs$bathymetry[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bathymetryVector"]
WDFWThisYearLbsAndCountsConvertedToLbs$bottomTemp[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomTempVector"]
WDFWThisYearLbsAndCountsConvertedToLbs$bottomOxy[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomOxyVector"]
WDFWThisYearLbsAndCountsConvertedToLbs$SSHVector[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "SSHVector"]
WDFWThisYearLbsAndCountsConvertedToLbs$Chla2m[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "Chla2mVector"]
WDFWThisYearLbsAndCountsConvertedToLbs$bottomArag[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomAragVector"]
WDFWThisYearLbsAndCountsConvertedToLbs$bottompH[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottompHVector"]
WDFWThisYearLbsAndCountsConvertedToLbs$bottomSalt[i] <- oceanMatrixThisYearThisMonth[ ocean.ind, "bottomSaltVector"]

#-----------
#   JUST PRINT SOME INFO TO SCREEN PERIODICALLY SO YOU KNOW IT IS RUNNING (SLOWLY)
if(i%%100==0)
{
print("row is")
print(i)
print("dist is")
print(WDFWThisYearLbsAndCountsConvertedToLbs$distToNearestGridpointM[i])
}
#---------------

}


return(WDFWThisYearLbsAndCountsConvertedToLbs)

}
