MatchEnvironmentLayerToLatLon <- function(LogbookDataFrameWithLatLon, HabitatFile,habitatHeaderName)
{
#----------
# Author: Isaac Kaplan (Feb 20 2019) 
# isaac.kaplan at noaa.gov
  
  # EXAMPLE: OregonLogbooksAndOceanConditionsDEPTH<-MatchEnvironmentLayerToLatLon(OregonLogbooksAndOceanConditions,"ngdc_bath_0to1200_m_pts_geo.csv","DEPTHM") 
  # EXAMPLE: OregonLogbooksAndOceanConditionsHABITAT<-MatchEnvironmentLayerToLatLon(OregonLogbooksAndOceanConditions,"grainsize_pts_geo.csv","GRAINSIZE") 
  
  # Matches habitat to logbook locations

  # INPUTS:
     # LogbookDataFrameWithLatLon is a dataframe from R that can contain logbook columns and also ocean variables. 
     #         Examples include OregonLogbooksAndOceanConditions or WashingtonLogbooksAndOceanConditions
     # HabitatFile is string name of CSV file with habitat data, such as "grainsize_pts_geo.csv"
  
  # RETURNS: LogbookDataFrameWithLatLonHabitat
  
  # WHERE TO CALL THIS: Call this from the main Script (runCrabLogbooks.R) right before you call the gams, and after the logbook points have
     # been combined with the ocean variables. 

#----------
print("Lat Lon matching to habitat")

source("PlotCoastForSardine.r")

HabitatLatLonZ <-read.csv(HabitatFile, header=TRUE)  

#----------------

require(geosphere)  # LOAD PACKAGE

#---------------
## reformat ocean lat,lon vectors into a matrix
habitatLonLatMatrix <- matrix(c(HabitatLatLonZ[,"LONGITUDE"], HabitatLatLonZ[ ,"LATITUDE"]), ncol = 2)

print("number of rows in logbook data: ")
print(nrow(LogbookDataFrameWithLatLon))

LogbookDataFrameWithLatLonHabitat <- LogbookDataFrameWithLatLon

for(i in 1:nrow(LogbookDataFrameWithLatLon)) {
#  NOW FORMAT LON AND LAT
## subset lat,lon for single crab point
logbookLonLat <- c(LogbookDataFrameWithLatLon$Lons[i],LogbookDataFrameWithLatLon$Lats[i])

# THIS IS THE MAIN FUNCTION FROM MIKE MALICK TO GET DISTANCES BETWEEN LOGBOOK AND ROM POINTS
## calc dist b/w CRAB point and all env obs
dist.i <- geosphere::distHaversine(logbookLonLat, habitatLonLatMatrix)

## find index of env obs w/ min distance
hab.ind <- which(dist.i == min(dist.i))[1]

## assign env obs w/ min distance to hake data

#-----------
# NOW ADD THE ENV VARIABLES AS NEW COLUMNS TO THE CRAB LOBGOOK DATA
#-----------
LogbookDataFrameWithLatLonHabitat$distToNearestHabPointM[i] <- dist.i[hab.ind]
LogbookDataFrameWithLatLonHabitat$NewHabLayer[i] <- HabitatLatLonZ[ hab.ind, habitatHeaderName]

#-----------
#   JUST PRINT SOME INFO TO SCREEN PERIODICALLY SO YOU KNOW IT IS RUNNING (SLOWLY)
if(i%%50==0) {
print("row is")
print(i)
print("dist is")
print(LogbookDataFrameWithLatLonHabitat$distToNearestHabPointM[i])
}
#---------------

}

return(LogbookDataFrameWithLatLonHabitat)

}
