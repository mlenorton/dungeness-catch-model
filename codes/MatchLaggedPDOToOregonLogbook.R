MatchLaggedPDOToOregonLogbook <- function(PDOAllYears,Oregon_Crab_Logbooks_Data_ThisYearOnly,startYearForFishingSeason,LagYrsPDO)
{
  #----------
  # Primary Author: Isaac Kaplan (Feb 20 2019) 
  # isaac.kaplan at noaa.gov
  # Secondary Author: Emily Norton
  # emilyln@uw.edu
  #------
  #  Matches PDO to logbook entries, and assigns PDO as new column of the logbook rows.
  #  Arguments:
  #       PDOAllYears: matrix with summed (Jan-July) PDO values for various years
  #       Oregon_Crab_Logbooks_Data_ThisYearOnly: single year of washington logbook data
  #       startYearForFishingSeason:  would be for instance 2017 for fishing year autum 2107--spring 2018.
  #                                   Important because this is used to find NC file we want to match to this logbook record.
  #       LagYrsPDO: number of years by which to lag the annual J-SCOPE conditions compared to the startYearForFishingSeason
  #
  #   Returns: OregonThisYearLbsAndOceanConditions,  each row has logbook record (including cath rate per trap) and associated PDO for that year
  
  #----------
  #-------------
  # NOTE THIS IS A SLIGHLY MODIFIED VERSION OF MatchLaggedPDOToLogbook(  ), which was for Washington.
  # SEE THAT FILE FOR BETTER COMMENTS.
  # THIS TWEAKED VERSION IS NECESSARY BECAUSE OREGON LOGBOOK DATA LOOKS DIFFERENT FROM WASHINGTON LOGBOOKS.
  #-------------

  #--------------
source("PlotCoastForSardine.r")

rowsForThisYear <- which(PDOAllYears[,"YEAR"]==startYearForFishingSeason-LagYrsPDO)
if(length(rowsForThisYear)==0) {
  stop("rows of ocean Jscope data for this fishing year are 0. Check if year (and months) match up in Jscope versus fishery data")
}

PDOThisYear <- PDOAllYears[rowsForThisYear, ]

require(tidyr)

OregonThisYearLbsAndOceanConditions <- Oregon_Crab_Logbooks_Data_ThisYearOnly
OregonThisYearLbsAndOceanConditions <- OregonThisYearLbsAndOceanConditions[which(!is.na(OregonThisYearLbsAndOceanConditions$AdjPounds)) ,]

#----------------------
# REMOVE ROWS THAT LACK LAT OR LON: 
rowsWithLonLat <- which(!is.na(OregonThisYearLbsAndOceanConditions$Lons) & !is.na(OregonThisYearLbsAndOceanConditions$Lats) )
OregonThisYearLbsAndOceanConditions <- OregonThisYearLbsAndOceanConditions[rowsWithLonLat, ]

#----------------

require(geosphere)

for(i in 1:nrow(OregonThisYearLbsAndOceanConditions)) {
 ## assign logbook obs to PDO based on year only (neither month nor spatial constraints need to be considered)

OregonThisYearLbsAndOceanConditions$LagPDO[i] <- PDOThisYear["PDO"]

if(i%%100==0) {
print("row is")
print(i)
}


}

# Rename these columns based on how many year lag they have...will need to run this script multiple times for different year lags
names(OregonThisYearLbsAndOceanConditions)[names(OregonThisYearLbsAndOceanConditions) == 'LagPDO'] <- paste0('PDO_', LagYrsPDO,'yrLag')


return(OregonThisYearLbsAndOceanConditions)

}
