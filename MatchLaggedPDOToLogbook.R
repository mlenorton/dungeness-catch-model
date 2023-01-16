MatchLaggedPDOToLogbook <- function(PDOAllYears,WDFW_Crab_Logbooks_Data_ThisYearOnly,startYearForFishingSeason,LagYrsPDO)
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
  #       WDFW_Crab_Logbooks_Data_ThisYearOnly: single year of washington logbook data
  #       startYearForFishingSeason:  would be for instance 2017 for fishing year autum 2107--spring 2018.
  #                                   Important because this is used to find NC file we want to match to this logbook record.
  #       LagYrsPDO: number of years by which to lag the annual J-SCOPE conditions compared to the startYearForFishingSeason
  #
  #   Returns: WDFWThisYearLbsAndCountsConvertedToLbs,  each row has logbook record (including cath rate per trap) and associated PDO for that year
  
  #----------
  
  source("PlotCoastForSardine.r")
  
  rowsForThisYear <- which(PDOAllYears[,"YEAR"]==startYearForFishingSeason-LagYrsPDO)
  
  if(length(rowsForThisYear)==0)
  {
    stop("Error: rows of PDO data for this fishing year are 0. Check if year matchs up in PDO file versus fishery data")
  }
  
  PDOThisYear <- PDOAllYears[rowsForThisYear, ]
  
  
  require(tidyr)
  #-----------
  # NOW ADD THE PDO VALUE AS NEW COLUMN TO THE CRAB LOBGOOK DATA 
  #-----------
  for(i in 1:nrow(WDFW_Crab_Logbooks_Data_ThisYearOnly))
  {
    #-----------
    # NOW ADD THE PDO VALUE AS NEW COLUMN TO THE CRAB LOBGOOK DATA 
    #-----------
    WDFW_Crab_Logbooks_Data_ThisYearOnly$PDO[i] <- PDOThisYear["PDO"]

    #-----------
    #   JUST PRINT SOME INFO TO SCREEN PERIODICALLY SO YOU KNOW IT IS RUNNING (SLOWLY)
    if(i%%100==0)
    {
      print("row is")
      print(i)
#      print("dist is")
#      print(WDFW_Crab_Logbooks_Data_ThisYearOnly$distToNearestGridpointM[i])
    }
    #---------------
    
  }

  # Rename these columns based on how many year lag they have
  names(WDFW_Crab_Logbooks_Data_ThisYearOnly)[names(WDFW_Crab_Logbooks_Data_ThisYearOnly) == 'PDO'] <- paste0('PDO_', LagYrsPDO,'yrLag')
 
  
  return(WDFW_Crab_Logbooks_Data_ThisYearOnly)
  
}