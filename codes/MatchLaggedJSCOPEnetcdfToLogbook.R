MatchLaggedJSCOPEnetcdfToLogbook <- function(LagoceanMatrixAllYears,WDFW_Crab_Logbooks_Data_ThisYearOnly,startYearForFishingSeason,LagYrs)
{
  #------
  # Primary Author: Isaac Kaplan (Feb 20 2019) 
  # isaac.kaplan at noaa.gov
  # Secondary Author: Emily Norton
  # emilyln@uw.edu
  
  #  Matches the ocean points up to logbook locations, and assigns annually-averaged, lagged ocean variables as new columns of the logbook rows.
  #  Arguments:
  #       LagoceanMatrixAllYears: very large matrix of data for every grid point of ROMS Jscope NC files (annually-averaged, lagged conditions)
  #       WDFW_Crab_Logbooks_Data_ThisYearOnly: single year of Washington logbook data
  #       startYearForFishingSeason:  would be for instance 2017 for fishing year autumn 2107--spring 2018.
  #                                   Important because this is used to find NC file we want to match to this logbook record.
  #       LagYrs: number of years by which to lag the annual J-SCOPE conditions compared to the startYearForFishingSeason
  #
  #   Returns: WDFWThisYearLbsAndCountsConvertedToLbs, each row has logbook record (including cath rate per trap) and associated ocean conditiosn predicted by ROMS nearest point for that year
  
  #----------
  print("STARTING TO MATCH JSCOPE netcdf FILE TO CRAB LOGBOOK FOR 1 YEAR")
  
  source("PlotCoastForSardine.r")
  
  rowsForThisYear <- which(LagoceanMatrixAllYears[,"yearROMSstart"]==startYearForFishingSeason-LagYrs)
  
  if(length(rowsForThisYear)==0)
  {
    stop("Error: rows of ocean Jscope data for this fishing year are 0. Check if year (and months) match up in Jscope versus fishery data")
  }
  

  LagoceanMatrixThisYear <- LagoceanMatrixAllYears[rowsForThisYear, ]
  
  
  require(tidyr)
     #----------------
    require(geosphere)  # LOAD PACKAGE
    #---------------
    #AIM NOW IS TO ADD COLUMNS OF BEST MATCHING ROMS JSCOPE OCEAN POINTS AND OCEAN VARIABLES (LAT, LON, AND COLUMNS FOR EACH VARIABLE)
  # THIS IS ADDED to the crab logbook data as new columns, sort of like adding ocean observations to what the fisher saw in the trap.
    #----------------
  for(i in 1:nrow(WDFW_Crab_Logbooks_Data_ThisYearOnly))
  {
      ## reformat ocean lat,lon vectors into a matrix
    LagoceanLonLatMatrix <- matrix(c(LagoceanMatrixThisYear[,"nc_lonVector"], LagoceanMatrixThisYear[ ,"nc_latVector"]), ncol = 2)
    
    #  NOW FORMAT LON AND LAT
    ## subset lat,lon for single crab point
    logbookLonLat <- c(WDFW_Crab_Logbooks_Data_ThisYearOnly$Lons[i],WDFW_Crab_Logbooks_Data_ThisYearOnly$Lats[i])

    
    # THIS IS THE MAIN FUNCTION FROM MIKE MALICK TO GET DISTANCES BETWEEN LOGBOOK AND ROM POINTS
    ## calc dist b/w CRAB point and all ocean obs
    dist.i <- geosphere::distHaversine(logbookLonLat, LagoceanLonLatMatrix)
    
    
    ## find index of ocean obs w/ min distance
    ocean.ind <- which(dist.i == min(dist.i))[1]
    
    ## assign ocean obs w/ min distance to crab data
    
    #-----------
    # NOW ADD THE OCEAN VARIABLES AS NEW COLUMNS TO THE CRAB LOBGOOK DATA 
    #-----------
    WDFW_Crab_Logbooks_Data_ThisYearOnly$LagbottomTemp[i] <- LagoceanMatrixThisYear[ ocean.ind, "bottomTempVector"]
    WDFW_Crab_Logbooks_Data_ThisYearOnly$LagbottomOxy[i] <- LagoceanMatrixThisYear[ ocean.ind, "bottomOxyVector"]
    WDFW_Crab_Logbooks_Data_ThisYearOnly$LagChla2m[i] <- LagoceanMatrixThisYear[ ocean.ind, "Chla2mVector"]
    
    
    #-----------
    #   JUST PRINT SOME INFO TO SCREEN PERIODICALLY SO YOU KNOW IT IS RUNNING (SLOWLY)
    if(i%%1000==0)
    {
      print("row is")
      print(i)
 #     print("dist is")
#      print(WDFW_Crab_Logbooks_Data_ThisYearOnly$distToNearestGridpointM[i])
    }
    #---------------
    
  }

  
  # Rename these columns based on how many year lag they have...will need to run this script multiple times for different year lags
  names(WDFW_Crab_Logbooks_Data_ThisYearOnly)[names(WDFW_Crab_Logbooks_Data_ThisYearOnly) == 'LagbottomTemp'] <- paste0('bottomTemp_', LagYrs,'yrLag')
  names(WDFW_Crab_Logbooks_Data_ThisYearOnly)[names(WDFW_Crab_Logbooks_Data_ThisYearOnly) == 'LagbottomOxy'] <- paste0('bottomOxy_', LagYrs,'yrLag')
  names(WDFW_Crab_Logbooks_Data_ThisYearOnly)[names(WDFW_Crab_Logbooks_Data_ThisYearOnly) == 'LagChla2m'] <- paste0('Chla2m_', LagYrs,'yrLag')

  
  return(WDFW_Crab_Logbooks_Data_ThisYearOnly)
  
}