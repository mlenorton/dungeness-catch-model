ExtractOceanMatrixFromOneYearNCFile<-function(thisLagNCFileName)
{
  # Primary Author: Isaac Kaplan 
  # isaac.kaplan at noaa.gov
  # Secondary Author: Emily Norton
  # emilyln@uw.edu
  
  # returns ocean matrix for one year, after extracing from NC ROMS file. 

  #------------------
  source("PlotCoastForSardine.r")
  source("PlotNCmap.r")
  
  require(ncdf4)
  require(reshape2)
  require(dplyr)

  yearROMSstart <-as.numeric(str_sub(thisLagNCFileName,-14,-11))

  #----------
  # Open a connection to this File
  #------------
  
  nc <- nc_open(thisLagNCFileName)
  
  #-----------
  
#  # (Optionally) Save the print(nc) dump to a text file (same name as the nc file with a txt extension)
#  {
#    sink(paste0(thisLagNCFileName, ".txt"))
#    #sink(paste0("data/", flist[1], ".txt"))
#    print(nc)
#    sink()
#  }
  
  # Get a list of the NetCDF's R attributes:
  attributes(nc)$names
  
  print(paste("The file has",nc$nvars,"variables,",nc$ndims,"dimensions and",nc$natts,"NetCDF attributes"))
  
  
  # Get a list of the nc variable names.
  attributes(nc$var)$names
  attributes(nc$var)$names[1]
  
  
  #-------------
  #GET VARIABLES
  #------------
  
  bathymetryVector<- as.vector(ncvar_get(nc,"h"))
  nc_latVector<- as.vector(ncvar_get( nc, "lat_rho"))
  nc_lonVector<- as.vector(ncvar_get( nc, "lon_rho"))
  bottomTempVector<- as.vector(ncvar_get(nc,"Bottom_Temp_Avg"))
  bottomOxyVector<- as.vector(ncvar_get(nc,"Bottom_Oxy_Avg"))
  bottomSaltVector<- as.vector(ncvar_get(nc,"Bottom_Salt_Avg"))
  SSHVector<- as.vector(ncvar_get(nc,"SSH_Zeta_Avg"))
  Chla2mVector<- as.vector(ncvar_get(nc,"Chla2m_Avg"))
  bottomAragVector<- as.vector(ncvar_get(nc,"Bottom_Arag_Avg"))
  bottompHVector<- as.vector(ncvar_get(nc,"Bottom_pH_Avg"))
  
  nc_close(nc)  # CLOSE THE NC file cleanly. 
  
  # DATA CHECK ----------------
 
  length(nc_latVector)
  print(paste(length(nc_latVector), "latitudes and", length(nc_lonVector), "longitudes"))
  
  #-------------------
  
  # PLOT SOME MAPS OF VARIABLES HERE (optional)
  
  #----------------
  
  #PlotNCmap(bathymetry, paste("bathymetry",thisMonthsNCFileName),nc_lon,nc_lat)
  #PlotNCmap(bottomOxy, paste("bottomOxy",thisMonthsNCFileName),nc_lon,nc_lat  )
  #PlotNCmap(bottomSalt, paste("bottomSalt",thisMonthsNCFileName),nc_lon,nc_lat  )
  #PlotNCmap(bottomTemp, paste("bottomTemp",thisMonthsNCFileName),nc_lon,nc_lat  )
  #PlotNCmap(SSH, paste("SSH",thisMonthsNCFileName),nc_lon,nc_lat  )
  #PlotNCmap(bottomArag, paste("bottomArag",thisMonthsNCFileName),nc_lon,nc_lat  )
  #PlotNCmap(bottompH, paste("bottompH",thisMonthsNCFileName),nc_lon,nc_lat  )
  #PlotNCmap(Chla2m, paste("Chla2m",thisMonthsNCFileName),nc_lon,nc_lat  )

  
 
  
  # Use this later for matching lat lon points from ROMS to logbooks, one lon/lat location per row
  LagoceanMatrix <- cbind(nc_latVector,nc_lonVector,bathymetryVector,bottomTempVector,bottomOxyVector,bottomSaltVector,SSHVector,Chla2mVector,bottomAragVector,bottompHVector,yearROMSstart)
  
}