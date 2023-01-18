PlotCoastForSardine <- function(thisPlotName,myXlim = c(-128.0, -123.0),myYlim=c(42.5, 50.5))

{
  ## PlotCoastForSardine.R
  # Primary Author: Isaac Kaplan (Dec 29 2008)
  # isaac.kaplan@noaa.gov
  # Secondary Author: Emily Norton
  # emilyln@uw.edu
  
  # This script uses PBSmapping to plot polygons along the US West Coast, and 
  # color them in. See the PBSmapping documentation for extensive info on the
  # package and available data. This function is called by scripts that match
  # logbooks to dynamic and lagged ocean conditions (e.g., 
  # ExtractOceanMatrixFromOneMonthNCFile.R, MatchLaggedPDOToLogbook.R,
  # MatchLaggedPDOToOregonLogbook.R, MatchJSCOPEnetcdfToLogbook.R, 
  # MatchLaggedJSCOPEnetcdfToLogbook.R, MatchJSCOPEnetcdfToOregonLogbook.R
  # MatchLaggedJSCOPEnetcdfToOregonLogbook.R)


require(PBSmapping)
palette("default")
data(nepacLL); # load the nepacLL data set for our region
png(filename=paste(thisPlotName,".png",sep=""))
plotMap(nepacLL, # plot the nepacLL data set
xlim=myXlim, # limit the region horizontally
ylim=myYlim, # limit the region vertically 
col=rgb(255, 255, 195, # set foreground colour
maxColorValue=255),
bg=rgb(224, 253, 254, # set background colour
maxColorValue=255),
tck=c(-0.03), # set tick mark length
cex = 1.5, # set font size
mgp=c(1.9, 0.4, 0),  # 1.9 0.7 0
border="black"); # adjust the axis label locations

}
