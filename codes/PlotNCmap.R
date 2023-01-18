PlotNCmap <- function (ncVariable, myTitle, nc_lon, nc_lat)
{
  
  ## PlotNCmap.R
  # Primary Author: Isaac Kaplan (May 1 2018)
  # isaac.kaplan@noaa.gov
  # Secondary Author: Emily Norton
  # emilyln@uw.edu
  
  # This is a plotting function (optionally) called by 
  # ExtractOceanMatrixFromOneYearNCFile.R
  
   # JUST PLOTS A GENERIC MAP OF US WEST COAST, USING OLD SARDINE MAPPING CODE 
  PlotCoastForSardine(myTitle)
  title(main= myTitle,cex.main=1.5)
  rbPal <-colorRampPalette(c('blue','red'))    
  
  # LOOK AT RANGE OF NCVARIABLE PASSED IN HERE
  thisMax<-max(ncVariable,na.rm = TRUE)
  thisMin<-min(ncVariable,na.rm = TRUE)
  thisMaxRounded <- round(thisMax,digits=1)
  thisMinRounded <- round(thisMin,digits=1)
  
   # COLOR IN THE MAP BASED ON VALUES FROM THE NCVARIABLE PASSED IN HERE. 
 thisRange<- (thisMax - thisMin)
  colorsForPredictions <- rbPal(10)[as.numeric(cut(ncVariable,breaks = seq(min(ncVariable,na.rm = TRUE),max(ncVariable,na.rm = TRUE),thisRange/10)))]
  points(nc_lon, nc_lat, pch=15,cex=0.2,col=colorsForPredictions)
  legend(x=-127.5,y=44.5,legend= c(thisMaxRounded,thisMinRounded),fill=c(rbPal(10)[10],rbPal(10)[1]))  
  dev.off()
  
}