.PBSfigCrabLogbook <- function(LogbookDataWithOceanConditions,usePredict,gam.model,yrLabel) {  

  ##  PBSfigCrabLogbook.R
  # Primary Author: Isaac Kaplan (May 1 2018)
  # isaac.kaplan@noaa.gov
  # Secondary Author: Emily Norton
  # emilyln@uw.edu
  
  # This function is used to generate figures and raw predictions for GAMs based 
  # on logbook data matched to J-SCOPE oceanography. This is a function called by
  # the prediction script "forecastGAM_fishingBehavSelect.R"
    
  library(RColorBrewer)
  clr <- .PBSclr(); 
  
  data(nepacLL,surveyData,envir=sys.frame(sys.nframe()));
  events <- LogbookDataWithOceanConditions
  events$X<-LogbookDataWithOceanConditions$Lons
  events$Y<-LogbookDataWithOceanConditions$Lats
  events$EID<-seq(1:dim(LogbookDataWithOceanConditions)[1])
  xl <- c(-125.5,-122);  yl <- c(42.5,48.5)

  # make a grid for crab fishing logbook area
  grid   <- makeGrid(x=seq(-125.5,-122,.1), y=seq(42.5,48.5,.1),
                     projection="LL", zone=9)
  # locate EventData in grid
  locData<- findCells(events, grid)
  
  #-----------
  # Define Z variable based on raw or predicted catch rate
  if( usePredict==0)
  {  
    events$Z <- events$KgsPerPot
  } else if (usePredict==1)
  {
    events$Z <- exp(predict.gam(gam.model,LogbookDataWithOceanConditions,exclude="s(CrabYearS)"))-0.01
    #events$Z <- exp(gam.model$fitted.values)-0.01
  } else if (usePredict==2)
  {
    # Difference betwwen predicted and true, so can see Mean error per point in plot
    events$Z <- exp(predict.gam(gam.model,LogbookDataWithOceanConditions,exclude="s(CrabYearS)"))-0.01-events$KgsPerPot   # Use average year term in forecast mode
  } else if (usePredict==3)
  {
    events$Z <- events$KgsPerPot    # note:  sum up the number of observations (KgsPerPot>0) below to get an estimate of the avg number of strings 
  } else if (usePredict==4)     
  {
    events$Z <- events$KgsPerPot    # note:  sum up the number of observations (KgsPerPot>0) below to get an estimate of the std of strings 
  } else if (usePredict==5)     
  {
    events$Z <- events$KgsPerPot    # note: we sum up the number of observations (KgsPerPot>0) below to get an estimate of the number of strings 
  } else if (usePredict==6)     # to plot percent difference between oberved and forecast
  {
    # Calculate the % difference after the average differences have been calculated
    events2 <- events
    events$Z <- exp(predict.gam(gam.model,LogbookDataWithOceanConditions,exclude="s(CrabYearS)"))-0.01-events$KgsPerPot
    events2$Z <- events$KgsPerPot
  } 
  
  #---------------------
  if(usePredict==3)  # count of data per cell 
  {
    events$Z[which(events$Z>0)]<-1
    pdata <- combineEvents(events, locData, FUN=sum)
    print('max of pdata from Counts per polygon')
    print(max(pdata))
  } else if (usePredict==4)  # avg annual number of strings per cell 
  {
    events$Z[which(events$Z>0)]<-1
    pdata <- combineEvents(events, locData, FUN=sum)
    pdata$Z <- pdata$Z/(length(unique(LogbookDataWithOceanConditions$CrabYearS)))
    print('max of pdata from Counts per polygon')
    print(max(pdata))  # Take mean catch per pot, per cell
  } else if (usePredict==6)      # calculate the % diff between obs + fore, after averaging over the grid cell:
  { 
    avgdiff  <- combineEvents(events, locData, FUN=mean)
    avgobs <- combineEvents(events2, locData, FUN=mean)
    pdata <- avgdiff
    pdata$Z <- 100*(avgdiff$Z)/avgobs$Z
    print(head(pdata))
  } else {
    pdata  <- combineEvents(events, locData, FUN=mean)
  }
  
  #-------------------
  # Identify low-count cells to mask (for confidentiality purposes)  
  if (usePredict ==0)
  {
    eventsForCounts <- events
    eventsForCounts$Z[which(eventsForCounts$Z>0)]<-1
    
    # Take the sum here, essentially as a count of strings. 
    pdataForCounts <- combineEvents(eventsForCounts, locData, FUN=sum)
    
    # Define minimum number of strings per grid cell in order to NOT mask it. 
    minCountThresh <- 10
    
    # subset the data to only those points we want to color in as mask 
    pdataForCountMask<-pdataForCounts[pdataForCounts$Z<minCountThresh,]
    
  } else if (usePredict==6)      
  { eventsForCounts <- events2
  eventsForCounts$Z[which(eventsForCounts$Z>0)]<-1
  
  # Take the sum here, essentially as a count of strings. 
  pdataForCounts <- combineEvents(eventsForCounts, locData, FUN=sum)
  
  # Define minimum number of strings per grid cell in order to NOT mask it. 
  minCountThresh <- 10
  
  # subset the data to only those points we want to color in as mask 
  pdataForCountMask<-pdataForCounts[pdataForCounts$Z<minCountThresh,]
    }# end of this identification of cells to mask

  #------------------
  print('usePredict =')
  print(usePredict)
  if(usePredict<2) # plots of kg per pot
  {  
    brks   <- c(0,3,6,9,12)   
    cpal <- brewer.pal(7,'Greens')
    lbrks <- length(brks)
    cols <- cpal[seq(1,length(cpal),2)]
    legendText <-"Kg per pot"
  } else if(usePredict==2) # plots of mean relative error
  {
    brks   <- c(-10,-2,-1,1,2,10)
    lbrks <- length(brks)
    cpal <- brewer.pal(9,"PuOr")
    cols <- cpal[seq(length(cpal),1,-2)]   # reverse coloring so purple is underpred and orange is overpredict
    legendText <-"Pred - Obs (kg/pot)"
  } else if (usePredict==3)
  {
    brks   <- c(0,5,50,100,500,2000)
    lbrks <- length(brks)
    cpal <- brewer.pal(9,"Purples")
    cols <- cpal[seq(1,length(cpal),2)]
    legendText <-"Number of strings"
  } else if (usePredict==4)
  {
    brks   <- c(0,5,50,100,500,2000)
    lbrks <- length(brks)
    cpal <- brewer.pal(9,"Purples")
    cols <- cpal[seq(1,length(cpal),2)]
    legendText <-"Avg # of strings"
  } else if (usePredict==5)    # for std of strings
  {
    brks   <- c(0,5,50,100,500)
    lbrks <- length(brks)
    cpal <- brewer.pal(9,"Purples")
    cpal <- brewer.pal(9,"Oranges")
    cols <- cpal[seq(1,length(cpal),2)]
    legendText <-"Stdev of strings"
  } else if (usePredict==6)    # for percent difference in forecast versus observed
  {
    print(min(pdata))
    brks <- c(-80, -50, -30, -10, 10, 30, 50, 70) 
    lbrks <- length(brks)
    cpal <- brewer.pal(7,"PuOr")
    cols <- cpal[seq(length(cpal),1,-1)] 
    legendText <-"% Bias"
  }
  
  
  pdata  <- makeProps(pdata, brks, "col", cols)
  par(mfrow=c(1,1),omi=c(0,0,0,0)) 
  par(mar=c(6,4,7,2))
  
  #------Plot-the-figure------
  plotMap(nepacLL, col=clr$land, bg=clr$sea, xlim=xl, ylim=yl, tck=-0.015,
          mgp=c(2,.5,0), cex=1.2, plt=c(.08,.98,.08,.98))
  addPolys(grid, polyProps=pdata)
  
  #----------------
  #Added July 19 2021
  # Masks the cells with less than minimum number of strings. 
  # Only do this if we are plotting raw data (usePredict 0) * OR * bias (b/c raw data could be estimated when forecast + bias are presented together)
  if( usePredict==0 | usePredict==6)   
  {
    colsForMask <- c("#969999")    # medium gray in hexidecimal format
    addPolys( grid, polyProps =pdataForCountMask, "col"=colsForMask,border="gray")
  }  
  

  if(usePredict==0) 
    {text(xl[1]+.9,yl[2]-.45,paste("Observed"),cex=0.9,adj=0)
    } else if(usePredict==1)
    {text(xl[1]+.9,yl[2]-.45,paste("Forecast"),cex=0.9,adj=0)
    } else if(usePredict==2)
    {text(xl[1]+.9,yl[2]-.45,paste("Difference"),cex=0.9,adj=0)
    }  else if (usePredict==3)
    {text(xl[1]+.9,yl[2]-.45,paste("Observed"),cex=0.9,adj=0)
    }   else if (usePredict==4)
    {text(xl[1]+.9,yl[2]-.45,paste("Observed"),cex=0.9,adj=0)
    }   else if (usePredict==5)
    {text(xl[1]+.9,yl[2]-.45,paste("Observed"),cex=0.9,adj=0)}
  
  # add a legend; right-justify the legend labels
  legend("bottomright", inset = 0.02, title=legendText, legend=paste(brks[1:(lbrks-1)],brks[2:lbrks], sep=" to "), fill=cols, box.lty=0)

## If you want to save plots automatically:  
#  dev.copy(png,paste(yrLabel,'Predict ',usePredict,'_Logbook.png',sep=""))
#  dev.copy(png, paste('CatchRates, ','usePredict is ' ,usePredict,".png"))
#  dev.off()
  

  pdata    # export pdata for analysis 
} # end function