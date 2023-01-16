.PBSfigCrabLogbookAno <- function(clm,LogbookDataWithOceanConditions,usePredict,gam.model,yrLabel) {  

  ##  PBSfigCrabLogbookAno.R
  # Primary Author: Isaac Kaplan (May 1 2018)
  # isaac.kaplan@noaa.gov
  # Secondary Author: Emily Norton
  # emilyln@uw.edu
  
  # This function is used to generate figures and anomaly predictions for GAMs based 
  # on logbook data matched to J-SCOPE oceanography. This is a function called by
  # the prediction script "forecastGAM_fishingBehavSelect.R"
  
  library(RColorBrewer)
  
  clr <- .PBSclr(); 
  data(nepacLL,surveyData,envir=sys.frame(sys.nframe()));
  events <- LogbookDataWithOceanConditions
  events$X<-LogbookDataWithOceanConditions$Lons
  events$Y<-LogbookDataWithOceanConditions$Lats
  events$EID<-seq(1:dim(LogbookDataWithOceanConditions)[1])
  clm$X <- clm$Lons   
  clm$Y <- clm$Lats
  clm$EID <- seq(1:dim(clm)[1])
  xl <- c(-125.5,-122);  yl <- c(42.5,48.5)

  # make a grid for crab fishing logbook area
  gwidth <- 0.1   
  grid   <- makeGrid(x=seq(-125.5,-122,gwidth), y=seq(42.5,48.5,gwidth),
                     projection="LL", zone=9)
 
  # locate EventData in grid
  locData<- findCells(events, grid)
  locDataclm <- findCells(clm, grid)     

  #-----------
  # Define Z variable based on raw (usePredict=0) or predicted (usePredict=1) catch rate
  if( usePredict==0)
  {  
    events$Z <- events$KgsPerPot
    clm$Z <- clm$KgsPerPot
  } else if (usePredict==1)
  {
    events$Z <- exp(predict.gam(gam.model,LogbookDataWithOceanConditions,exclude="s(CrabYearS)"))-0.01    # Use average year term for predict mode
    clm$Z <- exp(predict.gam(gam.model,clm,exclude="s(CrabYearS)"))-0.01   
  } else if (usePredict == 2)
  { # Difference betwwen predicted and true, so can see Mean error per point in plot
    events$Z <- exp(predict.gam(gam.model,LogbookDataWithOceanConditions,exclude="s(CrabYearS)"))-0.01-events$KgsPerPot
  } else if (usePredict==3)
  {
    events$Z <- events$KgsPerPot
  }
  #---------------------
  if(usePredict<3)  # Just taking mean catch per pot, per cell
  {
    pdata  <- combineEvents(events, locData, FUN=mean)
    pdataclm <- combineEvents(clm, locDataclm, FUN=mean)
  } else  # count of data per cell
  {
    events$Z[which(events$Z>0)]<-1
    pdata <- combineEvents(events, locData, FUN=sum)
    pdataclm <- combineEvents(clm, locDataclm, FUN=sum)
    print('max of pdata from Counts per polygon')
    print(max(pdata))
  }
  
  #------------------
 # Identify low-count cells to mask (for confidentiality purposes)
  if (usePredict ==0)
  {
    eventsForCounts <- events
    eventsForCounts$Z[which(eventsForCounts$Z>0)]<-1
    
    # Take the sum here, i.e. a count of strings. 
    pdataForCounts <- combineEvents(eventsForCounts, locData, FUN=sum)
    
    # Define minimum number of strings per grid cell in order to NOT mask it. 
    minCountThresh <- 10
    
    # subset the data to only those points we want to color in as mask 
    pdataForCountMask<-pdataForCounts[pdataForCounts$Z<minCountThresh,]
  }  
  
  
  #------------------
  print('usePredict =')
  print(usePredict)
  if(usePredict<2) # plots of kg per pot
  {  
    brks   <- c(0,3,6,9,12,50)   
    brks <- round(brks,2)
    lbrks <- length(brks)
    cpal <- brewer.pal(9,"PuOr")
    cols <- cpal[seq(length(cpal),1,-2)]   # reverse coloring so purple is underpred and orange is overpredict
    legendText <-"Kg per pot"
  } else if(usePredict==2) # plots of mean relative error
  {
    brks   <- c(-10,-2,-1,1,2,10)     
    brks <- round(brks,2)
    lbrks <- length(brks)
    cpal <- brewer.pal(9,"PuOr")
    cols <- cpal[seq(length(cpal),1,-2)] # reverse coloring so purple is underpred and orange is overpredict
    legendText <-"Pred - Obs (kg/pot)"
  } else if (usePredict == 3)
  {
    brks   <- c(0,5,50,100,500,2000)     
    brks <- round(brks,2)
    lbrks <- length(brks)
    cpal <- brewer.pal(9,"Purples")
    cols <- cpal[seq(1,length(cpal),2)]
    legendText <-"Number of strings"
  }
  
  ## To calculate anomaly, generate pdata for clm too (above) and then combine pdata and pdataclm prior to coloring and setting brks:
pdataOG <- pdata
ano <- pdata 
ano$Z <- NaN

uPID <- unique(locData$PID)
uSID <- unique(locData$SID)

for (r in 1:length(pdata[,1])) {   
  p <- pdata[r,"PID"]
  s <- pdata[r,"SID"]
  
  indsp <- which(pdataclm[,1]==p)
  indss <- which(pdataclm$SID==s)
  goodinds <- intersect(indsp,indss)
  
  if (length(goodinds)>0) {
  ano[r,"Z"] <- pdata[r,"Z"] - pdataclm[goodinds,"Z"]
  } else {
    ano[r,"Z"] <- NA
    
  }
  
}

#to make better breaks ---
 anoS <- sort(ano$Z)
 B1 <- anoS[(15/100*length(anoS))]
 B3 <- anoS[(85/100*length(anoS))]
 
 minano <- anoS[1]
 maxano <- anoS[length(anoS)]
 dx <- (B3-B1)/3
 
 brks <- c(minano,-2*dx, -1*dx, dx, 2*dx, maxano)
 if (brks[1]>=brks[2]) {
   brks[1]<- brks[2]-0.1
 }
 if (brks[6]<=brks[5]) {
   brks[6]<- brks[5]+0.1
 }
 brks <- round(brks,1)
 print(brks)
 
# ------
  
 pdata <- ano
 # basic fix for redundant break issue - not foolproof
 if (length(unique(brks))<length(brks)) {
   for (t in 1:length(unique(brks))) {
     inds <-which(brks == brks[t])
     if (length(inds)>1) {
       brks[inds[2]]<- brks[inds[2]]+0.1    ## this assumes that only two of the break points are identical, and it nudges the second one slightly higher
   }
 }
 }  
 
pdata  <- makeProps(pdata, brks, "col", cols)     

  par(mfrow=c(1,1),omi=c(0,0,0,0)) 
  
  #------Plot-the-figure------
  plotMap(nepacLL, col=clr$land, bg=clr$sea, xlim=xl, ylim=yl, tck=-0.015,
          mgp=c(2,.5,0), cex=1.2, plt=c(.08,.98,.08,.98))
  addPolys(grid, polyProps=pdata)
  
  #----------------
  # Masks the cells with less than minimum number of strings. 
  # Only do this if we are plotting raw data (usePredict 0)
  if( usePredict==0)
  {
    colsForMask <- c("#969999")    # medium gray in hexidecimal format
    addPolys( grid, polyProps =pdataForCountMask, "col"=colsForMask,border="gray")
  }  # end of masking (recoloring) of these cells

  
  if(usePredict==0) 
  {text(xl[1]+.9,yl[2]-.45,paste("Obs - Clm"),cex=0.9,adj=0)
  } else if(usePredict==1)
  {text(xl[1]+.9,yl[2]-.45,paste("Fore - Clm"),cex=0.9,adj=0)
  } 
  legend("bottomright", inset = 0.02, title=legendText, legend=paste(brks[1:(lbrks-1)],brks[2:lbrks], sep=" to "), fill=cols, box.lty=0)

## Uncomment if you want to save .pngs automatically
#  dev.copy(png,paste(yrLabel,'Predict_',usePredict,'_Logbook_Ano.png',sep=""))
#  dev.off()
#  dev.copy(png, paste('CatchRates, ','usePredict is ' ,usePredict,".png"))
  
  
## Return anomalies to calculate anomaly correlation coefficient:
  ano     # export the entire ano dataframe so we have PID and SID too

} # end function

 
