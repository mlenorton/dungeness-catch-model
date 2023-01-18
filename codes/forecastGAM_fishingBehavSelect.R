## forecastGAM_fishingBehavSelect.R
# Primary Author: Emily Norton (Dec 15 2020)
# emilyln@uw.edu
# Secondary Author: Isaac Kaplan
# isaac.kaplan@noaa.gov

# This script predicts crab CPUE on the PBS mapping grid when different 
# fishing behaviors are selected (i.e. "fishingBehavSelect"). The user may also 
# choose to subsample (perhaps repeatedly) those fishing behaviors. This script 
# should be run after "match_allFishingData_forecastYears.R" is complete.

# ------------------------------------------------------------------------------

# The following packages are required for this script:
require(mgcv)
require(dplyr)
require(PBSmapping)
require(basicTrendline)
source("PBSfigCrabLogbookAno.R")
source("PBSfigCrabLogbook.R")


figdir <- "../Output/AllFishingObs_SubsamplingApproach/"


fishingBehavSelect <- 1      # select the type of fishing behavior to try: 1 = use observed fishing behavior from the "forecast" year, 
                             # 2 = use fishing behavior from one year prior, 
                             # 3 = use fishing behavior from all years prior (option to subset and subsample multiple times)

resampn <- 1   # the number of times you'd like to re-sample the dataframe to generate predictions (>=1)
indsn <- 30000   # number of rows to subsample each time (avg number of obs for 2009-2015 (when both WA + OR fishing was reported) was ~31400); if you don't want to subsample, use 1
years <-  c(2016:2018)     # the years you'd like to forecast for

obsfig <- 'T'       # 'T' to generate figs for the observed CPUE anomaly; 'F' to skip over this section of code
predfig <- 'T'      # 'T' to generate figs for the predicted (i.e. forecast) CPUE anomaly; 'F' to skip over this section of code


# Load data frame with PBS mapping coordinates for all grid cells which have been fished historically (2007-2018)
load("../Output/PBSMappingPIDSIDdf.RData")   # makes it easier to calculate skill for forecast 

# Load Observed fishing behaviors ()
load(file="../Output/WA_OR_LogbookDataWithOceanConditionsHistoricalFuture_reallyreadyforGAM.RData")   # load dfs with all of the observed fishing behaviors + corresponding ocn conditions from 2007-2018: LogbookDataWithOceanConditionsHistorical,LogbookDataWithOceanConditionsFuture



for (f in 1:length(years)) {   

  forecastYear <- years[f]  
  print(forecastYear)
  
# Pre-allocate dfs for each new forecast year
ano_fore <- PBSmapPIDSIDdf
raw_fore <- PBSmapPIDSIDdf
ano_fore_WA <- ano_fore
ano_fore_OR <- ano_fore
avg_ano_fore <- array(NA,c(length(years),resampn))
avg_ano_fore_WA <- avg_ano_fore
avg_ano_fore_OR <- avg_ano_fore

ano_obs <- array(NA,c(320,3))     #ano_fore
raw_obs <- array(NA,c(320,3))
ano_obs_WA <- ano_obs
ano_obs_OR <- ano_obs
avg_ano_obs <- NA
avg_raw_obs <- NA


#Load matrix with all historical fishing obs matched to J-SCOPE ocn for the forecast year:
load(file=paste("../Output/AllFishingObs_SubsamplingApproach/forecastLaggedPDOLogbookDataWithOceanConditions_forecastYear", forecastYear, ".Rdata",sep="")) 

subLogbookData <- na.omit(forecastLaggedPDOLogbookDataWithOceanConditions)


# Subsample the data and make predictions
for (r in c(1:resampn)) {  

if (forecastYear == 2018) {	# since we only have 2018 data for WA, let's just subsample WA fishing behaviors
subLogbookData <- subLogbookData[which(subLogbookData$state=='WA'),]     
} else if (forecastYear == 2007| forecastYear ==2008) {  # since we only have 2007 and 2008 data for OR, let's just subsample OR fishing behaviors
  subLogbookData <- subLogbookData[which(subLogbookData$state=='OR'),]     
    }
 
	# Determine the kind of fishing behavior we're going to use for forecasting, based on the "fishingBehavSelect" we chose above:
	if (fishingBehavSelect == 1) {
	  yearsforfore <- forecastYear
	} else if (fishingBehavSelect == 2) {
	  yearsforfore <- forecastYear-1
	} else if (fishingBehavSelect == 3) {
	  yearsforfore <- c(2009:(forecastYear-1))    # start with 2009 since we only have fishing behavior for WA + OR starting this year
	}
	
 goodmats <- NA 
    for (y in c(1:length(yearsforfore))) {
        temp <- subLogbookData[which(subLogbookData$CrabYearS==yearsforfore[y]),] 
	      goodmats <- rbind(goodmats,temp)
    }
  subLogbookData <- goodmats[2:dim(goodmats)[1],]    # get rid of first NA row   

 rm(goodmats)

 if (indsn>1){	
 subinds <- sample(c(1:dim(subLogbookData)[1]),indsn)    # Sample without replacement for now
 LogbookDataWithOceanConditionsThisYR <- subLogbookData[subinds,]
 } else {
   LogbookDataWithOceanConditionsThisYR <- subLogbookData
 }
 
LogbookDataWithOceanConditionsThisYR <- LogbookDataWithOceanConditionsThisYR %>%
  mutate(CrabYearS = as.numeric(forecastYear))  # for the forecast, the CrabYearS should be the same as the forecastYear (rather than the year the fishing obs were grabbed from) - just temporarily

#-----------------------------
# NOW MAKE GAM-PREDICTED PLOTS, PER YEAR

if (obsfig == 'T') {      # if we want to make the obs figure.... load the fishing obs from all years, etc.


# Just on the first iteration, make the figure showing the anomaly of the observations compared to "climatological" catch for 2007-2015
if (r==1 & forecastYear<2019){    
 if (forecastYear>2006 & forecastYear<2016){ 
	ObsThisYR <- LogbookDataWithOceanConditionsHistorical[which(LogbookDataWithOceanConditionsHistorical$CrabYearS == forecastYear),]
 } else if (forecastYear>2015 & forecastYear<2019){
        ObsThisYR <- LogbookDataWithOceanConditionsFuture[which(LogbookDataWithOceanConditionsFuture$CrabYearS == forecastYear),]
 }

  # MAKE MAP of Anomaly for observations (for 2007-2018)
   dev.new()
   #tiff(paste(figdir, 'Map_CPUE_Obs_Ano_',forecastYear,'.tiff',sep=""), width=3.5, height=7, units='in', res=300)
   ano_obs <- .PBSfigCrabLogbookAno(LogbookDataWithOceanConditionsHistorical,ObsThisYR,usePredict = 0,yrLabel = as.character(forecastYear)) # note: we need to keep PID and SID since they won't be the same between obs and fore (predicted with subsampled fishing obs)  
   avg_ano_obs <- mean(ano_obs$Z,na.rm=T)   # calculate average CPUE anomaly for the obs over the entire domain
   # dev.off()
  
## Also grab raw values for the obs (so we can calculate correlation coefficient below, in addition to ACC):
    raw_obs <- .PBSfigCrabLogbook(ObsThisYR,usePredict = 0,yrLabel = as.character(forecastYear))
    avg_raw_obs <- mean(raw_obs$Z,na.rm=T)   # calculate average CPUE raw values for the obs over the entire domain
    # dev.off()

 
}
#save(ano_obs,avg_ano_obs,file=paste("../Output/AllFishingObs_SubsamplingApproach/PBSmapping_SubsampledPredictions_ano_obs_", forecastYear, ".Rdata",sep=""))
}    # end of if (obsfig == T)

   
   # MAKE MAP FIGURE. NOTE NOW USING USEPREDICT =1. 
if (predfig == 'T') {
       #dev.new()
 
       # Load the GAM we're going to use for forecasting
       load(file="../Output/gamDiS_SDL_less.RData")   # gamDiS_SDL_less
       gammod <- gamDiS_SDL_less     

       dev.new()
       #	tiff(paste(figdir, 'Map_CPUE_Fore_Ano_',forecastYear,'_ObsFishingBehav_bestGAM.tiff',sep=""), width=3.5, height=7, units='in', res=300)         
       ano_yr <-  .PBSfigCrabLogbookAno(LogbookDataWithOceanConditionsHistorical,LogbookDataWithOceanConditionsThisYR,usePredict = 1,gam.model=gammod,yrLabel = as.character(forecastYear))
       #dev.off() 
       dev.new()
       raw_yr <-  .PBSfigCrabLogbook(LogbookDataWithOceanConditionsThisYR,usePredict = 1,gam.model=gammod,yrLabel = as.character(forecastYear))
       #dev.off() 

ano_fore[,r+2] <- NA   

for (l in 1:dim(ano_yr)[1]) {
  pid <- ano_yr$PID[l]
  sid <- ano_yr$SID[l]
  pinds <- which(ano_fore$PID==pid)
  sinds <- which(ano_fore$SID==sid)
  rowind <- intersect(pinds,sinds)
    
    if (length(rowind)==1) {
      ano_fore[rowind,r+2] <- ano_yr$Z[l]
      raw_fore[rowind,r+2] <- raw_yr$Z[l]     # since we're using the same fishing obs here, we should be able to do the same for the raw forecast as the ano fore
    }
}     #Below calculate the average/std of these Z values (columns 3:resampn+2) 

	     avg_ano_fore[f,r] <- mean(ano_yr$Z,na.rm=T) 
	         
	  #save(ano_fore,avg_ano_fore,gammod, file=paste("../Output/Ano_fore_",forecastYear,"_PBSmap_output.Rdata"))
rm(LogbookDataWithOceanConditionsThisYR)     
 }  # if (predfig == t)
}   # for loop through resampn (r)


## Calculate ACC - to compare observed and forecast CPUE - but first, will need to average CPUE for each grid cell (and calc std) for predicted fishing since multiple subsamplings were conducted 

if (obsfig=='F') {
       load(file=paste("../Output/AllFishingObs_SubsamplingApproach/PBSmapping_SubsampledPredictions_ano_obs_",forecastYear,".Rdata",sep=""))
}


require(dplyr)

# calculate the average across multiple subsampling attempts (if necessary) and rename col - do this for the ano and raw dfs in parallel
maxcol <- paste('V',resampn+2,sep='')     

ano_foreSum <- ano_fore %>% 
	      rowwise() %>% 
	      mutate(m = mean(c_across(V3:maxcol),na.rm=T)) %>%    ## must modify based on sumbsampn, e.g., if reps = 30, use V3:V32; if reps = 50, use V3:V52, etc.
        mutate(sd =sd(c_across(V3:maxcol),na.rm=T)) 	

ano_foreSum <- ano_foreSum %>%
  dplyr::select(PID, SID, m)

raw_foreSum <- raw_fore %>%
        rowwise() %>%
        mutate(m = mean(c_across(V3:maxcol),na.rm=T)) %>%    # if reps = 30, use V3:V32; if reps = 50, use V3:V52, etc.
        mutate(sd =sd(c_across(V3:maxcol),na.rm=T))     # going to need to modify based on sumbsampn

raw_foreSum <- raw_foreSum %>%
        dplyr::select(PID, SID, m)

# Calculate anomaly correlation coefficient (ACC) for the anomaly values and correlation coefficient (r) for the raw values
# first prep dfs so that the obs and fores line up (i.e. make sure PID and SID are the same)
ano_obs_thisyr <- ano_obs
ano_fore_thisyr <- ano_foreSum
ano_obs_fore_thisyr <- ano_obs_thisyr
ano_obs_fore_thisyr <- ano_obs_fore_thisyr %>%
    rename(Zobs=Z)
ano_obs_fore_thisyr$Zfore <- NA

for (n in 1:dim(ano_obs_thisyr)[1]) {
  pid <- ano_obs_thisyr$PID[n]
  sid <- ano_obs_thisyr$SID[n]
    pinds <- which(ano_fore_thisyr$PID==pid)
    sinds <- which(ano_fore_thisyr$SID==sid)
    rowind <- intersect(pinds,sinds)
    
    if (length(rowind)==1) {
      ano_obs_fore_thisyr$Zfore[n] <- ano_fore_thisyr$m[rowind]
    }
}

goodrows <- which(!is.na(ano_obs_fore_thisyr$Zfore))   # grab only the rows for which there are both obs and fore anomalies
accval <- cor.test(ano_obs_fore_thisyr$Zfore[goodrows],ano_obs_fore_thisyr$Zobs[goodrows])
acc_thisyr <- accval$estimate

print('ACC is')
print(acc_thisyr)

# do the same thing for the raw vals to calculate correlation coeff.
raw_obs_thisyr <- raw_obs
raw_fore_thisyr <- raw_foreSum
raw_obs_fore_thisyr <- raw_obs_thisyr
raw_obs_fore_thisyr <- raw_obs_fore_thisyr %>%
    rename(Zobs=Z)
raw_obs_fore_thisyr$Zfore <- NA


for (n in 1:dim(raw_obs_thisyr)[1]) {
  pid <- raw_obs_thisyr$PID[n]
  sid <- raw_obs_thisyr$SID[n]
    pinds <- which(raw_fore_thisyr$PID==pid)
    sinds <- which(raw_fore_thisyr$SID==sid)
    rowind <- intersect(pinds,sinds)

    if (length(rowind)==1) {
      raw_obs_fore_thisyr$Zfore[n] <- raw_fore_thisyr$m[rowind]
    }
}

goodrows <- which(!is.na(raw_obs_fore_thisyr$Zfore))   # grab only the rows for which there are both obs and fore anomalies
rccval <- cor.test(raw_obs_fore_thisyr$Zfore[goodrows],raw_obs_fore_thisyr$Zobs[goodrows])
rcc_thisyr <- rccval$estimate

print('R corr coeff is')
print(rcc_thisyr)

}    # for loop through f years
