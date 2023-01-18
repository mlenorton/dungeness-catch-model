## TestGAMPerformance_PredVsExplSkill_Zwolinski.R
# Primary Author: Isaac Kaplan 
# isaac.kaplan@noaa.gov
# Secondary Author: Emily Norton (Nov 6 2019)
# emilyln@uw.edu

# This code was originally adapted from Zwolinski et al. 2011 to test 
# predictive and explanatory skill of the GAM. We subset our data two ways 
# (i.e., by year and random proportion) to develop cross-validation groups.

# ------------------------------------------------------------------------------

# This script requires the following packages:
library(ncdf4)
library('tidyr')
library('geosphere')
library('tidyverse')
library('maptools')

require(mgcv)
setwd('~/Documents/DungenessCPUEModel/Rcode/')

load(file="../Output/WA_OR_LogbookDataWithOceanConditionsHistoricalFuture_readyforGAM.RData")

LogbookDataWithOceanConditions <- LogbookDataWithOceanConditionsHistorical

#-----------------
# PART I: Examine skill on yearly basis
#-----------------

# PART Ia: CALCULATE * PREDICTIVE* SKILL FOLLOWING ZWOLINSKI

# FIT THE GAM, OMITTING ONE DATA SET FOR EACH ITERATION OF THE LOOP.  

rp <- numeric(0)   #r value for the predictive performance test

# LOOP OVER EACH OF THE UNIQUE YEARS
UniqueYrs <- unique(LogbookDataWithOceanConditions$CrabYearS)

# NOTE THAT Due to  TIE IN THE WILCOXON TEST, do not use exact solution in the Wilcox test ( "exact = FALSE" in wilcox.test)

for(i in 1:length(UniqueYrs))  {
  print('CALCULATING PREDICTIVE SKILL FOR ONE SURVEY')
  
   # Omit one year's data. 
   crabDataOmit1Yr <-  LogbookDataWithOceanConditions[LogbookDataWithOceanConditions$CrabYearS != UniqueYrs[i],]
  
   # Fit GAM, omitting one year's data. 
   gamO1 <- gam(formula= LnKgsPerPot ~  s(bathymetry,k=3)+s(SoakTime,k=3)+s(dayinseason,k=3)+s(CrabYearS,k=8)+s(bottomArag,k=3)+s(bottomSalt,k=3)+s(bottomTemp,k=3)+s(bottomOxy,k=3)+s(Chla2m,k=3)+s(bottomTemp_4yrLag,k=3)+s(bottomOxy_4yrLag,k=3)+s(Chla2m_4yrLag,k=3)+s(bottomTemp_3yrLag,k=3)+s(bottomOxy_3yrLag,k=3)+te(Lats,Lons,k=3), family=gaussian, data = crabDataOmit1Yr, scale=-1,select=TRUE)   #static + dynamic + lagged vars

      # Predict  one sample year using the newly fit GAM
  
   p <- predict(gamO1 ,newdata= LogbookDataWithOceanConditions[LogbookDataWithOceanConditions$CrabYearS ==
                UniqueYrs[i],] , type="response", exclude="s(CrabYearS)")   # using the average year effect term 
   # Look at true observations of presence./absence at that year.   
   y <- LogbookDataWithOceanConditions$LnKgsPerPot[LogbookDataWithOceanConditions$CrabYearS == UniqueYrs[i]]

 # Calculate correlation coefficient between observed and predicted CPUE
 c <- cor.test(y,p)

rp <- c(rp, c$estimate) 

}  # end loop over sample years


mean(rp)


#------------------
# PART Ib:CALCULATE *EXPLANATORY* SKILL FOLLOWING ZWOLINSKI   
#-----------------
# Note: In wilcox.test, exact set to fALSE only because this was necessary for * Predictive * skill calcs

# DON'T NEED TO FIT THE GAM SEPARATELY FOR EACH ITERATION OF THE LOOP (because we aren't omitting any data)
re <- numeric(0)   #r value for the explanatory performance test

#   Static + Dyn + Lag Ocn - i.e. BEST MODEL
game <- gam(formula= LnKgsPerPot ~  s(bathymetry,k=3)+s(SoakTime,k=3)+s(dayinseason,k=3)+s(CrabYearS,k=9)+s(bottomArag,k=3)+s(bottomSalt,k=3)+s(bottomTemp,k=3)+s(bottomOxy,k=3)+s(Chla2m,k=3)+s(bottomTemp_4yrLag,k=3)+s(bottomOxy_4yrLag,k=3)+s(Chla2m_4yrLag,k=3)+s(bottomTemp_3yrLag,k=3)+s(bottomOxy_3yrLag,k=3)+te(Lats,Lons,k=3), family=gaussian, data = LogbookDataWithOceanConditions, scale=-1,select=TRUE)   


# LOOP OVER EACH OF THE YEARS
for(i in 1:length(UniqueYrs))  {
  print('CALCULATING EXPLANATORY SKILL FOR ONE SURVEY')

  # Predict  one year based on full GAM from above
  pe <- predict(game, newdata = LogbookDataWithOceanConditions[LogbookDataWithOceanConditions$CrabYearS ==
                 UniqueYrs[i],] , type="response", exclude="s(CrabYearS)")   # using the average year effect term 
  
  # Look at true observations of presence./absence at that year. 
  y <- LogbookDataWithOceanConditions$LnKgsPerPot[LogbookDataWithOceanConditions$CrabYearS == UniqueYrs[i]]
  
  # Calculate correlation coefficient between observed and predicted CPUE
  ce <- cor.test(y,pe)
  
  re <- c(re, ce$estimate) 

}  # end loop over years

mean(re)

#------------------------
# PART II: EVALUATE SKILL ON RANDOM SUBSET OF DATA
#------------------------

# Part IIa - Predictive Skill
#------------------------
# leave out a portion of data and calculate predictive skill 
# Only LEARNING SET data is used to fit GAM; then fitted gam is tested against TEST DATA. 
# Note below you can set numReplicates (e.g., 10) and fractionForTestSet (e.g., 0.20). 

# Set parameters that dictate  how to set up the resampling: 

numReplicates <- 10   #  CHOOSE THE NUMBER OF REPLICATE TEST DATA SETS WE WANT
fractionForTestSet <- 0.20  # CHOOSE THE FRACTION OF THE DATA TO RESERVE AS TEST DATA

#-------
numPoints <- length(LogbookDataWithOceanConditions$CrabYearS)

rp20 <- numeric(0)   #r value for the leave-20%-out "learning" performance test

#-------

#  LOOP OVER EACH replicate data set, which will leave out a portion of the data poitns. 
#  Note THAT Due to TIE IN THE WILCOXON TEST, do not use exact solution in the Wilcox test ("exact = FALSE" in wilcox.test)

for(i in 1:numReplicates) {
  print('CALCULATING PREDICTIVE SKILL FOR ONE batch of Test Data')
  
  # Subset data to omit a portion of it 
  reshuffledPointNumbers <- sample(1:numPoints, numPoints, replace=F)  # just a trick to reshuffle the points
  numPointsForTestSet  <- round(fractionForTestSet*numPoints)  
  
  crabDataTestSet <- LogbookDataWithOceanConditions[reshuffledPointNumbers[1:numPointsForTestSet] ,]
  crabDataLearningSet <- LogbookDataWithOceanConditions[ reshuffledPointNumbers[(numPointsForTestSet+1):numPoints] ,]  
  
  # Fit GAM only on the Learning data set
  gaml <- gam(formula= LnKgsPerPot ~  s(bathymetry,k=3)+s(SoakTime,k=3)+s(dayinseason,k=3)+s(CrabYearS,k=9)+s(bottomArag,k=3)+s(bottomSalt,k=3)+s(bottomTemp,k=3)+s(bottomOxy,k=3)+s(Chla2m,k=3)+s(bottomTemp_4yrLag,k=3)+s(bottomOxy_4yrLag,k=3)+s(Chla2m_4yrLag,k=3)+s(bottomTemp_3yrLag,k=3)+s(bottomOxy_3yrLag,k=3)+te(Lats,Lons,k=3), family=gaussian, data = crabDataLearningSet, scale=-1,select=TRUE)   #static + dynamic + lagged vars
  
  # Predict TEST DATA , from the  GAM that you just trained with the LEARNING SET data
  p20<- predict(gaml, newdata = crabDataTestSet, type="response", exclude="s(CrabYearS)")   # using the average year effect term 
  
  #Look at true observations of CPUE in that test data: 
  y <- crabDataTestSet$LnKgsPerPot
  
  # Calculate correlation coefficient between the observed and predicted
  cp20 <- cor.test(y,p20)
  
  rp20 <- c(rp20 , cp20$estimate)  
  
}  # end loop over say the 10 replicates. 


mean(rp20)

#-----
# Part IIb - Explanatory Skill (leaving out a portion of your data):

re20 <- numeric(0)   #r value for the leave-20%-out "learning" performance test

# LOOP OVER EACH OF THE DATA SUBSETS  # previously was for each of the FOUR SAMPLE DAYS NORTH OF 43 DEGREES (AND IN OUR MODEL DOMAIN) IN 2009. 
for(i in 1:numReplicates)  {
  print('CALCULATING EXPLANATORY SKILL FOR % omitted data')
  
  # Use best fitting GaM from above to predict our TEST set 
  crabDataTestSet <- LogbookDataWithOceanConditions[reshuffledPointNumbers[1:numPointsForTestSet] ,]
  pe20<- predict(game, newdata = crabDataTestSet, type="response", exclude="s(CrabYearS)")   # using the average year effect term 
  
  # Look at true observations of CPUE for the TEST set 
  y <- crabDataTestSet$LnKgsPerPot
  
  # Calculate correlation coefficient between observed and predicted CPUE
  ce20 <- cor.test(y,pe20)

  re20 <- c(re20, ce20$estimate)   

}  # end loop over sample days


mean(re20)


