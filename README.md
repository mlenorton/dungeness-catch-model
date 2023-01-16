# dungeness-catch-model
Match ocean conditions to Dungeness crab catch per unit effort (CPUE) to model catch rates  

These scripts support Norton et al. (in prep) [doi:XXX]. The parent script for fitting the generalized additive models (GAMs) is "runCrabLogbooks.R". This calls several subroutines, included here. Additionally, there are scripts to test the performance of the GAMs: 1) TestGAMPerformance_PredVsExplSkill_Zwolinski.R - examines the predictive and explanatory skill of the GAMs in predict mode; and 2) match_allFishingData_forecastYears.R and forecastGAM_fishingBehavSelect.R - when run together, they can generate forecast maps for a specific year when different fishing behaviors are selected. Finally, crabLogbooksSpatialExploration.R is used to evaluate spatial autocorrelation of the GAM predictor variables. The function "point.in.SpatialPolygons" is publicly available: https://rdrr.io/cran/prevR/src/R/point.in.SpatialPolygons.prevR.r. We relied heavily on the model performance testing outlined by Zwolinski et al. (2011).  

Zwolinski, J.P., Emmett, R.L., and Demer, D.A. 2011. Predicting habitat to optimize sampling of Pacific sardine (Sardinops sagax). ICES Journal of Marine Science, 68(5): 867-879. doi:10.1093/icesjms/fsr038  

We used R version 3.6.3 for this project. The following R libraries (version) were also used:  
ape (5.4)  
basicTrendline  (2.0.3)  
dplyr (1.0.0)  
geoR (1.8.1)  
geosphere (1.5.10)  
ggplot2 (3.3.1)  
gstat (2.0.8)  
lubridate (1.7.9)  
maptools (1.0.1)  
mgcv (1.8.31)  
ncdf4 (1.17)  
nlme (3.1.148)  
PBSmapping (2.72.1)  
RColorBrewer (1.1.2)   
reshape2 (1.4.4)  
rgdal (1.5.12)  
sp (1.4.2)  
stringr (1.4.0)   
tidyr (1.1.0)  
tidyverse (1.3.0)  
