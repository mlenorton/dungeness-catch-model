AssignDaysSinceStartDatePerArea <- function(WDFW1YearCrabLbsAndOceanConditions, crabYear,forecastYear)
{
  #------
  # Primary Author: Isaac Kaplan (Feb 20 2019) 
  # isaac.kaplan at noaa.gov
  # Secondary Author: Emily Norton
  # emilyln@uw.edu
  
  # Assigns the start date for which each crab area opened (in Washington) based on the fishing location, and then calculates day in season for that fishing observation based on that start date
  #  Arguments:
  #       WDFW1YearCrabLbsAndOceanConditions: single year of Washington logbook data
  #       crabYear:  would be for instance 2017 for fishing year autumn 2107--spring 2018.
  #       forecastYear: the year you'd like to use this forecast for (usually same as the "crabYear" but wouldn't have to be)
  #
  #   Returns: WDFW1YearCrabLbsAndOceanConditions, each row has logbook record (including cath rate per trap) and associated opening date and day in season
  #------

  
library(tidyverse)
library(lubridate)
library(stringr)
library(rgdal)
library(maptools)

# adding column of crabYear to dataframe. crabYear is 2008 for fishing year that would start in Dec 2008 for instance. 
  
  WDFW1YearCrabLbsAndOceanConditions$crabYear <- as.character(crabYear)
  WDFW1YearCrabLbsAndOceanConditions$forecastYear <- forecastYear
  

############################################################
#                                                          #
#                      Opening dates : from Kate Richerson code      #
#                                                          #
############################################################

wa_dates_post05<-read_csv("../Data/SeasonOpenings/WAestimatedopenings.csv") %>% 
  gather(season,open_date,-`Management Areas`) %>% 
  rename(catch_area=`Management Areas`) %>% 
  mutate(open_date=str_replace(open_date,"Open:","")) %>% 
  mutate(open_date=as.Date(open_date,format=" %b %d, %Y")) %>% 
  mutate(catch_area=tolower(catch_area)) %>%  
  filter(catch_area!="area_3") %>% #since area 2 and 3 have always opened at the same time, combine them into one opening date
  mutate(catch_area=replace(catch_area, catch_area=="area_2", "area_2_3")) 

latcuts2<-c(-Inf, 46+28/60, 47+40.50/60, 48+26/60, Inf) 
latnames2<-c("area_1", "area_2_3", "area_4","north")

# NOTE WILLAPA OPENS WITH AREA 1, HENCE IT IS BOXED SEPARATELY HERE: 
#draw a rough polygon around willapa bay
#cape showalter 46.741677, -124.094060
#ledbetter point 46.659584, -124.061027 (approx)
willapa_coords = matrix(c(46.741677, -124.094060, #cape showalter
                          46.659584, -124.061027, #ledbetter
                          46.325372, -124.034934, #southwest corner approx
                          46.336158, -123.591018, #southeast corner approx
                          46.792201, -123.666549), #northeast cornder approx
                        ncol = 2, byrow = TRUE)

willapa_poly<-Polygon(willapa_coords)
willapa_sp_poly<- SpatialPolygons(list(Polygons(list(willapa_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))


WDFW1YearCrabLbsAndOceanConditions<-  WDFW1YearCrabLbsAndOceanConditions %>%
  mutate(catch_area=as.character(cut(Lats,breaks=latcuts2,labels=latnames2))) %>%
  mutate(catch_area=ifelse(point.in.SpatialPolygons(Lats, Lons, willapa_sp_poly)==TRUE,"area_1",catch_area))  %>%  #ADD IN WILLAPA HARBOR TO AREA 1--these coords make an approxiamte box around willapa
  left_join(wa_dates_post05, c("crabYear"="season","catch_area"))  #%>% 

#-------------
# tests
print("open date")
print(WDFW1YearCrabLbsAndOceanConditions$open_date[1])
as.POSIXct(as.Date(WDFW1YearCrabLbsAndOceanConditions$open_date,format='%Y-%m-%d'))
#---------------

WDFW1YearCrabLbsAndOceanConditions$dayinseason <-as.numeric(difftime(as.POSIXct(as.Date(WDFW1YearCrabLbsAndOceanConditions$Set.Date,format='%m/%d/%Y')),as.POSIXct(as.Date(WDFW1YearCrabLbsAndOceanConditions$open_date,format='%Y-%m-%d')),units="days"))
WDFW1YearCrabLbsAndOceanConditions <- WDFW1YearCrabLbsAndOceanConditions %>% filter(dayinseason>=0)


return(WDFW1YearCrabLbsAndOceanConditions)

}
