#################################################################
#################################################################

### DIURNAL PATTERNS OF VISUALLY OBSERVED HUMPBACK WHALES FROM A LAND STATION OFF MAUI, HAWAIʻI

# Written by Anke Kügler
# Last updated by Anke Kügler, 06/10/2022

# This script analyses sightings of humpback whales obtained during land-based visual surveys off Olowalu, Maui, Hawaii
# during the 2017-2020 whale seasons to investigate diurnal patterns of occurrence relative to three bathymetric features.
# Latitude/Longitude coordinates of pods were calculated in the Program PYTHAGORAS (Gailey & Ortega-Ortiz, 2000) after being
# fixed with a surveyor's theodolite ('fixes'). Some whales were not able to be located with the theodolite, and compass bearing
# and reticle distance were determined with binoculars and subsequently converted to Latitude/Longitude coordinates in this script.

# Data were obtained during hourly 30-minute scans between ca. 8:00 HST and 14:00 HST per survey day, weather permitting.

# For the scope of this study, only adult pods (pods not containing a calf) are analysed due to the large distance cutoff
# and the decreasing detection probability of calves with distance to the land station.

# Input data: excel spreadsheet containing fix and non-fix sighting data per year
#             bathymetry data for Maui

#Generalized Additive Models (GAMs) are fit for three bathymetric metrics (distance to shore, depth, distance to -200 m isobath) and 
#   day of the season
#   hour of the day
#   season (representing respective whale seasons from December through April)

#Functions:
# pkgTest() is a function that tests if a R package already exists on a system and if not, installs and loads it
# pkgTest.github() does the same as pkgTest() for packages that are only available through github
# gam_plots_params() and gam_paraplots_params() define settings for the GAM effect plots

#Note: all calls to save data/plots to file have been commented out for public distribution.


################################################################
################################################################

#clear all objects just in case
rm(list=ls())

################################################################
#Load required packages

pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    
    if(!require(x,character.only = TRUE)) {stop("Package not found")
    }else {
      require(x, character.only = T)
  }}
}

pkgTest.github <- function(x){
  
  if(!require(x,character.only = T))
  {
    pkgTest('devtools')
    pkgTest('remotes')
    
    if(x=='bangarang'){remotes::install_github("ericmkeen/bangarang")}
    if(x=='ungeviz'){remotes::install_github("wilkelab/ungeviz")}
    
    if(!require(x,character.only = T)) stop("Package not found")
    else {
      require(x, character.only = T, quietly=T)
    }
  }}    

# ---

pkgTest("readxl") #to read in data from excel files

pkgTest.github('bangarang') #Calculate Lat/Long for non-fix sightings
pkgTest('marmap') #get depths and distance to isobaths for sightings

pkgTest('dplyr') #to do all kinds of data summary (e.g. for effort, descriptive statistics)
  options(dplyr.summarise.inform = FALSE) # Suppress dplyr::summarise info because they're heckin' annoying

pkgTest('chron') #date and time manipulations
pkgTest('lubridate') #date and time manipulations

pkgTest('mgcv') #to fit and analyse GAMs

pkgTest('ggplot2') #create pretty plots
pkgTest.github('ungeviz') #draw point as horizontal line in GAM parameter plot
pkgTest('patchwork') #combine multiple ggplots into one plot


################################################################
#Functions

#GAM plots settings

gam_plots_params<-function(season=F, diel=F, signif=T){
  list(
    if(signif)geom_ribbon(aes(ymin=fit-se, ymax=fit+se, group=1), 
                fill = "grey75", alpha = 0.6),
    if(signif)geom_line(color='black'),
    if(!signif)geom_ribbon(aes(ymin=fit-se, ymax=fit+se, group=1), 
                          fill = "lightcoral", alpha = 0.3),
    if(!signif)geom_line(color='#ba6a68'),
    theme_bw(),
    if(signif)theme(axis.text = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust=0.5),
          panel.border=element_rect(color="black", size=0.75), 
          panel.grid = element_blank(),
          axis.line.x.top = element_line(color="black", size=0.75)),
    if(!signif)theme(axis.text = element_text(size = 14, color = "#ba6a68"),
                    axis.title = element_text(size = 14, color = "#ba6a68"),
                    axis.text.y = element_text(size = 14, color = "#ba6a68", angle = 90, hjust=0.5),
                    panel.border=element_rect(color="#ba6a68", size=0.75), 
                    panel.grid = element_blank(),
                    axis.line.x.top = element_line(color="#ba6a68", size=0.75),
                    axis.ticks= element_line(color="#ba6a68")),
    if(season)scale_x_continuous(limits=c(0,145), breaks=c(0,20,40,60,80,100,120, 140), labels=c('0','20','40','60','80','100','120','140')),
    if(diel)  scale_x_continuous(breaks=seq(from=8, to=14, by=1)/24, labels=c('8', '9', '10', '11', '12', '13', '14'))
  )
}

gam_paraplots_params<-function(diel=F, connect=F, signif=T){
  list(
    if(signif)geom_errorbar(aes(ymin=y-se, ymax=y+se), 
                  width=0.75, position=position_dodge(width=0.9), stat = "identity", color='grey65'),
    if(!signif)geom_errorbar(aes(ymin=y-se, ymax=y+se), 
                            width=0.75, position=position_dodge(width=0.9), stat = "identity", color='lightcoral', alpha=0.3),
    if(connect)geom_line(color='grey50', linetype='dashed', group=1),
    if(signif)geom_hpline(stat='identity', size=1, width=0.75, color='black'),
    if(!signif)geom_hpline(stat='identity', size=1, width=0.75, color='#ba6a68'),
    theme_bw(),
    if(signif)theme(axis.text = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black", angle = 90, hjust=0.5),
          panel.border=element_rect(color="black", size=0.75), 
          panel.grid = element_blank(),
          axis.line.x.top = element_line(color="black", size=0.75)),
    if(!signif)theme(axis.text = element_text(size = 14, color = "#ba6a68"),
          axis.title = element_text(size = 14, color = "#ba6a68"),
          axis.text.y = element_text(size = 14, color = "#ba6a68", angle = 90, hjust=0.5),
          panel.border=element_rect(color="#ba6a68", size=0.75), 
          panel.grid = element_blank(),
          axis.line.x.top = element_line(color="#ba6a68", size=0.75),
          axis.ticks= element_line(color="#ba6a68")),
  if(diel)scale_x_discrete(labels=c('8:00', '9:00', '10:00', '11:00', '12:00', '13:00', '14:00'))
   )
}

################################################################


################################################################
#######################  Analysis  #############################
################################################################

#### 1) GET DATA ----

#Bathymetry data for the study area

## not necessary to run; raster csv file provided with data ##

#read in bathymetry data saved as text from ArcGIS and convert to raster
#pkgTest("raster")
#nc.brick <- brick(file.choose())
#nc.df <- as.data.frame(nc.brick[[1]], xy=T)
#names(nc.df)[3]<-'z'
#write.csv(nc.df, 'Land/maui_bathy_5m.csv', row.names=F)

###

#import as bathymetry with marmap (long loading time; .Rda file of bathymetry data provided for faster loading)

#maui.bathy<-read.bathy('LandObservations/maui_bathy_5m.csv', header=T)
#saveRDS(maui.bathy, file='LandObservations/maui_bathy_5m.Rda')
maui.bathy<-readRDS('LandObservations/maui_bathy_5m.Rda')


#---

#Land-based data from Olowalu station

#get data folder path
#dat.dir<-choose_dir()
dat.dir<-'LandObservations/'

###

olowalu.2017<-as.data.frame(read_excel(paste(dat.dir, "Kuegler_Maui_LandObservations_2017-2020.xlsx", sep=''), sheet = "Fixes_2017"))
olowalu.2018<-as.data.frame(read_excel(paste(dat.dir, "Kuegler_Maui_LandObservations_2017-2020.xlsx", sep=''), sheet = "Fixes_2018"))
olowalu.2019<-as.data.frame(read_excel(paste(dat.dir, "Kuegler_Maui_LandObservations_2017-2020.xlsx", sep=''), sheet = "Fixes_2019"))
olowalu.2020<-as.data.frame(read_excel(paste(dat.dir, "Kuegler_Maui_LandObservations_2017-2020.xlsx", sep=''), sheet = "Fixes_2020"))

#add year
olowalu.2017$year<-"2017"
olowalu.2018$year<-"2018"
olowalu.2019$year<-"2019"
olowalu.2020$year<-"2020"

#get column indices for subsetting
i.2017<-which(names(olowalu.2017) %in% c('Date', 'Time', 'year', 'Scan', 'Whale', 'Lat', 'Long', 'Distance'))
i.2018<-which(names(olowalu.2018) %in% c('Date', 'Time', 'year', 'Scan', 'Whale', 'Lat', 'Long', 'Distance'))
i.2019<-which(names(olowalu.2019) %in% c('Date', 'Time', 'year', 'Scan', 'Whale', 'Lat', 'Long', 'Distance'))
i.2020<-which(names(olowalu.2020) %in% c('Date', 'Time', 'year', 'Scan', 'Whale', 'Lat', 'Long', 'Distance'))

#combine into one data frame
olowalu.land<-rbind(olowalu.2017[,i.2017], olowalu.2018[,i.2018], olowalu.2019[,i.2019], olowalu.2020[,i.2020])

#---

#non-fix data from Olowalu station

olowalu.nonfix.2017<-as.data.frame(read_excel(paste(dat.dir, "Kuegler_Maui_LandObservations_2017-2020.xlsx", sep=''), sheet = "NonFixes_2017"))
olowalu.nonfix.2018<-as.data.frame(read_excel(paste(dat.dir, "Kuegler_Maui_LandObservations_2017-2020.xlsx", sep=''), sheet = "NonFixes_2018"))
olowalu.nonfix.2019<-as.data.frame(read_excel(paste(dat.dir, "Kuegler_Maui_LandObservations_2017-2020.xlsx", sep=''), sheet = "NonFixes_2019"))
olowalu.nonfix.2020<-as.data.frame(read_excel(paste(dat.dir, "Kuegler_Maui_LandObservations_2017-2020.xlsx", sep=''), sheet = "NonFixes_2020"))

#add year
olowalu.nonfix.2017$year<-"2017"
olowalu.nonfix.2018$year<-"2018"
olowalu.nonfix.2019$year<-"2019"
#olowalu.nonfix.2020$year<-"2020" #Note: olowalu.nonfix.2020 is empty because there were no nonfixes, so trying to add year creates an error

#get column indices for subsetting
i.2017.nonfix<-which(names(olowalu.nonfix.2017) %in% c('Date', 'Time', 'year', 'Scan', 'Whale', 'Reticle', 'Bearing'))
i.2018.nonfix<-which(names(olowalu.nonfix.2018) %in% c('Date', 'Time', 'year', 'Scan', 'Whale', 'Reticle', 'Bearing'))
i.2019.nonfix<-which(names(olowalu.nonfix.2019) %in% c('Date', 'Time', 'year', 'Scan', 'Whale', 'Reticle', 'Bearing'))
i.2020.nonfix<-which(names(olowalu.nonfix.2020) %in% c('Date', 'Time', 'year', 'Scan', 'Whale', 'Reticle', 'Bearing'))

#combine into one data frame
olowalu.land.nonfix<-rbind(olowalu.nonfix.2017[,i.2017.nonfix], olowalu.nonfix.2018[,i.2018.nonfix], olowalu.nonfix.2019[,i.2019.nonfix], olowalu.nonfix.2020[,i.2020.nonfix])

#---

#some house-keeping: 

#exclude off-effort fixes and non-pod fixes
olowalu.land<-olowalu.land[olowalu.land$Scan!="off-effort",] 
olowalu.land<-olowalu.land[olowalu.land$Whale!="other",] 

#convert scan factor to numeric 
olowalu.land$Scan<-as.numeric(olowalu.land$Scan)

#adjust sign for Longitude
olowalu.land$Long<-olowalu.land$Long*(-1) 

#exclude whales that do not have GPS coordinates
olowalu.land<-olowalu.land[complete.cases(olowalu.land),] 

#convert distance to numeric
olowalu.land$Distance<-as.numeric(olowalu.land$Distance)

#fix erroneous Time timestamps (artifact when reading in from Excel spreadsheet)
olowalu.land$Time<-times(strftime(olowalu.land$Time, format="%H:%M:%S", tz="UTC"))
olowalu.land.nonfix$Time<-times(strftime(olowalu.land.nonfix$Time, format="%H:%M:%S", tz="UTC"))

#standardize pod class names/fix spelling
olowalu.land$Whale[olowalu.land$Whale=='whale']<-'Whale' 
olowalu.land$Whale[olowalu.land$Whale=='Comp'| olowalu.land$Whale=='Com Pod' | olowalu.land$Whale=='Comp pod']<-'Comp Pod' 
olowalu.land$Whale[olowalu.land$Whale=='Mc']<-'M/c' 

#### ----


# 2) CALCULATE LAT/LONG FOR NON-FIX SIGHTINGS ---

#initialize columns
olowalu.land.nonfix$Lat<-NA
olowalu.land.nonfix$Long<-NA
olowalu.land.nonfix$Distance<-NA

col.reticle<-which(names(olowalu.land.nonfix)=='Reticle')
col.bearing<-which(names(olowalu.land.nonfix)=='Bearing')
col.Lat<-which(names(olowalu.land.nonfix)=='Lat')
col.Long<-which(names(olowalu.land.nonfix)=='Long')
col.Dist<-which(names(olowalu.land.nonfix)=='Distance')

#function to calculate Lat/Long
get_whalePos<-function(bear, reticle){whalemap(X=-156.6315,Y=20.829,bearing=bear,reticle=reticle,eye.height=83.7, mag = 9.51, vessel.hdg=50, deg.per.ret = 0.279, toplot=F)}

#Loop through Nonfix DF and calculate for each sighting
for(j in 1:length(olowalu.land.nonfix$Date)){
  
  bear<-olowalu.land.nonfix[j,col.bearing]
  reticle<-olowalu.land.nonfix[j,col.reticle]
  
  position<-get_whalePos(bear, reticle)
  
  olowalu.land.nonfix[j,col.Long]<-position[[1]]
  olowalu.land.nonfix[j,col.Lat]<-position[[2]]
  olowalu.land.nonfix[j,col.Dist]<-position[[3]]
  
}

#convert distance from km to m
olowalu.land.nonfix$Distance<-olowalu.land.nonfix$Distance*1000

#include measuring method
olowalu.land$method<-'fix'
olowalu.land.nonfix$method<-'nonfix'

#---

#combine nonfix and fix data into single DF
olowalu.land<-rbind(olowalu.land, olowalu.land.nonfix[,-c(col.reticle,col.bearing)])
olowalu.land<-olowalu.land[order(olowalu.land$Date, olowalu.land$Time),] #order by date and time
rownames(olowalu.land)<-NULL #reset rownumbers

#### ----

# 2) SOME MORE HOUSE-KEEPING AND ADDITIONAL COVARIATES ---

#calculate hour from time
olowalu.land$Hour<-as.POSIXct(as.character(olowalu.land$Time), format = "%H:%M:%S")
olowalu.land$Hour<-round_date(olowalu.land$Hour, unit="hour")
olowalu.land$Hour<-times(strftime(olowalu.land$Hour, "%H:%M:%S"))

#---

#convert scan and year to factor (important for GAM later)
olowalu.land$Scan<-factor(olowalu.land$Scan, ordered=T)
olowalu.land$year<-as.factor(olowalu.land$year)
#olowalu.land$year<-factor(olowalu.land$year, ordered=T, levels = c("2017", "2018", "2019", "2020"))

#create a calf-pod factor
olowalu.land$calf.pod<-apply(olowalu.land['Whale'],1, function(x) as.integer(x=='M/c' | x=='M/c/E')) 

#---

#exclude days with only one scan

survey.summary<-as.data.frame(olowalu.land%>%
                      group_by(Date)%>%
                      summarize(n.Scans=length(unique(Scan)))
                      )

i<-which(survey.summary$n.Scans==1)

olowalu.land<-olowalu.land[!olowalu.land$Date %in% survey.summary$Date[i],]

#get dates with only 1 scan
excluded.days<-as.Date(survey.summary$Date[i])

#---

#summarize effort before subsetting

survey.effort.month<-as.data.frame(olowalu.land %>%
                                     group_by(year, month(Date)) %>%
                                     summarise(
                                       n=length(unique(Date))))

survey.effort.h<-as.data.frame(olowalu.land %>%
                                 group_by(year, Hour) %>%
                                 summarise(
                                   n=length(unique(Date))))

#write.csv(survey.effort.month, 'Results/Maui_LandObservations_effort_month.csv', row.names=F)
#write.csv(survey.effort.h, 'Results/Maui_LandObservations_effort_hour.csv', row.names=F)


#---

#reduce fixes to a distance of 10 km from land station

olowalu.land<-olowalu.land[olowalu.land$Distance<10000,]

#take out calf pods, we're only intersted in non-calf pods for this study

olowalu.land<-olowalu.land[olowalu.land$calf.pod==0,]

# ---

#Get associated bathymetric metrics for each sighting: depth, distance to shore, distance to -200 isobath

#Note: very long processing times

#a) Distance to Shore (shore=isobath 0) ---

dist_shore<-round(dist2isobath(maui.bathy, olowalu.land[,c('Long', 'Lat')], isobath=0, locator=F))
olowalu.land$dist_shore<-dist_shore[,1]

#b) Depth ---

depths<-get.depth(maui.bathy, olowalu.land[,c('Long', 'Lat')], locator=F)
olowalu.land$depth<-depths[,3]

#c) Distance to -200m isobath ---

dist_200<-round(dist2isobath(maui.bathy, olowalu.land[,c('Long', 'Lat')], isobath=-200, locator=F))
olowalu.land$dist_200<-dist_200[,1]

# ---

#reference whale season relative to December 1st

olowalu.land$julian<-as.numeric(format(as.Date(olowalu.land$Date), "%j"))

for (n in 1:nrow(olowalu.land)){
  
  season.start<-as.numeric(format(as.Date(paste(olowalu.land$year[n],'-12-01', sep='')),'%j'))
  season.start.diff<-ifelse(leap_year(olowalu.land$Date[n]),366-as.numeric(format(as.Date(paste(olowalu.land$year[n],'-12-01', sep='')),'%j')),365-as.numeric(format(as.Date(paste(olowalu.land$year[n],'-12-01', sep='')),'%j')))
  
  if (olowalu.land$julian[n]>=season.start){
    olowalu.land$day.season[n]<-olowalu.land$julian[n]-season.start
  }else{
    olowalu.land$day.season[n]<-olowalu.land$julian[n]+season.start.diff
  }
  
  n<-n+1
}	


### ---


# 3) SUMMARIZE DATA (descriptive Statistics) ---

#descriptive statistics by year and across years
olowalu.land.summary<-as.data.frame(olowalu.land %>%
                             group_by(year) %>%
                             summarise(n=n(),
                                       median.dist_shore=median(dist_shore, na.rm=T),
                                       IQR.dist_shore=IQR(dist_shore, na.rm=T),
                                       median.depth=median(depth, na.rm=T),
                                       IQR.depth=IQR(depth, na.rm=T),
                                       median.dist_200=median(dist_200, na.rm=T),
                                       IQR.dist_200=IQR(dist_200, na.rm=T)) %>%
                               ungroup() %>%
                               mutate_at(vars(year), funs(as.character(.))) %>%
                               bind_rows(summarise(year="2017-2020",olowalu.land,
                                                   n=n(),
                                                   median.dist_shore=median(dist_shore, na.rm=T),
                                                   IQR.dist_shore=IQR(dist_shore, na.rm=T),
                                                   median.depth=median(depth, na.rm=T),
                                                   IQR.depth=IQR(depth, na.rm=T),
                                                   median.dist_200=median(dist_200, na.rm=T),
                                                   IQR.dist_200=IQR(dist_200, na.rm=T)))
                              )

#write.csv(olowalu.land.summary, 'Results/Maui_LandObservations_summary.csv', row.names=F)


##correct for scans with no whale sightings (create "NA" sightings) to get accurate n averages per scan

n.whales<-as.data.frame(xtabs(~Date+Hour, data=olowalu.land))
no.whales<-n.whales[which(n.whales$Freq==0),]
no.whales$timestamp<-as.POSIXct(paste(no.whales$Date, no.whales$Hour, sep=' '))

#find scans with no effort
effort<-read.csv(paste(dat.dir,'Kuegler_Maui_LandObservations_2017-2020_effort.csv', sep=''), header=T)
names(effort)<-c('Date', '08:00:00', '09:00:00', '10:00:00', '11:00:00', '12:00:00', '13:00:00', '14:00:00')
effort$Date<-as.Date(effort$Date)
effort<-reshape2::melt(effort, id.vars='Date', variable.name='Hour', value.name='effort')
effort$Hour<-times(strftime(as.POSIXct(as.character(effort$Hour), format = "%H:%M:%S"), "%H:%M:%S"))
no.effort<-effort[which(effort$effort==0),]
no.effort$timestamp<-as.POSIXct(paste(no.effort$Date, no.effort$Hour, sep=' '))

#which scans (effort) existed but had zero sightings?
no.whales.but.effort<-no.whales[!no.whales$timestamp %in% no.effort$timestamp,]
no.whales.but.effort$year<-year(no.whales.but.effort$Date)
#correct year for December
no.whales.but.effort$year[which(as.numeric(format(as.Date(no.whales.but.effort$Date), "%j"))>334)]<-no.whales.but.effort$year[which(as.numeric(format(as.Date(no.whales.but.effort$Date), "%j"))>334)]+1

#convert year to factor
no.whales.but.effort$year<-factor(no.whales.but.effort$year, ordered=T)

#create dummy NA's for the other variables
xx<-c('Time', 'Scan', 'Whale', 'Lat', 'Long', 'Distance', 'method', 'depth', 'calf.pod', 'depth', 'dist_shore', 'dist_shore', 'dist_200', 'julian', 'day.season')
no.whales.but.effort<-cbind(no.whales.but.effort, setNames( lapply(xx, function(x) x=NA), xx) )
no.whales.but.effort<-no.whales.but.effort[,c('Date','Time', 'Scan', 'Whale', 'Lat', 'Long', 'Distance', 'year', 'method', 'Hour', 'calf.pod',
                                              'depth', 'dist_shore', 'dist_200', 'julian', 'day.season')]

olowalu.land.complete<-rbind(olowalu.land,no.whales.but.effort)

#descriptive statistics by hour
olowalu.land.summary.byhour<-as.data.frame(olowalu.land.complete %>%
                                             group_by(year, Hour)%>%
                                             mutate(id = factor(Date)) %>%
                                             summarise(n=n(),
                                                       n.median=median(table(Date)),
                                                       n.IQR=IQR(table(Date)),
                                                       median.dist_shore=median(dist_shore, na.rm=T),
                                                       IQR.dist_shore=IQR(dist_shore, na.rm=T),
                                                       median.depth=median(depth, na.rm=T),
                                                       IQR.depth=IQR(depth, na.rm=T),
                                                       median.dist_200=median(dist_200, na.rm=T),
                                                       IQR.dist_200=IQR(dist_200, na.rm=T)) %>%
                                             ungroup() %>%
                                             mutate_at(vars(year), funs(as.character(.))) %>%
                                             bind_rows( summarise(year="2017-2020",olowalu.land.complete[!(olowalu.land.complete$year=='2018'&olowalu.land.complete$Hour=='08:00:00'),]%>%group_by(Hour),
                                                                 n=n(),
                                                                 n.median=median(table(Date)),
                                                                 n.IQR=IQR(table(Date)),
                                                                 median.dist_shore=median(dist_shore, na.rm=T),
                                                                 IQR.dist_shore=IQR(dist_shore, na.rm=T),
                                                                 median.depth=median(depth, na.rm=T),
                                                                 IQR.depth=IQR(depth, na.rm=T),
                                                                 median.dist_200=median(dist_200, na.rm=T),
                                                                 IQR.dist_200=IQR(dist_200, na.rm=T))) 
                                              )

#write.csv(olowalu.land.summary.byhour, 'Results/Maui_LandObservations_summary_hour.csv', row.names=F)


# 4) STATISTICAL ANALYSIS (Generalized Additive Models) ----

gam.dist_shore<-gam(dist_shore ~ s(day.season) + s(Hour, k=7) + year, data=olowalu.land, method='REML', select=T)
gam.depth<-gam(depth ~ s(day.season) + s(Hour, k=7) + year, data=olowalu.land, method='REML', select=T)
gam.dist_200<-gam(dist_200 ~ s(day.season) + s(Hour, k=7) + year, data=olowalu.land, method='REML', select=T)

#gam.dist_shore.full<-gam(dist_shore ~ s(day.season) + ti(day.season, Hour) + s(Hour, k=7) + year, data=olowalu.land, method='REML', select=T)
#gam.depth.full<-gam(depth ~ s(day.season) + ti(day.season, Hour) + s(Hour, k=7) + year, data=olowalu.land, method='REML', select=T)
#gam.dist_200.full<-gam(dist_200 ~ s(day.season) + ti(day.season, Hour) + s(Hour, k=7) + year, data=olowalu.land, method='REML', select=T)

#no significant ti interactions:
#dist_shore: edf=0.04823, F=0.002, p=0.692851
#depth: edf=0.009323, F=0.000 p=0.670484
#dist_200: edf=0.95321, F=0.120, p=0.1284

#get signifcances
anova.dist_shore<-anova(gam.dist_shore)
anova.depth<-anova(gam.depth)
anova.dist_200<-anova(gam.dist_200)

# ---

#Create effect plots. Because plot.gam() doesn't allow for a lot of adjustments of the look of the plots,
#we're going to create them from scratch.

# (a) - pull out all smoothers (seasonal and diel)

sm.shore <- plot(gam.dist_shore, pages=1, scale=0)  # plot.gam returns a list of n elements, one per plot
sm.depth <- plot(gam.depth, pages=1, scale=0)  
sm.200 <- plot(gam.dist_200, pages=1, scale=0)  

#get x, fit, and se for each plot and save to DF
get_GAMsmoothers<-function(smoothers,n){
  sm.coef<-smoothers[[n]]
  sm.fit<-as.data.frame(sm.coef$fit)
  names(sm.fit)<-'fit'
  sm.fit$x<-as.vector(sm.coef$x)
  sm.fit$se<-as.vector(sm.coef$se)
  return(sm.fit)
}

# Distance to shore

sm.shore.seasonal<-get_GAMsmoothers(sm.shore,1)
sm.shore.diel<-get_GAMsmoothers(sm.shore,2)

# Depth

sm.depth.seasonal<-get_GAMsmoothers(sm.depth,1)
sm.depth.diel<-get_GAMsmoothers(sm.depth,2)

# Distance to -200 m

sm.200.seasonal<-get_GAMsmoothers(sm.200,1)
sm.200.diel<-get_GAMsmoothers(sm.200,2)


###---

# (b) - get parametric terms (season)

para.shore <- termplot(gam.dist_shore, se = TRUE, plot = FALSE)
para.depth <- termplot(gam.depth, se = TRUE, plot = FALSE)
para.200 <- termplot(gam.dist_200, se = TRUE, plot = FALSE)

###---

# (c) - create individual plots and save as objects

# Distance to Shore

gam.plot.shore.seasonal<-ggplot(data=sm.shore.seasonal, aes(x=x, y=fit)) +
  gam_plots_params(season=T, signif=ifelse(anova.dist_shore$s.table[1,4]<=0.05, TRUE, FALSE)) +
  labs(x='Day of the season', y="Partial effect for distance to shore") +
  scale_y_continuous(breaks=c(0, 500, 1000,1500), labels=c('0', '500','1000', ''))

gam.plot.shore.hour<-ggplot(data=sm.shore.diel, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, signif=ifelse(anova.dist_shore$s.table[2,4]<=0.05, TRUE, FALSE)) +
  labs(x='Hour', y="Partial effect for distance to shore") +
  scale_y_continuous(breaks=c(-300, -200, -100, 0, 100, 200), labels=c('', '-200','', '0', '', '200')) 

gam.paraplot.shore.seasons<-ggplot(data=para.shore$year, aes(x=x, y=y)) +
  gam_paraplots_params(signif=ifelse(anova.dist_shore$pTerms.table[3]<=0.05, TRUE, FALSE)) +
  labs(x='', y="Partial effect for distance to shore")

# Depth

gam.plot.depth.seasonal<-ggplot(data=sm.depth.seasonal, aes(x=x, y=fit)) +
  gam_plots_params(season=T, signif=ifelse(anova.depth$s.table[1,4]<=0.05, TRUE, FALSE)) +
  labs(x='Day of the season', y="Partial effect for depth")

gam.plot.depth.hour<-ggplot(data=sm.depth.diel, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, signif=ifelse(anova.depth$s.table[2,4]<=0.05, TRUE, FALSE)) +
  labs(x='Hour', y="Partial effect for depth")

gam.paraplot.depth.seasons<-ggplot(data=para.depth$year, aes(x=x, y=y)) +
  gam_paraplots_params(signif=ifelse(anova.depth$pTerms.table[3]<=0.05, TRUE, FALSE)) +
  labs(x='', y="Partial effect for depth") 

# Distance to -200 m isobath

gam.plot.200.seasonal<-ggplot(data=sm.200.seasonal, aes(x=x, y=fit)) +
  gam_plots_params(season=T, signif=ifelse(anova.dist_200$s.table[1,4]<=0.05, TRUE, FALSE)) +
  labs(x='Day of the season', y="Partial effect for distance to -200 isobath")


gam.plot.200.hour<-ggplot(data=sm.200.diel, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, signif=ifelse(anova.dist_200$s.table[2,4]<=0.05, TRUE, FALSE)) +
  labs(x='Hour', y="Partial effect for distance to -200 isobath")

gam.paraplot.200.seasons<-ggplot(data=para.200$year, aes(x=x, y=y)) +
  gam_paraplots_params(signif=ifelse(anova.dist_200$pTerms.table[3]<=0.05, TRUE, FALSE)) +
  labs(x='', y="Partial effect for distance to -200 isobath") 

### ---

#combine subplots into one plot and save to file

gam.land<-(gam.plot.shore.seasonal | gam.plot.shore.hour | gam.paraplot.shore.seasons) /
  (gam.plot.depth.seasonal | gam.plot.depth.hour | gam.paraplot.depth.seasons) /
  (gam.plot.200.seasonal | gam.plot.200.hour | gam.paraplot.200.seasons)

#ggsave(gam.land, filename='Plots/Maui_LandObservations_GAM.png', width=18.5, height=10.5, units='in',  bg = "white")    


################################################################

#### ----

#End of File

# -------

################################################################

