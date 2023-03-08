#################################################################
#################################################################

### DIEL PATTERNS OF INDIVIDUAL HUMPBACK WHALE SINGERS OFF MAUI, HAWAIʻI

# Written by Anke Kügler
# Last updated by Anke Kügler, 06/10/2022
# Contact: anke.kuegler@gmail.com

# This script analyses tracks of individual singers localized with vector-sensors (Directional Autonomous Seafloor Acoustic Recorders, DASARs)
# off Maui, Hawaiʻi in April 2020 to investigate diel patterns of the number of singers, spacing distances among singers, and movement behavior
# while singing.

# Input data: csv files containing coordinates of individual singers (tracks) every 60 s over 3 hours.

#             Note: Coordinates are in UTM space relative to the position of DASAR Y
#                   Tracks were manually compared across 24 hours and TrackIDs adjusted if multiple tracks belonged to the same individual

#Generalized Additive Models (GAMs) are fit for hourly averages of numbers of singers and minimum spacing distances and
#   day of the season
#   hour of the day

#A binomial GAM is fit for the probability of a singer to spend the majority (>70%) per hour remaining stationary, traveling, or both and
#   day of the season
#   hour of the day

#Functions:
# pkgTest() is a function that tests if a R package already exists on a system and if not, installs and loads it
# gam_plots_params() defines settings for the GAM effect plots
# fix_timestamp<-function() is a function that fixes erroneous timestamps
# get_tracks() is a function that reads in all 3-hour track files per data, combines all tracks, removes tracks shorter than 5 minutes
#   and localization outliers, smoothes the track with a running average, and for each localization of each track, calculates distances 
#   to the DASAR recorders, calculates the swimming speed and assigns a stationary/travel state based on the swimming speed.

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
    if(!require(x,character.only = TRUE)) stop("Package not found")
    else {
      require(x)
    }  
  }
}

# ---

pkgTest('dplyr') #to nicely summarize and mutate data in dataframes by groups; tidyverse alternative to aggregate
  options(dplyr.summarise.inform = FALSE) # Suppress dplyr::summarise info because they're heckin' annoying

pkgTest('anytime') #to fix timestamps with missing midnight hour:min:sec
pkgTest('lubridate') #date and time manipulations
pkgTest('chron') #date and time manipulations

pkgTest('stringr')

pkgTest('zoo') #for smoothing of tracks with rollapply

pkgTest('terra') #to project DASAR lat/long coordinates to UTM

pkgTest('mgcv') #to fit and analyse GAMs

pkgTest('ggplot2') #nice plots
pkgTest('patchwork') #to arrange multiple ggplots nicely into one

################################################################
#Functions

#GAM plots settings

gam_plots_params<-function(seasonal=F, diel=F,site=NULL, signif=T){
  list(
    if(signif)geom_ribbon(aes(ymin=fit-se, ymax=fit+se, group=1), 
                          fill = "grey75", alpha = 0.6),
    if(signif)geom_line(color='black'),
    if(!signif)geom_ribbon(aes(ymin=fit-se, ymax=fit+se, group=1), 
                           fill = "lightcoral", alpha = 0.3),
    if(!signif)geom_line(color='#ba6a68'),
    theme_bw(),
    if(signif)theme(axis.text = element_text(size = 12, color = "black"),
                    axis.title = element_text(size = 12, color = "black"),
                    axis.text.y = element_text(size = 12, color = "black", angle = 90, hjust=0.5),
                    panel.border=element_rect(color="black", size=0.75), 
                    panel.grid = element_blank(),
                    axis.line.x.top = element_line(color="black", size=0.75)),
    if(!signif)theme(axis.text = element_text(size = 12, color = "#ba6a68"),
                     axis.title = element_text(size = 12, color = "#ba6a68"),
                     axis.text.y = element_text(size = 12, color = "#ba6a68", angle = 90, hjust=0.5),
                     panel.border=element_rect(color="#ba6a68", size=0.75), 
                     panel.grid = element_blank(),
                     axis.line.x.top = element_line(color="#ba6a68", size=0.75),
                     axis.ticks= element_line(color="#ba6a68")),
    if(diel)scale_x_continuous(breaks=seq(from=0, to=23, by=3)/24, labels=c('0', '3', '6', '9', '12', '15', '18', '21')),
    if(diel)labs(x="Hour", title=site),
    if(seasonal)labs(x="Julian day", title="")
  )
}

fix_timestamp<-function(df){
  
  df.dates<-rep(NA, nrow(df))
  df.times<-rep(NA, nrow(df))
  
  timestamp_str<-str_split(df$timestamp, ' ')
  
  for (i in 1:nrow(df)){
    df.dates[i]<-timestamp_str[[i]][1]
    df.times[i]<-timestamp_str[[i]][2]
  }
  
  df.dates<-str_pad(df.dates, 8, pad = "0")
  df.dates<-paste(df.dates, '20', sep='')
  
  df.times<-str_pad(df.times, 5, pad = "0")
  df.times<-paste(df.times, ':00', sep='')
  
  df$timestamp<-as.POSIXct(paste(df.dates, df.times),  format='%m/%d/%Y %H:%M:%S')
  
  return(df)
}

### --- ---

get_tracks<-function(tracks_date, dat_dir, DASARs, whale_speed=NULL, new_trackID=F, ...){
  
# 1) get the tracks data from file(s) ---
  
  #initialize dataframe
  df<-data.frame()

  #find the folder from the provided general path and date name
  #structure needs to be: .../path to folder with result subfolders/yyyymmdd/individual_tracks.csv
  dat_dir<-paste(dat_dir, as.character(tracks_date), '/', sep='')
  
  #get all individual tracks files for that day
  files<-list.files(dat_dir, pattern='.csv')
 
  #have the tracks been cleaned up? They now have an "edited" suffix
  if(new_trackID==T){  
    files<-subset(files, grepl('_edited', files, fixed=T))
  }else{
    files<-subset(files, !grepl('_edited', files, fixed=T))
  }
  
  #read in each file and save to one dataframe
  for(i in 1:length(files)){
    dat<-read.csv(paste(dat_dir, files[i], sep=''))
    df<-rbind(df, dat)
  }
  
 # ---
  
  
# 2) fix timestamp if necessary ---
  
  #(when editing csv's, excel sometimes messes up the timestamp and saves in a different format
  #that's not understood by R)
  
  if(is.na(as.POSIXct(df$timestamp[1], format="%d-%b-%Y %H:%M:%S"))){
  df.dates<-rep(NA, nrow(df))
  df.times<-rep(NA, nrow(df))
  
  timestamp_str<-str_split(df$timestamp, ' ')
  
  for (i in 1:nrow(df)){
    df.dates[i]<-timestamp_str[[i]][1]
    df.times[i]<-timestamp_str[[i]][2]
  }
  
#  df.dates<-str_pad(df.dates, 8, pad = "0")
  df.dates<-paste(df.dates, '20', sep='')
  
  df.times<-str_pad(df.times, 5, pad = "0")
  df.times<-paste(df.times, ':00', sep='')
  
  df$timestamp<-as.POSIXct(paste(df.dates, df.times),  format='%m/%d/%Y %H:%M:%S', tz='UTC')
  }else{
    df$timestamp<-as.POSIXct(df$timestamp, format="%d-%b-%Y %H:%M:%S", tz='UTC')
  }

# ---
  
  
# 3) order by time within each track ---
  
  if(new_trackID){
    df<-df[order(df$TrackID_new, df$timestamp),]
  }else{
    df<-df[order(df$TrackID, df$timestamp),]    
  }
  
  row.names(df)<-NULL
  
  
# 4) remove rows with NA localization
  
  df<- df[!is.na(df$x),]
  
  # ---
  
  
# 5) remove tracks < 5 min
  
  #calculate length in minutes per track (based on TrackID)
  if(new_trackID){
    track.times<-as.data.frame(df%>%
                            group_by(TrackID_new)%>%
                            mutate(duration=n())%>%
                            summarise(duration=duration))
  }else{
    track.times<-as.data.frame(df%>%
                                 group_by(TrackID)%>%
                                 mutate(duration=n())%>%
                                 summarise(duration=duration))
  }
  
  #take out rows of tracks with duration less than 5 minutes
  
  df<-df[track.times$duration>5,]

# ---    
  
  
# 6) remove duplicates (based on TrackID and timestamp)
  
  #(this can be an artifact of subsequently merging split TrackIDs and briefly overlapping localizations)
  
  if(new_trackID){
    df<-df[!duplicated(df[c("timestamp","TrackID_new")]),]
  }else{
    df<-df[!duplicated(df[c("timestamp","TrackID")]),]
  }

# ---
  
  
# 7) convert localization coordinates to UTM easting/northing (relative to DASAR X, in km) ---
  
  df$easting<-df$x+DASARs[1,'easting']
  df$northing<-df$y + DASARs[1,'northing']
  
# ---
  
  
# 8) take out localization outliers
  
  #a localization is considered an outlier if its closest distance to the track is outside
  #the 95% quantile of all neighbor-neighbor distances of the track
  
  #calculate distance to nearest neighbor within each track
  groups<-split(df, df$TrackID_new)
  
  answer <-lapply(groups, function(df){
    #calculate distance matrix
    d<-dist(df[ , c("easting", "northing")], upper=TRUE, diag = TRUE)
    #convert to matrix
    dm<-as.matrix(d)
    
    #find the smallest value in each row less self
    closest <-sapply(1:nrow(dm), function(i){
      min(dm[i,-i])})
    
    #are you an outlier or not?
    df$outlier <-as.integer(closest>quantile(closest,probs=0.95))
    df
  })
  
  df<-dplyr::bind_rows(answer)
  
  #remove outliers from dataframe
  df<-df[df$outlier==0,]
  df<-df[,names(df)!='outlier'] #we don't need this column anymore
  
# ---
  
  
# 9) convert TrackIDs to factors ---
  
  df$TrackID<-as.factor(df$TrackID)
  if(new_trackID){df$TrackID_new<-as.factor(df$TrackID_new)}

# ---  
    
# 10) smooth tracks ---
  
  #zoo::rollyapply uses a running average of n for each x and y separately, smoothes track as function of time
  
  if(new_trackID){  
    df.fitted.sma<-as.data.frame(df %>%
                                   group_by(TrackID_new) %>%
                                   summarize(
                                      timestamp=timestamp,
                                      easting.SMA=rollapply(easting, 5, FUN=function(x) mean(x, na.rm=TRUE), by=1 , by.column=TRUE, partial=TRUE, align='right', fill=NA),
                                      northing.SMA=rollapply(northing, 5, FUN=function(x) mean(x, na.rm=TRUE), by=1 , by.column=TRUE, partial=TRUE, align='right', fill=NA)
                                   ))
  }else{
    df.fitted.sma<-as.data.frame(df %>%
                                   group_by(TrackID) %>%
                                   summarize(
                                      timestamp=timestamp,
                                      easting.SMA=rollapply(easting, 5, FUN=function(x) mean(x, na.rm=TRUE), by=1 , by.column=TRUE, partial=TRUE, align='right', fill=NA),
                                      northing.SMA=rollapply(northing, 5, FUN=function(x) mean(x, na.rm=TRUE), by=1 , by.column=TRUE, partial=TRUE, align='right', fill=NA)
                                   ))
  }
  
#  names(df.fitted.sma)[which(names(df.fitted) %in% c('easting', 'northing'))]<-c('easting.SMA', 'northing.SMA')

  #merge dfs back together
  if(new_trackID){    
    df<-merge(df,df.fitted.sma, by=c('TrackID_new','timestamp'), all.x=T)    
  }else{
    df<-merge(df,df.fitted.sma, by=c('TrackID','timestamp'), all.x=T)    
  }
  
# ---  
  
  
# 11) calculate distance to DASARs (in km)
  
  #note: euclidean distance using cartesian coordinates
  
  df$d_x<-sqrt((df$easting-DASARs$easting[DASARs$id=='X'])^2+(df$northing-DASARs$northing[DASARs$id=='X'])^2)
  df$d_y<-sqrt((df$easting-DASARs$easting[DASARs$id=='Y'])^2+(df$northing-DASARs$northing[DASARs$id=='Y'])^2)
  df$d_z<-sqrt((df$easting-DASARs$easting[DASARs$id=='Z'])^2+(df$northing-DASARs$northing[DASARs$id=='Z'])^2)

  df$d_x.sma<-sqrt((df$easting.SMA-DASARs$easting[DASARs$id=='X'])^2+(df$northing.SMA-DASARs$northing[DASARs$id=='X'])^2)
  df$d_y.sma<-sqrt((df$easting.SMA-DASARs$easting[DASARs$id=='Y'])^2+(df$northing.SMA-DASARs$northing[DASARs$id=='Y'])^2)
  df$d_z.sma<-sqrt((df$easting.SMA-DASARs$easting[DASARs$id=='Z'])^2+(df$northing.SMA-DASARs$northing[DASARs$id=='Z'])^2)
  
# ---
  
  
# 13) Get distance and swim speed between two subsequent localizations (in km and m/s)
  
  #use raw distance and distance between sma smoothed coordinates
  #calculate speed between subsequent localizations and distance averaged over 5 timesteps

  #function to calculate euclidean distance in cartesian coordinates
  get_d<-function(x, x.lag,y, y.lag){
    sqrt((x-x.lag)^2+(y-y.lag)^2)
  }
  
  #calculate subsequent distances
  #calculate time difference between subsequent localizations (in seconds)
  #calculate speed
  if(new_trackID){
    df<-as.data.frame(df%>%
                        group_by(TrackID_new) %>%
                        mutate(d=get_d(easting, lag(easting), northing, lag(northing)),
                               d.sma=get_d(easting.SMA, lag(easting.SMA), northing.SMA, lag(northing.SMA)),
                               delta_t=as.numeric(difftime(timestamp, dplyr::lag(timestamp, order_by=timestamp), units ="secs")),
                               v=d*1000/delta_t,
                               v.sma=d.sma*1000/delta_t,
                        ))  
  }else{
    df<-as.data.frame(df%>%
                        group_by(TrackID) %>%
                        mutate(d=get_d(easting, lag(easting), northing, lag(northing)),
                               d.sma=get_d(easting.sma, lag(easting.sma), northing.fitted, lag(northing.fitted)),
                               delta_t=as.numeric(difftime(timestamp, dplyr::lag(timestamp, order_by=timestamp), units ="secs")),
                               v=d*1000/delta_t,
                               v.sma=d.sma*1000/delta_t,
                               ))
  }
    

# ---  
  
  
# 14) assign stationary/travel state
  
  if(is.null(whale_speed)){
  whale_speed<-2*1000/3600} #2km/h, see Tenorio et al. 2022, Frankel et al. 1995, Noad & Cato 2007
  
  df$state<-apply(df['v'],1, function(x) ifelse(x>whale_speed,'travel','stationary')) 
  df$state.sma<-apply(df['v.sma'],1, function(x) ifelse(x>whale_speed,'travel','stationary')) 
  
# ---
  
  
#return df to save to object for subsequent analysis  
return(df)  
   
}


################################################################


################################################################
#######################  Analysis  #############################
################################################################

# 1) INITIALIZE SOME DATA ----


#DASAR locations (Lat/Long)
DASARs<-data.frame(id=c('X', "Y", "Z"),
                   Long=c(-156.6274,-156.6457,-156.6676),
                   Lat=c(20.8081,20.8292,20.8512))

#convert DASAR coordinates to UTM northing/easting
DASARs<-cbind(DASARs,as.data.frame(terra::project(as.matrix(DASARs[,c('Long', 'Lat')]), "+proj=longlat", "+proj=utm +zone=4 +units=m"))/1000)

colnames(DASARs)<-c("id", "Long", "Lat", "easting", "northing") #name columns


################################################################


# 2) GET DASAR TRACKS ----

#where do all my track localizations live?
dat_dir<-'SingerTracks/'

#process data per day
tracks.20200404<-get_tracks(20200404, dat_dir, DASARs, new_trackID=T)
tracks.20200406<-get_tracks(20200406, dat_dir, DASARs, new_trackID=T)
tracks.20200410<-get_tracks(20200410, dat_dir, DASARs, new_trackID=T)
tracks.20200414<-get_tracks(20200414, dat_dir, DASARs, new_trackID=T)
tracks.20200416<-get_tracks(20200416, dat_dir, DASARs, new_trackID=T)
tracks.20200418<-get_tracks(20200418, dat_dir, DASARs, new_trackID=T)
tracks.20200423<-get_tracks(20200423, dat_dir, DASARs, new_trackID=T)
tracks.20200429<-get_tracks(20200429, dat_dir, DASARs, new_trackID=T)

#combine all days into one dataframe
tracks<-rbind(tracks.20200404, 
              tracks.20200406,
              tracks.20200410,
              tracks.20200414,
              tracks.20200416,
              tracks.20200418, 
              tracks.20200423, 
              tracks.20200429)

### ---

# some house-keeping

#round timestamp to hour
tracks$timestamp.h<-floor_date(tracks$timestamp, unit="hour")

#subset to only whales within 6km of DASAR Y (middle of array)
tracks.6km<-subset(tracks, d_y.sma<=6)

#subset to only whales within every half hour
tracks.6km.30min<-subset(tracks.6km, (minute(timestamp)>=00 & minute(timestamp)<=30))





# 3) SUMMARIZE EFFORT ----

tracks.effort<-as.data.frame(tracks %>%
                               group_by(Date=as.Date(timestamp, tz='UTC'))%>%
                               summarise(n=length(unique(TrackID_new))
                               )           
)

tracks.6km.effort<-as.data.frame(tracks.6km %>%
                                   group_by(Date=as.Date(timestamp, tz='UTC'))%>%
                                   summarise(n.6km=length(unique(TrackID_new))
                                   )           
)

tracks.effort<-merge(tracks.effort, tracks.6km.effort)

#write.csv(tracks.effort, 'Results/Maui_SingerTracks_effort.csv', row.names=F)


################################################################


# 4) CALCULATE DATA OF INTEREST PER HOUR ----

# a - get number of whales within 6km every hour
# b - get distances to nearest neighbor
# c - get number of stationary/traveling/both states per hour



### ---



# (a) - get number of whales within 6km every hour

tracks.n<-as.data.frame(tracks.6km.30min %>% 
                           group_by(timestamp.h) %>%
                           summarise(whales.n=length(unique(TrackID_new))
                           ))


#get date and hour from timestamp
tracks.n$Hour<-floor_date(tracks.n$timestamp.h, unit="hour")
tracks.n$Hour<-times(strftime(tracks.n$Hour, "%H:%M:%S", tz='UTC'))
tracks.n$Date<-as.Date(tracks.n$timestamp.h, tz='UTC')

#convert date to julian day
tracks.n$julian = as.numeric(format(as.Date(tracks.n$timestamp.h, tz='UTC'), "%j"))


### ---

# (b) - get distances to nearest neighbor

#calculate mean location for each track and hour
tracks.dist<-as.data.frame(tracks.6km.30min %>% 
                          group_by(timestamp.h, TrackID_new) %>%
                          summarise(easting=mean(easting, na.rm=T),
                                    northing=mean(northing,na.rm=T),
                                    easting.sma=mean(easting.SMA, na.rm=T),
                                    northing.sma=mean(northing.SMA,na.rm=T)
                          ))


### ---

groups<-split(tracks.dist, tracks.dist$timestamp.h)

answer <-lapply(groups, function(tracks.dist){
  #calculate distance matrix
  d<-dist(tracks.dist[ , c("easting", "northing")], upper=TRUE, diag = TRUE)
  d.sma<-dist(tracks.dist[ , c("easting.sma", "northing.sma")], upper=TRUE, diag = TRUE)

  #convert to matrix
  dm<-as.matrix(d)
  dm.sma<-as.matrix(d.sma)

  #find the smallest value in each row less self
  closest <-sapply(1:nrow(dm), function(i){min(dm[i,-i])})
  closest.sma <-sapply(1:nrow(dm.sma), function(i){min(dm.sma[i,-i])})

  tracks.dist$min_d <-closest
  tracks.dist$min_d.sma <-closest.sma
  tracks.dist
})

tracks.dist<-dplyr::bind_rows(answer)

#clean up Inf for hours with only one track
tracks.dist$min_d[tracks.dist$min_d==Inf]<-NA
tracks.dist$min_d.sma[tracks.dist$min_d.sma==Inf]<-NA

#get date and hour from timestamp
tracks.dist$Hour<-times(strftime(tracks.dist$timestamp.h, "%H:%M:%S", tz='UTC'))
tracks.dist$Date<-as.Date(tracks.dist$timestamp.h,tz='UTC')

#convert date to julian day
tracks.dist$julian = as.numeric(format(as.Date(tracks.dist$timestamp.h, tz='UTC'), "%j"))


### ---


# (c) - get proportion of stationary/traveling/both states per hour

tracks.state<-as.data.frame(tracks.6km.30min %>% 
                          group_by(timestamp.h, TrackID_new) %>%
                          summarise(#how many minutes was whale present per hour
                                    samples=n(),
                                    #how many minutes were spent stationary
                                    stationary.t=sum(state.sma=='stationary', na.rm=T),
                                    #how many minutes were spent traveling
                                    travel.t=sum(state.sma=='travel', na.rm=T)
                          )%>%
                          mutate(#was whale predominantly stationary, traveling, or both during a given hour
                                stationary_percent=stationary.t/samples, 
                                travel_percent=travel.t/samples,
                                stationary_binary=as.integer(stationary_percent>0.7),
                                travel_binary=as.integer(travel_percent>0.7),
                                both_binary=as.integer((0.3<travel_percent & travel_percent<0.7) | (0.3<stationary_percent & stationary_percent<0.7)),
                          ))


#convert dataframe to long format
tracks.state.binary<-reshape2::melt(tracks.state[,c('timestamp.h', 'travel_binary', 'stationary_binary', 'both_binary')], id='timestamp.h')
names(tracks.state.binary)<-c('timestamp.h','behavior','state')

#get date and hour from timestamp
tracks.state.binary$Hour<-times(strftime(tracks.state.binary$timestamp.h, "%H:%M:%S", tz='UTC'))
tracks.state.binary$Date<-as.Date(tracks.state.binary$timestamp.h,tz='UTC')

#convert date to julian day
tracks.state.binary$julian<-as.numeric(format(as.Date(tracks.state.binary$timestamp.h, tz='UTC'), "%j"))


### ---

# - pool by hour -

tracks.state.sum<-as.data.frame(tracks.state %>%
                                  group_by(timestamp.h) %>%
                                  summarise(#how many whales per hour
                                            n=n(), 
                                            #how many were stationary
                                            stationary.n=sum(stationary_binary),
                                            
                                            #how many were traveling
                                            travel.n=sum(travel_binary),

                                            #how many were both
                                            both.n=sum(both_binary),

                                            #calculate percent of animals in different states out of all animals per hour
                                            stationary.p=stationary.n/n,
                                            travel.p=travel.n/n,
                                            both.p=both.n/n,
                               ))

###


################################################################

# 4) SUMMARIZE DATA (descriptive Statistics) ----


tracks.n.summary<-as.data.frame(tracks.n %>%
                                  group_by(Hour)%>%
                                  summarise(whales.n.median=median(whales.n),
                                            whales.n.IQR=IQR(whales.n)
                                  ))       


### ---

tracks.dist.summary<-as.data.frame(tracks.dist %>%
                                     group_by(Hour=times(strftime(timestamp.h, "%H:%M:%S", tz='UTC')))%>%
                                     summarise(min_d.median=median(min_d, na.rm=T),
                                               min_d.sma.IQR=IQR(min_d.sma, na.rm=T)
                                      ))       

### ---

tracks.state.summary<-as.data.frame(tracks.state.sum %>%
                                      group_by(Hour=times(strftime(timestamp.h, "%H:%M:%S", tz='UTC')))%>%
                                      summarise(travel.median=median(travel.p, na.rm=T),
                                                travel.IQR=IQR(travel.p, na.rm=T),
                                                stationary.median=median(stationary.p, na.rm=T),
                                                stationary.IQR=IQR(stationary.p, na.rm=T),
                                                both.median=median(both.p, na.rm=T),
                                                both.IQR=IQR(both.p, na.rm=T)
                                      ))       

### ---

tracks.summary<-merge(tracks.n.summary, tracks.dist.summary)
tracks.summary<-merge(tracks.summary, tracks.state.summary)

#write.csv(tracks.summary, 'Results/Maui_SingerTracks_summary.csv', row.names=F)


################################################################

# 5) STATISTICAL ANALYSIS (Generalized Additive Models) ----


gam.tracks.n<-gam(whales.n ~ s(julian, k=8) + s(Hour, bs='cc'), data=tracks.n, method='REML', select=T)

anova.gam.n<-anova(gam.tracks.n)

### ---

gam.tracks.dist<-gam(min_d.sma ~ s(julian, k=7) + s(Hour, bs='cc'), data=tracks.dist, method='REML', select=T)
anova.gam.dist<-anova(gam.tracks.dist)

### ---

gam.tracks.state<-gam(state~ s(julian, k=8, by=behavior) + s(Hour, bs='cc', by=behavior), data=tracks.state.binary, method='REML', family='binomial', select=T)
anova.gam.state<-anova(gam.tracks.state)


#---

#Create effect plots. Because plot.gam() doesn't allow for a lot of adjustments of the look of the plots,
#we're going to create them from scratch.

# (a) - pull out all smoothers 

smoothers.tracks.n <- plot(gam.tracks.n, pages=1, scale=0)  # plot.gam returns a list of n elements, one per plot
smoothers.tracks.dist <- plot(gam.tracks.dist, pages=1, scale=0)  
smoothers.tracks.state <- plot(gam.tracks.state, pages=1, scale=0) 

#get x, fit, and se for each plot and save to DF
get_GAMsmoothers<-function(smoothers,n){
  sm.coef<-smoothers[[n]]
  sm.fit<-as.data.frame(sm.coef$fit)
  names(sm.fit)<-'fit'
  sm.fit$x<-as.vector(sm.coef$x)
  sm.fit$se<-as.vector(sm.coef$se)
  return(sm.fit)
}

# Number of singers

sm.julian.n<-get_GAMsmoothers(smoothers.tracks.n,1)
sm.diel.n<-get_GAMsmoothers(smoothers.tracks.n,2)

# Minimum distance

sm.julian.dist<-get_GAMsmoothers(smoothers.tracks.dist,1)
sm.diel.dist<-get_GAMsmoothers(smoothers.tracks.dist,2)

# Stationary/Traveling/Both

sm.julian.stationay<-get_GAMsmoothers(smoothers.tracks.state,2)
sm.diel.stationary<-get_GAMsmoothers(smoothers.tracks.state,5)

sm.julian.travel<-get_GAMsmoothers(smoothers.tracks.state,1)
sm.diel.travel<-get_GAMsmoothers(smoothers.tracks.state,4)

sm.julian.both<-get_GAMsmoothers(smoothers.tracks.state,3)
sm.diel.both<-get_GAMsmoothers(smoothers.tracks.state,6)


###---

# (b) - create individual plots and save as objects

# Number of singers

gam.plot.n.julian<-ggplot(data=sm.julian.n, aes(x=x, y=fit)) +
  gam_plots_params(seasonal=T, signif=ifelse(anova.gam.n$s.table[1,4]<=0.05, TRUE, FALSE)) +
  labs(y="Partial effect for number of whales") +
  scale_y_continuous(labels=prettyNum)

gam.plot.n.diel<-ggplot(data=sm.diel.n, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, signif=ifelse(anova.gam.n$s.table[2,4]<=0.05, TRUE, FALSE)) +
  labs(y="Partial effect for number of whales") +
  scale_y_continuous(labels=prettyNum)

# Minimum distance

gam.plot.dist.julian<-ggplot(data=sm.julian.dist, aes(x=x, y=fit)) +
  gam_plots_params(seasonal=T, signif=ifelse(anova.gam.dist$s.table[1,4]<=0.05, TRUE, FALSE)) +
  labs(y="Partial effect for min neighbor dist.") +
  scale_y_continuous(labels=prettyNum)

gam.plot.dist.diel<-ggplot(data=sm.diel.dist, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, signif=ifelse(anova.gam.dist$s.table[2,4]<=0.05, TRUE, FALSE)) +
  labs(y="Partial effect for min neighbor dist.") +
  scale_y_continuous(labels=prettyNum)

# Stationary 

gam.plot.stationary.julian<-ggplot(data=sm.julian.stationay, aes(x=x, y=fit)) +
  gam_plots_params(seasonal=T, signif=ifelse(anova.gam.state$s.table[2,4]<=0.05, TRUE, FALSE)) +
  labs(y='Partial effect', title='Stationary') +
  scale_y_continuous(labels=prettyNum)

gam.plot.stationary.diel<-ggplot(data=sm.diel.stationary, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, signif=ifelse(anova.gam.state$s.table[5,4]<=0.05, TRUE, FALSE)) +
  labs(y="Partial effect", title='Stationary') +
  scale_y_continuous(labels=prettyNum)

# Traveling

gam.plot.travel.julian<-ggplot(data=sm.julian.travel, aes(x=x, y=fit)) +
  gam_plots_params(seasonal=T, signif=ifelse(anova.gam.state$s.table[1,4]<=0.05, TRUE, FALSE)) +
  labs(y='Partial effect', title='Traveling') +
  scale_y_continuous(labels=prettyNum)  

gam.plot.travel.diel<-ggplot(data=sm.diel.travel, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, signif=ifelse(anova.gam.state$s.table[4,4]<=0.05, TRUE, FALSE)) +
  labs(y="Partial effect", title='Traveling') +
  scale_y_continuous(breaks=c(-0.001, 0, 0.001), labels=prettyNum)

# Both Stationary and Traveling

gam.plot.both.julian<-ggplot(data=sm.julian.both, aes(x=x, y=fit)) +
  gam_plots_params(seasonal=T, signif=ifelse(anova.gam.state$s.table[3,4]<=0.05, TRUE, FALSE)) +
  labs(y='Partial effect', title='Stationary+Traveling') +
  scale_y_continuous(labels=prettyNum)


gam.plot.both.diel<-ggplot(data=sm.diel.both, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, signif=ifelse(anova.gam.state$s.table[6,4]<=0.05, TRUE, FALSE)) +
  labs(y="Partial effect", title='Stationary+Traveling') +
  scale_y_continuous(labels=prettyNum)

### ---

#combine subplots into one plot and save to file

gam.plots.tracks<-(gam.plot.n.julian | gam.plot.n.diel) /
                  (gam.plot.dist.julian | gam.plot.dist.diel) / 
                  (gam.plot.stationary.julian | gam.plot.stationary.diel)/
                  (gam.plot.travel.julian | gam.plot.travel.diel)/
                  (gam.plot.both.julian | gam.plot.both.diel)

#ggsave(gam.plots.tracks, filename='Plots/Maui_SingerTracks_GAM.png', width=12, height=14, units='in',  bg = "white")    


################################################################

#### ----

#End of File

# -------

################################################################
