#################################################################
#################################################################

### DIEL PATTERNS OF HUMPBACK WHALE CHORUSING OFF MAUI, HAWAIʻI

# Written by Anke Kügler
# Last updated by Anke Kügler, 12/05/2023
# Contact: anke.kuegler@gmail.com
# Corresponds with analyses and figures in Kügler et al. (in Review)

# This script analyses hourly averages of acoustic recordings (RMS SPL) 
# collected from autonomous passive acoustic recorders (Ecological Acoustic Recorders, EARs)
# off Maui, Hawaiʻi between 2016 and 2021 to investigate diel patterns of male humpback whale (HW) chorusing.
# In addition, ambient noise in the adjacent frequency band to HW chorusing is analysed to investigate
# possible impacts of environmental noise. 

# Input data: text files containing RMS SPL values for each 30-second recording per deployment and site

# Levels are analysed in two frequency bands:
# 1: 0-1.5 kHz obtained from subsampled (3 kHz) recordings, associated with humpback whale chorusing
# 2: 1.56-3.23 kHz obtained as octave bands from the fullband recordings, adjacent to HW band

#Generalized Additive (Mixed) Models (GAM(M)s) are fit for hourly averages for both frequency bands and 
#   day of the season
#   hour of the day
#   season (representing respective whale seasons from December through April)
#   site

#Generalized Additive (Mixed) Models are fit for hourly averages for a subset of the data for both frequency bands and
#   day of the season
#   hour of the day
#   season (representing respective whale seasons from December through April)
#   site
#to illustrate differences in diel patters during the beginning, peak, and end of the season corresponding with differing
#whale densities

#Generalized Additive Models (GAMs) are fit for hourly averages in August/September 2016 at one site for both frequency bands and 
#   day of the season
#   hour of the day
#to illustrate seasonal and diel patterns during the off-season. 
#There was no additional off-season data available to include in these analyses.

#Functions:
# pkgTest() is a function that tests if a R package already exists on a system and if not, installs and loads it
# get_rms.EAR() reads in the raw RMS SPL data from file, creates date and hour timestamps, corrects dates for
#   leap years if necessary, and averages (mean or median) the data per hour per day.
# get_GAMsmoothers() returns the x, fit, and se used to create smoother effect plots in plot.gam()
# gam_plots_params(), gam_paraplots_params(), and gam_plots_sub.diel_params() define settings for the GAM effect plots
# polarplots_params() and octaves_plots_params() define settings for the visualization plots
# scale_fun() and normalize_fun() are functions to standardize the data


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

pkgTest('stringr')

pkgTest('chron') #date and time manipulations
pkgTest('lubridate') #date and time manipulations

pkgTest('dplyr') #to nicely summarize and mutate data in dataframes including by groups; tidyverse alternative to aggregate
  options(dplyr.summarise.inform = FALSE) # Suppress dplyr::summarise info because they're heckin' annoying
pkgTest('reshape') #data frame manipulation

pkgTest('mgcv') #to fit and analyse GAMs

pkgTest('suncalc') #to get sunrise/sunset times based on Lat/Long and date

#necessary packages for plotting
pkgTest('ggplot2')
pkgTest('scales') #nicer ggplot axes
pkgTest('geomtextpath')
pkgTest('ungeviz') #draw point as horizontal line in GAM parameter plot
pkgTest('ggpubr') #ggarrange (combine multiple ggplots into one plot)
pkgTest('patchwork') #combine multiple ggplots into one plot


################################################################
#Functions

#read in EAR RMS SPL files and compute hourly medians per day
get_rms.EAR<-function(file=NULL, fs=NULL, year=NULL, site=NULL, round_hour=F, sub=F,...){
  
  #get optionally provided arguments
  opt_args <- list(...)
  
  ###
  
  #GET DATA AND PARAMETERS:
  
  if(is.null(file)==T | class(file)!='character'){ 
    cat('Please choose the file with the RMS SPL output\n')
    file<-file.choose()
  }
  
  rms<-read.table(file, header=F)
  
  
  ###  
  
  #COLUMN NAMES:
  
  #get all octave frequency bands if they were calculated and name columns
  if(length(names(rms))>4){
    band5<-fs/2
    band4<-band5/2
    band3<-band4/2
    band2<-band3/2
    band1<-band2/2
    bands<-c(paste(prettyNum(band4, digits=3, drop0trailing=T),"-", prettyNum(band5, digits=3, drop0trailing=T), ' kHz', sep=""), 
             paste(prettyNum(band3, digits=3, drop0trailing=T),"-", prettyNum(band4, digits=3, drop0trailing=T), ' kHz', sep=""), 
             paste(prettyNum(band2, digits=3, drop0trailing=T),"-", prettyNum(band3, digits=3, drop0trailing=T), ' kHz', sep=""), 
             paste(prettyNum(band1, digits=3, drop0trailing=T),"-", prettyNum(band2, digits=3, drop0trailing=T), ' kHz', sep=""), 
             paste("0","-", prettyNum(band1, digits=3, drop0trailing=T), ' kHz', sep=""))
    
    names(rms) <- c("File", "Date", "Time", "Fullband", rev(bands))
    all.bands<-c('Fullband',rev(bands))
  } else{
    bands<-NULL
    names(rms)<-c("File", "Date", "Time", paste("0", "-", prettyNum(fs/2, digits=3, drop0trailing=T), ' kHz',sep=""))
    all.bands<-paste("0", "-", prettyNum(fs/2, digits=3, drop0trailing=T), ' kHz',sep="")
  }
  
  #Take out NA columns. This safeguards if octave bands couldn't be calculated (contain NaNs)
  rms<-rms[ , apply(rms, 2, function(x) !any(is.na(x)))]
  
  ###
  
  #CHANGE DATE AND TIME FORMAT:
  
  #add year: test if a new year starts during the recording period and adjust accordingly
  if("101" %in% rms[,'Date']==T){
    
    i<-min(which(rms[,'Date']=="101"))
    rms[1:(i-1),'Date']<-paste(as.character(year), rms[1:(i-1),'Date'], sep="-")
    rms[i:length(rms[,'Date']),'Date']<-paste(as.character(year+1), rms[i:length(rms[,'Date']),'Date'], sep="-")
    
  } else{
    rms[,'Date']<-paste(as.character(year), rms[,'Date'], sep="-")
  } 
  
  #add deliminator between month and day
  str_sub(rms[,'Date'], nchar(rms[,'Date'])-1, nchar(rms[,'Date'])-2)<-"-"
  
  #convert to date
  rms$Date<-as.Date(rms$Date, tz="UTC")
  
  
  #---
  
  #adjust time format: 
  
  #add seconds
  rms$Time<-paste0(rms$Time,'00')
  #add leading zeros for minutes and hour
  rms$Time<-str_pad(rms$Time, 6, pad = "0") 
  #add deliminator between minutes and seconds
  str_sub(rms$Time, nchar(rms$Time)-1, nchar(rms$Time)-2)<-":"
  #add deliminator between hours and minutes
  str_sub(rms$Time,nchar(rms$Time)-4, nchar(rms$Time)-5)<-":"
  
  #convert timestamps to time/POSIX objects
  rms$Time<-times(strftime(strptime(rms$Time, "%H:%M:%S", tz='UTC'), "%H:%M:%S", tz="UTC")) 
  
  #---
  
  #Check and adjust for leap years if necessary
  
  for(n in 1:length(unique(year(rms$Date)))){
    
    check.year<-unique(year(rms$Date))[n]
    
    if(as.Date(paste(check.year,"03-01", sep='-'), tz="UTC") %in% rms[,'Date']==T && leap.year(check.year)==T){
      
      i.leap.1<-min(which(rms[,'Date']==paste(check.year,"03-01", sep='-')))
      i.leap.2<-min(which(rms[,'Date']==paste(check.year,"03-02", sep='-')))
      
      #insert YEAR/02/29 for date in all rows that are YEAR/03/01
      rms[i.leap.1:(i.leap.2-1),'Date']<-as.Date(paste(as.character(check.year), '2-29', sep="-"), tz="UTC")
      
      #decrease date of all subsequent days by 1 to shift dates
      rms[i.leap.2:nrow(rms),'Date']<-rms[i.leap.2:nrow(rms),'Date'] - 1
    }
  }
  
  #---  
  
  ###
  
  #SUMMARIZE DATA:
  
  if (round_hour==T){
    
    rms$Date<-round_date(rms$Date, unit="hour") #round to full hour
  }
  

  #pool by hour
  
  rms$Hour<-as.POSIXct(as.character(rms$Time), format = "%H:%M:%S", tz="UTC")
  
  if (sub==T){
    #subset to only whales within every half hour
    rms<-subset(rms, (minute(Hour)>=00 & minute(Hour)<=30))
  }
  
  rms$Hour<-round_date(rms$Hour, unit="hour")
  rms$Hour<-times(strftime(rms$Hour, "%H:%M:%S", tz="UTC"))
  
  rms.h<-aggregate(.~Hour+Date, data=rms, median, na.action=na.omit) 
  
  
  ####
  
  
  #RETURN FROM FUNCTION
  
  #Save to object/list for subsequent analysis
  
  return(rms.h[,-c(3,4)]) #take out "File"  and Time columns
  
} 

#---

#Function to get x, fit, and se for each smooth plot of a GAM and save to DF
get_GAMsmoothers<-function(smoothers,n){
  sm.coef<-smoothers[[n]]
  sm.fit<-as.data.frame(sm.coef$fit)
  names(sm.fit)<-'fit'
  sm.fit$x<-as.vector(sm.coef$x)
  sm.fit$se<-as.vector(sm.coef$se)
  return(sm.fit)
}

#---

#GAM plots settings

gam_plots_params<-function(seasonal=F, diel=F,site=NULL, hw=F, env=F, se=T, secondary=F){
  list(
    if(se)geom_ribbon(aes(ymin=fit-se, ymax=fit+se, group = 1), 
                fill = "grey75", alpha = 0.6),
    geom_line(color='black'),
    theme_bw(),
    theme(axis.text = element_text(size = 16, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black", angle = 90),
          plot.title=element_text(size = 16, color = "black")),
    theme(panel.border=element_rect(color="black", size=0.75), 
          panel.grid = element_blank()),
    theme(axis.line.x.top = element_line(color="black", size=0.75)),
    labs(y="Partial effect"),
    if(diel)scale_x_continuous(breaks=seq(from=0, to=23, by=3)/24, labels=c('0', '3', '6', '9', '12', '15', '18', '21')),
    if(diel) labs(x="Hour", title=site),
    if(seasonal)labs(x="Day of the season", title=""),
    if(env & !secondary) scale_y_continuous(limits=c(-3.5, 11), breaks=c(-2, 0, 2, 4, 6, 8, 10)),
    if(seasonal)scale_y_continuous(limits=c(-14, 7), breaks=c(-10, -5, 0, 5)),
    if(diel & hw)scale_y_continuous(limits=c(-2.5, 2.5), breaks=c(-2, -1, 0, 1,2)),
    if(secondary)theme(axis.line.y.right = element_line(color = 'red'), 
                       axis.ticks.y.right = element_line(color = 'red'),
                       axis.text.y.right = element_text(color = 'red'), 
                       axis.title.y.right = element_blank())
  )
}

#---

gam_paraplots_params<-function(site=F, season=F, hw=F, env=F){
  list(
    if(site)geom_errorbar(aes(ymin=site.y-site.se, ymax=site.y+site.se), 
                          width=0.75, position=position_dodge(width=0.9), stat = "identity", color='grey65'),
    if(season)geom_errorbar(aes(ymin=season.y-season.se, ymax=season.y+season.se), 
                            width=0.75, position=position_dodge(width=0.9), stat = "identity", color='grey65'),
    geom_hpline(stat='identity', size=1, width=0.75),
    theme_bw(),
    theme(axis.text = element_text(size = 16, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black", angle = 90),
          plot.title=element_text(size = 16, color = "black")),
    theme(panel.border=element_rect(color="black", size=0.75), 
          panel.grid = element_blank()),
    theme(axis.line.x.top = element_line(color="black", size=0.75)),
    labs(x='', y="Partial effect", title=''),
    #scale_x_discrete(labels = function(x) str_wrap(x, width = 10)),
    if(site)scale_x_discrete(labels=c(paste0("<span style='font-size: 14pt'>Olowalu</span><br><span style='font-size: 13pt'>(shallow)</span>"),
                     "<span style='font-size: 16pt'>Maui7</span>", 
                     paste0("<span style='font-size: 14pt'>Olowalu</span><br><span style='font-size: 13pt'>(deep)</span>"), 
                     "<span style='font-size: 16pt'>Maui6</span>",
                     "<span style='font-size: 16pt'>Pali</span>")),
    if(site)theme(axis.text.x=ggtext::element_markdown()),
#    if(env) scale_y_continuous(limits=c(-3.5, 11), breaks=c(-2, 0, 2, 4, 6, 8, 10)), #, labels=c('','0','','4','','8','')
    if(hw) scale_y_continuous(limits=c(-2, 9), breaks=c(-2,0,2,4,6,8)),
    if(env) scale_y_continuous(limits=c(-2, 9), breaks=c(-2,0,2,4,6,8))
  )
}

#---

gam_plots_sub.diel_params<-function(ribbon=T, int=F){
  list(
    geom_rect(data=sun.times, aes(xmin = 0, xmax = sunrise, ymin = -Inf, ymax = Inf),alpha = 0.42, fill = "grey90", inherit.aes=F),
    geom_rect(data=sun.times, aes(xmin = sunset, xmax = Inf, ymin = -Inf, ymax = Inf),alpha = 0.42, fill = "grey90", inherit.aes=F),
    if(ribbon)geom_ribbon(aes(ymin=fit-se, ymax=fit+se, fill=site, group=1), alpha = 0.21),
    if(!int)geom_line(aes(color=site), linewidth=1.5),
    if(int)  geom_line(aes(color=interaction(site,fband)), linewidth=1.5),
    theme_bw(),
    labs(x='Hour', y="Partial effect"),
    theme(axis.text = element_text(size = 16, color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black", angle = 90),
          panel.border=element_rect(color="black", size=0.75), 
          panel.grid = element_blank(),
          axis.line.x.top = element_line(color="black", size=0.75),
          legend.position='none',
          strip.background=element_blank(),
          strip.text=element_text(size=16)),
    scale_x_continuous(breaks=seq(from=0, to=23, by=3)/24, 
                       labels=c('0', '3', '6', '9', '12', '15', '18', '21')),
    scale_color_manual(values=c('Maui6'='#377b26','Olowalu (shallow)'='#93d275',
                                'Maui6.HW'='#377b26','Olowalu (shallow).HW'='#93d275',
                                'Maui6'=muted('#104566'), 'Olowalu (shallow)'='#729bc9',
                                'Maui6.Env'=muted('#104566'), 'Olowalu (shallow).Env'='#729bc9')),
    scale_fill_manual(values=c('Maui6'='#377b26','Olowalu (shallow)'='#93d275',
                               'Maui6.HW'='#377b26','Olowalu (shallow).HW'='#93d275',
                               'Maui6'=muted('#104566'), 'Olowalu (shallow)'='#729bc9',
                               'Maui6.Env'=muted('#104566'), 'Olowalu (shallow).Env'='#729bc9')),
    facet_grid(rows=vars(factor(site,levels=c('Olowalu (shallow)', 'Maui6'))), 
               cols=vars(factor(seasontime, levels=c('early','peak','late'))),
               scales='free')
  )
}

#---

#Polarplots and octave band comparison plots settings

polarplots_params<-function(hw=F, env=F){
  list(
    theme_minimal(),
      facet_wrap(~factor(site,levels=c('Olowalu (shallow)', 'Maui7', 'Olowalu (deep)', 'Pali', 'Maui6')), nrow=1),
      geom_hline(yintercept = seq(0, 1, by = 0.25), colour = "grey75", size = 0.2),
      geom_vline(xintercept = as.numeric(seq(as.POSIXct("00:00:00", format = "%H:%M:%S", tz='UTC'),as.POSIXct("23:00:00",format = "%H:%M:%S",  tz='UTC'),'1 hour')), colour = "grey75", size = 0.2),
      geom_bar(stat='identity'),
      if(hw)scale_fill_gradient2(position="bottom" , low = "#2D4E1D", mid = muted("#2D4E1D"), high = "#93d275", 
                                 midpoint = median(rms.hw.h$SPL.norm.median.scaled)),
      if(env)  scale_fill_gradient2(position="bottom" , low = "#104566", mid = muted('#104566'), high = "turquoise3",  ##105966
                                    midpoint = median(rms.env.h$SPL.norm.median.scaled)),
      scale_x_datetime(breaks=seq(as.POSIXct("00:00:00", format = "%H:%M:%S", tz='UTC'),as.POSIXct("23:00:00",format = "%H:%M:%S",  tz='UTC'),'1 hour'), 
                       labels = date_format("%H"), timezone='UTC',
                       expand=c(0,0)),
      labs(x="", y="", title=""),  
      theme(legend.position='none',
            axis.text.y = element_blank(),
            axis.text.x = element_text(vjust=0.25, size=14, color='black'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1, "lines"),
            strip.text = element_text (size=16, color='black', margin = margin (t=0, r=0, b=0.5, l=0, 'lines')),
            axis.line.x=element_blank()),
      scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1)),
      coord_curvedpolar(start=(-(15*pi)/(180)/2))    
  )
}

#---

octaves_plots_params<-function(seasonal=F, diel=F, site=NULL){
  list(
    if(seasonal)  geom_line(aes(y=ave(value, variable, FUN = function(Z) zoo::rollmean(Z, 10, fill=NA, align='center')), 
                                color=variable), na.rm= TRUE, size=1),
    if(diel)  geom_line(aes(color=variable), size=1),
    theme_classic(),
    theme(legend.position='bottom', legend.text=element_text(size=19), legend.title=element_blank()),
      theme(axis.text = element_text(size=19, color="black"),
            plot.title = element_text(size=19, color="black")),
      if(seasonal)theme(axis.title.y = element_text(size = 19, color="black")),
      if(diel)theme(axis.title.y = element_text(size = 19, color="black")),
      if(diel)theme(axis.title.x = element_text(size = 19, margin = margin(t = 10, r = 0, b = 0, l = 0))),
      if(seasonal)scale_x_date(date_breaks='1 month', date_labels = '%b', expand = expansion(mult=c(0.05, 0))),
      if(diel)scale_x_datetime(breaks=seq(as.POSIXct("00:00:00", format = "%H:%M:%S", tz='UTC'),as.POSIXct("23:00:00",format = "%H:%M:%S",  tz='UTC'),'3 hour'), labels = function(z) gsub("^0", "", strftime(z, "%H", tz='UTC')), timezone='UTC', expand = expansion(mult=c(0.05, 0))),
      scale_y_continuous(limits=c(89,117), breaks=c(90, 95, 100, 105, 110, 115)),
      guides(colour = guide_legend(nrow = 1)),
      scale_colour_manual(values=cols),
      if(seasonal)labs(x="", y=lab.y.day, title=site),
      if(diel)labs(x="Hour of Day", y=lab.y.h, title=site)
  )
}

#---

#scaling functions
scale_fun <- function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}
normalize_fun <- function(x){x/mean(x, na.rm=T)}


################################################################

################################################################
#######################  Analysis  #############################
################################################################

#### 1) GET EAR DATA ----

#where does my data live?
dat.dir<-'RMSSPL/'

# get data in the 0-1.5 Hz frequency band (humpbacks) ---

rms.maui6.2017<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2016-12-22to2017-05-04)RMS_SPL_3kHz.txt", sep=''),3, 2016)
rms.maui6.2018<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2017-12-07to2018-04-30)RMS_SPL_3Khz.txt", sep=''),3, 2017)
rms.maui6.2019<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2018-12-01to2019-04-24)RMS_SPL_3kHz.txt", sep=''),3, 2018)
rms.maui6.2020<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2019-11-15to2020-05-16)RMS_SPL_3kHz.txt", sep=''),3, 2019)
rms.maui6.2021<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2020-12-07to2021-05-26)RMS_SPL_3kHz.txt", sep=''),3, 2020)

###

rms.olowalu.deep.2017<-get_rms.EAR(paste(dat.dir, "Olowalu_deep_EAR(2016-12-22to2017-05-04)RMS_SPL_3kHz.txt", sep=''),3, 2016)
rms.olowalu.deep.2019<-get_rms.EAR(paste(dat.dir, "Olowalu_deep_EAR(2018-12-01to2019-04-25)RMS_SPL_3kHz.txt", sep=''),3, 2018)
rms.olowalu.deep.2020<-get_rms.EAR(paste(dat.dir, "Olowalu_deep_EAR(2019-11-06to2020-03-30)RMS_SPL_3kHz.txt", sep=''),3, 2019)
rms.olowalu.deep.2021<-get_rms.EAR(paste(dat.dir, "Olowalu_deep_EAR(2020-12-01to2021-05-26)RMS_SPL_3kHz.txt", sep=''),3, 2020)

###

rms.olowalu.shallow.2016<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2016-07-30to2017-01-08)RMS_SPL_3kHz.txt", sep=''),3, 2016)
  names(rms.olowalu.shallow.2016)[3]<-'0-1.5 kHz'
rms.olowalu.shallow.2017<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2017-01-09to2017-05-17)RMS_SPL_3kHz.txt", sep=''),3, 2017)
rms.olowalu.shallow.2018<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2017-12-07to2018-05-17)RMS_SPL_3kHz.txt", sep=''),3, 2017)
rms.olowalu.shallow.2019<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2018-11-02to2019-03-14)RMS_SPL_3kHz.txt", sep=''),3, 2018)
rms.olowalu.shallow.2020<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2019-11-22to2020-05-09)RMS_SPL_3kHz.txt", sep=''),3, 2019)
rms.olowalu.shallow.2021<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2020-12-11to2021-06-18)RMS_SPL_3kHz.txt", sep=''),3, 2020)

###

rms.pali.2017<-get_rms.EAR(paste(dat.dir, "Pali_EAR(2016-12-22to2017-05-04)RMS_SPL_3kHz.txt", sep=''),3, 2016)

rms.maui7.2017<-get_rms.EAR(paste(dat.dir, "Maui7_EAR(2017-01-04to2017-05-04)RMS_SPL_3kHz.txt", sep=''),3, 2017)


#---

# get data in the 1.5-3 Hz frequency band (non-humpback environmental) ---

#Note: raw RMS SPL results contain fullband and octave band levels, but we're only interested in 1.5-3 Hz
#      hence subsetting the dataframe after reading in and processing the data.

rms.maui6.2017.env<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2016-12-22to2017-05-04)RMS_SPL_50kHz.txt", sep=''),50, 2016)
rms.maui6.2018.env<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2017-12-07to2018-04-30)RMS_SPL_50Khz.txt", sep=''),50, 2017)
rms.maui6.2019.env<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2018-12-01to2019-04-24)RMS_SPL_50kHz.txt", sep=''),50, 2018)
rms.maui6.2020.env<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2019-11-15to2020-05-16)RMS_SPL_25kHz.txt", sep=''),25, 2019)
rms.maui6.2021.env<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2020-12-07to2021-05-26)RMS_SPL_25kHz.txt", sep=''),25, 2020)

rms.maui6.2017.env<-rms.maui6.2017.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.maui6.2018.env<-rms.maui6.2018.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.maui6.2019.env<-rms.maui6.2019.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.maui6.2020.env<-rms.maui6.2020.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.maui6.2021.env<-rms.maui6.2021.env[,c('Date', 'Hour', '1.56-3.12 kHz')]

###

rms.olowalu.deep.2017.env<-get_rms.EAR(paste(dat.dir, "Olowalu_deep_EAR(2016-12-22to2017-05-04)RMS_SPL_50kHz.txt", sep=''),50, 2016)
rms.olowalu.deep.2019.env<-get_rms.EAR(paste(dat.dir, "Olowalu_deep_EAR(2018-12-01to2019-04-25)RMS_SPL_50kHz.txt", sep=''),50, 2018)
rms.olowalu.deep.2020.env<-get_rms.EAR(paste(dat.dir, "Olowalu_deep_EAR(2019-11-06to2020-03-30)RMS_SPL_25kHz.txt", sep=''),25, 2019)
rms.olowalu.deep.2021.env<-get_rms.EAR(paste(dat.dir, "Olowalu_deep_EAR(2020-12-01to2021-05-26)RMS_SPL_25kHz.txt", sep=''),25, 2020)

rms.olowalu.deep.2017.env<-rms.olowalu.deep.2017.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.deep.2019.env<-rms.olowalu.deep.2019.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.deep.2020.env<-rms.olowalu.deep.2020.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.deep.2021.env<-rms.olowalu.deep.2021.env[,c('Date', 'Hour', '1.56-3.12 kHz')]

###

rms.olowalu.shallow.2016.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2016-07-30to2017-01-08)RMS_SPL_50kHz.txt", sep=''),50, 2016)
rms.olowalu.shallow.2017.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2017-01-09to2017-05-17)RMS_SPL_50kHz.txt", sep=''),50, 2017)
rms.olowalu.shallow.2018.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2017-12-07to2018-05-17)RMS_SPL_25kHz.txt", sep=''),25, 2017)
rms.olowalu.shallow.2019.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2018-11-02to2019-03-14)RMS_SPL_25kHz.txt", sep=''),25, 2018)
rms.olowalu.shallow.2020.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2019-11-22to2020-05-09)RMS_SPL_25kHz.txt", sep=''),25, 2019)
rms.olowalu.shallow.2021.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2020-12-11to2021-06-18)RMS_SPL_25kHz.txt", sep=''),25, 2020)

rms.olowalu.shallow.2016.env<-rms.olowalu.shallow.2016.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.shallow.2017.env<-rms.olowalu.shallow.2017.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.shallow.2018.env<-rms.olowalu.shallow.2018.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.shallow.2019.env<-rms.olowalu.shallow.2019.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.shallow.2020.env<-rms.olowalu.shallow.2020.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.shallow.2021.env<-rms.olowalu.shallow.2021.env[,c('Date', 'Hour', '1.56-3.12 kHz')]

###

rms.pali.2017.env<-get_rms.EAR(paste(dat.dir, "Pali_EAR(2016-12-22to2017-05-04)RMS_SPL_50kHz.txt", sep=''),50, 2016)

rms.pali.2017.env<-rms.pali.2017.env[,c('Date', 'Hour', '1.56-3.12 kHz')]

###

rms.maui7.2017.env<-get_rms.EAR(paste(dat.dir, "Maui7_EAR(2017-01-04to2017-05-04)RMS_SPL_50kHz.txt", sep=''),50, 2017)

rms.maui7.2017.env<-rms.maui7.2017.env[,c('Date', 'Hour', '1.56-3.12 kHz')]


#---

# combine into one data frame per site ---
rms.maui6<-rbind(rms.maui6.2017, rms.maui6.2018, rms.maui6.2019, rms.maui6.2020, rms.maui6.2021)
rms.olowalu.deep<-rbind(rms.olowalu.deep.2017, rms.olowalu.deep.2019, rms.olowalu.deep.2020,rms.olowalu.deep.2021)
rms.olowalu.shallow<-rbind(rms.olowalu.shallow.2016[,1:3], rms.olowalu.shallow.2017, rms.olowalu.shallow.2018, rms.olowalu.shallow.2019, rms.olowalu.shallow.2020, rms.olowalu.shallow.2021) 
rms.pali<-rms.pali.2017
rms.maui7<-rms.maui7.2017

###

rms.maui6.env<-rbind(rms.maui6.2017.env, rms.maui6.2018.env, rms.maui6.2019.env, rms.maui6.2020.env, rms.maui6.2021.env)
rms.olowalu.deep.env<-rbind(rms.olowalu.deep.2017.env, rms.olowalu.deep.2019.env, rms.olowalu.deep.2020.env, rms.olowalu.deep.2021.env)
rms.olowalu.shallow.env<- rbind(rms.olowalu.shallow.2016.env, rms.olowalu.shallow.2017.env,rms.olowalu.shallow.2018.env, rms.olowalu.shallow.2019.env,rms.olowalu.shallow.2020.env, rms.olowalu.shallow.2021.env)
rms.pali.env<-rms.pali.2017.env
rms.maui7.env<-rms.maui7.2017.env

#---

# add sites ---
rms.maui6$site<-'Maui6'
rms.olowalu.deep$site<-'Olowalu (deep)'
rms.olowalu.shallow$site<-'Olowalu (shallow)'
rms.pali$site<-'Pali'
rms.maui7$site<-'Maui7'

###

rms.maui6.env$site<-'Maui6'
rms.olowalu.deep.env$site<-'Olowalu (deep)'
rms.olowalu.shallow.env$site<-'Olowalu (shallow)'
rms.pali.env$site<-'Pali'
rms.maui7.env$site<-'Maui7'


#---

#combine dataframes
rms.hw<-rbind(rms.maui6, rms.olowalu.deep, rms.olowalu.shallow, rms.pali, rms.maui7)
rms.env<-rbind(rms.maui6.env, rms.olowalu.deep.env, rms.olowalu.shallow.env, rms.pali.env, rms.maui7.env)

#rename "0-1.5 kHz"/"1.56-3.12 kHz" to "SPL" for easier indexing later
names(rms.hw)[names(rms.hw) == '0-1.5 kHz'] <- 'SPL'
names(rms.env)[names(rms.env) == '1.56-3.12 kHz'] <- 'SPL'

### --- --- --- ---

#convert site to factor (important for GAM later)
rms.hw$site<-as.factor(rms.hw$site)
rms.env$site<-as.factor(rms.env$site)

#get julian day and season for HW frequency band
rms.hw$julian<-as.numeric(format(as.Date(rms.hw$Date), "%j"))

rms.hw$season<-year(rms.hw$Date)

#reference whale season relative to December 1st

for (n in 1:nrow(rms.hw)){
  
  season.start<-as.numeric(format(as.Date(paste(rms.hw$season[n],'-12-01', sep='')),'%j'))
  season.start.diff<-ifelse(leap_year(rms.hw$Date[n]),365-as.numeric(format(as.Date(paste(rms.hw$season[n],'-12-01', sep='')),'%j')),364-as.numeric(format(as.Date(paste(rms.hw$season[n],'-12-01', sep='')),'%j')))
  
  if (rms.hw$julian[n]>=season.start){
    rms.hw$day.season[n]<-rms.hw$julian[n]-season.start
  }else{
    rms.hw$day.season[n]<-rms.hw$julian[n]+season.start.diff
  }
  
  n<-n+1
}	

#correct year for December
rms.hw$season[which(rms.hw$julian>334)]<-rms.hw$season[which(rms.hw$julian>334)]+1

#get julian day and season for non-humpback whale environmental frequency band

rms.env$julian<-as.numeric(format(as.Date(rms.env$Date), "%j"))

rms.env$season<-year(rms.env$Date)

for (n in 1:nrow(rms.env)){
  
  season.start<-as.numeric(format(as.Date(paste(rms.env$season[n],'-12-01', sep='')),'%j'))
  season.start.diff<-ifelse(leap_year(rms.env$Date[n]),365-as.numeric(format(as.Date(paste(rms.env$season[n],'-12-01', sep='')),'%j')),364-as.numeric(format(as.Date(paste(rms.env$season[n],'-12-01', sep='')),'%j')))
  
  if (rms.env$julian[n]>=season.start){
    rms.env$day.season[n]<-rms.env$julian[n]-season.start
  }else{
    rms.env$day.season[n]<-rms.env$julian[n]+season.start.diff
  }
  
  n<-n+1
}	

#correct year for December
rms.env$season[which(rms.env$julian>334)]<-rms.env$season[which(rms.env$julian>334)]+1

#get off-season for Olowalu shallow
rms.hw.offseason<-subset(rms.hw, month(Date) %in% c(8,9))
rms.env.offseason<-subset(rms.env, month(Date) %in% c(8,9))

#subset to whale season= December 1st to April 30
rms.hw<-subset(rms.hw, month(Date)==12 | (month(Date)>=1 & month(Date)<=4))
rms.env<-subset(rms.env, month(Date)==12 | (month(Date)>=1 & month(Date)<=4))

#---

#subset to show different patterns early vs peak vs late season
rms.hw$seasontime<-NA
rms.hw$seasontime[month(rms.hw$Date)==12]<-'early' #subset to December 15 - January 15 instead?
rms.hw$seasontime[month(rms.hw$Date)==4]<-'late'
rms.hw$seasontime[month(rms.hw$Date)==2]<-'peak'
rms.hw.sub<-rms.hw[!is.na(rms.hw$seasontime),]
rms.hw.sub$seasontime<-factor(rms.hw.sub$seasontime, levels=c('early', 'peak', 'late'), ordered=T)

rms.env$seasontime<-NA
rms.env$seasontime[month(rms.env$Date)==12]<-'early'
rms.env$seasontime[month(rms.env$Date)==4]<-'late'
rms.env$seasontime[month(rms.env$Date)==2]<-'peak'
rms.env.sub<-rms.env[!is.na(rms.env$seasontime),]
rms.env.sub$seasontime<-factor(rms.env.sub$seasontime, levels=c('early', 'peak', 'late'), ordered=T)

#take out Maui7 and Pali, because we don't have December data anyways
rms.hw.sub<-rms.hw.sub[rms.hw.sub$site %in% c('Maui6', 'Olowalu (deep)', 'Olowalu (shallow)'),]
rms.env.sub<-rms.env.sub[rms.env.sub$site %in% c('Maui6', 'Olowalu (deep)', 'Olowalu (shallow)'),]


#get all dates
sun.times.sub<-as.data.frame(rms.hw.sub%>%
                           group_by(Date)%>%
                           summarize(Date=unique(Date),
                                     season=unique(season),
                                     seasontime=unique(seasontime)))

#get sunrise and sunset times for all dates of the subset data.
#use lat/long of Olowalu (deep) as the reference

#initialize columns for sunrise/sunset times
sun.times.sub$sunrise<-NA
sun.times.sub$sunset<-NA

#Get sunrise and sunset time per date
for(n in 1:nrow(sun.times.sub)){ 

  sun.times.sub$sunrise[n]<-as.character(with(sun.times.sub[n,], getSunlightTimes(date = Date, lat = 20.809, lon = -156.656, keep = c("sunrise", "sunset"), tz = 'HST')[1,'sunrise']))
  sun.times.sub$sunset[n]<-as.character(with(sun.times.sub[n,], getSunlightTimes(date = Date, lat = 20.809, lon = -156.656, keep = c("sunrise", "sunset"), tz = 'HST')[1,'sunset']))
}
  
#convert back to POSIXct object
sun.times.sub$sunrise<-as.POSIXct(sun.times.sub$sunrise, tz='HST')
sun.times.sub$sunset<-as.POSIXct(sun.times.sub$sunset, tz='HST')

#remove date information
sun.times.sub$sunrise<-times(strftime(sun.times.sub$sunrise, "%H:%M:%S", tz="HST"))
sun.times.sub$sunset<-times(strftime(sun.times.sub$sunset, "%H:%M:%S", tz="HST"))

#average sunrise and sunset times per seasontime (early/peak/late)
sun.times<-as.data.frame(sun.times.sub%>%
                           group_by(season,seasontime)%>%
                           summarize(sunrise=median(sunrise),
                                     sunset=median(sunset))%>%
                           group_by(seasontime)%>%
                           summarize(sunrise=mean(sunrise),
                                     sunset=mean(sunset)))
sun.times<-merge(sun.times, c('Maui6', 'Olowalu (shallow)'))
names(sun.times)[4]<-'site'

#---

#get sunrise and sunset times for off-season
sun.times.offseason<-getSunlightTimes(date = as.Date('2016-08-30', tz='HST'), lat = 20.809, lon = -156.656, keep = c("sunrise", "sunset"), tz = 'HST')[1,c('sunrise','sunset')]
sun.times.offseason$sunrise<-times(strftime(sun.times.offseason$sunrise, "%H:%M:%S", tz="HST"))
sun.times.offseason$sunset<-times(strftime(sun.times.offseason$sunset, "%H:%M:%S", tz="HST"))


################################################################


#### 2) SUMMARIZE DATA  ----

#Descriptive statistics pooled across year per month and site and pooled across month and year per site

#rms.hw.summary.month<-as.data.frame(rms.hw%>%
#                                group_by(site, month(Date), Hour)%>%
#                                summarize(SPL.median=median(SPL, na.rm=T, na.action=na.pass),
#                                          SPL.mean=mean(SPL, na.rm=T, na.action=na.pass),
#                                          SPL.sd=sd(SPL, na.rm=T),
#                                          SPL.IQR=IQR(SPL, na.rm=T))
#                              )

#rms.hw.summary.hour<-as.data.frame(rms.hw%>%
#                                group_by(site, Hour)%>%
#                                summarize(SPL.median=median(SPL, na.rm=T, na.action=na.pass),
#                                          SPL.mean=mean(SPL, na.rm=T, na.action=na.pass),
#                                          SPL.sd=sd(SPL, na.rm=T),
#                                          SPL.IQR=IQR(SPL, na.rm=T)))


rms.hw.summary<-as.data.frame(rms.hw%>%
                                group_by(site) %>%
                                summarise(min=min(SPL, na.rm=T),
                                          max=max(SPL, na.rm=T),
                                          median=median(SPL, na.rm=T, na.action.pass),
                                          IQR=IQR(SPL, na.rm=T)))

rms.hw.summary.month<-as.data.frame(rms.hw%>%
                                      group_by(site, month=month(Date), Hour)%>%
                                      summarize(SPL.median=median(SPL, na.rm=T, na.action=na.pass),
                                                SPL.IQR=IQR(SPL, na.rm=T)) %>%
                                      ungroup() %>%
                                      mutate_at(vars(month), funs(as.character(.))) %>%
                                      bind_rows(summarise(month="season",rms.hw%>%group_by(site,Hour),
                                                          SPL.median=median(SPL, na.rm=T, na.action=na.pass),
                                                          SPL.IQR=IQR(SPL, na.rm=T)))%>%
                                      ungroup() %>%
                                      group_by(site,month) %>%
                                      summarise(min=min(SPL.median, na.rm=T), 
                                               max=max(SPL.median, na.rm=T), 
                                               delta=max(SPL.median)-min(SPL.median))  
                                )

#rms.hw.summary<-as.data.frame(rms.hw.summary%>%
#                                group_by(site, year)%>%
#                                summarise(min=min(SPL.median), 
#                                          max=max(SPL.median), 
#                                          delta=max(SPL.median)-min(SPL.median))
#                              )

#write.csv(rms.hw.summary, '../Results/Maui_RMSSPL_summary_hw.csv', row.names=F)
#write.csv(rms.hw.summary.month, '../Results/Maui_RMSSPL_summary_month_hw.csv', row.names=F)

###

rms.env.summary<-as.data.frame(rms.env%>%
                                group_by(site) %>%
                                summarise(min=min(SPL, na.rm=T),
                                          max=max(SPL, na.rm=T),
                                          median=median(SPL, na.rm=T, na.action.pass),
                                          IQR=IQR(SPL, na.rm=T)))

rms.env.summary.month<-as.data.frame(rms.env%>%
                                      group_by(site, month=month(Date), Hour)%>%
                                      summarize(SPL.median=median(SPL, na.rm=T, na.action=na.pass),
                                                SPL.IQR=IQR(SPL, na.rm=T)) %>%
                                      ungroup() %>%
                                      mutate_at(vars(month), funs(as.character(.))) %>%
                                      bind_rows(summarise(month="season",rms.env%>%group_by(site,Hour),
                                                          SPL.median=median(SPL, na.rm=T, na.action=na.pass),
                                                          SPL.IQR=IQR(SPL, na.rm=T)))%>%
                                      ungroup() %>%
                                      group_by(site,month) %>%
                                      summarise(min=min(SPL.median, na.rm=T), 
                                                max=max(SPL.median, na.rm=T), 
                                                delta=max(SPL.median)-min(SPL.median)) 
                                      )

#write.csv(rms.env.summary, 'Results/Maui_RMSSPL_summary_env.csv', row.names=F)
#write.csv(rms.env.summary.month, 'Results/Maui_RMSSPL_summary_month_env.csv', row.names=F)


################################################################


#### 3) STATISTICAL ANALYSIS ----

#Fit generalized additive models (GAMs) for both frequency bands

# Humpback whale band

gam.rms.hw<-gam(SPL ~ s(day.season, k=15) + s(Hour, bs='cc',by=site) + site + as.factor(season), data=rms.hw, method='REML', select=T)
#gam.rms.full.hw<-gam(SPL ~ s(day.season) + te(day.season, Hour, by=site, bs=c('tp', 'cc')) + s(Hour, bs='cc',by=site) + site + as.factor(season), data=rms.hw, method='REML', select=T)
#We do not use the full model to reduce the number of model terms and avoid overfitting

#check if model is appropriate
gam.check(gam.rms.hw, old.style=T)

#get signifcances
anova(gam.rms.hw)

#check for autocorrelation 
acf(residuals(gam.rms.hw))
pacf(residuals(gam.rms.hw))

#There is considerable autocorrelation, so we're refitting the model as a GAMM:

gamm.rms.hw<-gamm(SPL ~ s(day.season, k=18) + s(Hour, bs='cc',by=site) + site + as.factor(season), 
                                    data=rms.hw, 
                                    method='REML',
                                    correlation=corAR1(form=~1|Date+site))

#check if model is appropriate
gam.check(gamm.rms.hw$gam, old.style=T)

#check for autocorrelation 
acf(resid(gamm.rms.hw$lme, type='normalized'))
pacf(resid(gamm.rms.hw$lme, type='normalized'))

#get signifcances
anova(gamm.rms.hw$gam)


# Non-humpback whale environmental band

gam.rms.env<-gam(SPL ~ s(day.season, k=14) + s(Hour, bs='cc', by=site) + site + as.factor(season), data=rms.env, method='REML', select=T)
#gam.rms.env.full<-gam(SPL ~ s(day.season) + te(day.season, Hour, by=site, bs=c('tp', 'cc')) + s(Hour, bs='cc',by=site) + site + as.factor(season), data=rms.env, method='REML', select=T)
#We do not use the full model to reduce the number of model terms and avoid overfitting

#check if model is appropriate
gam.check(gam.rms.env, old.style=T)

#get signifcances
anova(gam.rms.env)

#check for autocorrelation 
acf(residuals(gam.rms.hw))
pacf(residuals(gam.rms.hw))

#There is considerable autocorrelation, so we're refitting the model as a GAMM:

gamm.rms.env<-gamm(SPL ~ s(day.season, k=18) + s(Hour, bs='cc',by=site) + site + as.factor(season), 
                  data=rms.env, 
                  method='REML',
                  correlation=corAR1(form=~1|Date+site))

#check if model is appropriate
gam.check(gamm.rms.hw$gam, old.style=T)

#check for autocorrelation 
acf(resid(gamm.rms.env$lme, type='normalized'))
pacf(resid(gamm.rms.env$lme, type='normalized'))

#get signifcances
anova(gamm.rms.env$gam)


#---

#Create effect plots. Because plot.gam() doesn't allow for a lot of adjustments of the look of the plots,
#we're going to create them from scratch.

# (a) - pull out all smoothers 

#smoother.dat.hw <- plot(gam.rms.hw, pages=1, scale=0)  # plot.gam returns a list of n elements, one per plot
smoother.dat.hw <- plot(gamm.rms.hw$gam, pages=1, scale=0) 
#smoother.dat.env <- plot(gam.rms.env, pages=1, scale=0)  
smoother.dat.env <- plot(gamm.rms.env$gam, pages=1, scale=0)  

#get x, fit, and se for each effect plot

# Humpback whale band

sm.seasonal.hw<-get_GAMsmoothers(smoother.dat.hw,1)
sm.diel.maui6.hw<-get_GAMsmoothers(smoother.dat.hw,2)
sm.diel.maui7.hw<-get_GAMsmoothers(smoother.dat.hw,3)
sm.diel.olo_s.hw<-get_GAMsmoothers(smoother.dat.hw,4)
sm.diel.olo_d.hw<-get_GAMsmoothers(smoother.dat.hw,5)
sm.diel.pali.hw<-get_GAMsmoothers(smoother.dat.hw,6)


# Non-humpback whale environmental band

sm.seasonal.env<-get_GAMsmoothers(smoother.dat.env,1)
sm.diel.maui6.env<-get_GAMsmoothers(smoother.dat.env,2)
sm.diel.maui7.env<-get_GAMsmoothers(smoother.dat.env,3)
sm.diel.olo_s.env<-get_GAMsmoothers(smoother.dat.env,4)
sm.diel.olo_d.env<-get_GAMsmoothers(smoother.dat.env,5)
sm.diel.pali.env<-get_GAMsmoothers(smoother.dat.env,6)


###---

# (b) - get parametric terms

#para.dat.hw <- as.data.frame(termplot(gam.rms.hw, se = TRUE, plot = FALSE))
para.dat.hw <- as.data.frame(termplot(gamm.rms.hw$gam, se = TRUE, plot = FALSE))
#para.dat.env <- as.data.frame(termplot(gam.rms.env, se = TRUE, plot = FALSE))
para.dat.env <- as.data.frame(termplot(gamm.rms.env$gam, se = TRUE, plot = FALSE))


###---

# (c) - create individual plots and save as objects

# Humpback whale band

gam.plot.seasonal.hw<-ggplot(data=sm.seasonal.hw, aes(x=x, y=fit)) +
  gam_plots_params(seasonal=T, hw=T) 

gam.plot.diel.maui6.hw<-ggplot(data=sm.diel.maui6.hw, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, hw=T, site='Maui6')

gam.plot.diel.maui7.hw<-ggplot(data=sm.diel.maui7.hw, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, hw=T, site='Maui7')

gam.plot.diel.olo_d.hw<-ggplot(data=sm.diel.olo_s.hw, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, hw=T, site='Olowalu (deep)')

gam.plot.diel.olo_s.hw<-ggplot(data=sm.diel.olo_d.hw, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, hw=T, site='Olowalu (shallow)')

gam.plot.diel.pali.hw<-ggplot(data=sm.diel.pali.hw, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, hw=T, site='Pali')

gam.paraplot.sites.hw<-ggplot(data=para.dat.hw, aes(x=site.x, y=site.y)) +
  gam_paraplots_params(site=T, hw=T) 

gam.paraplot.seasons.hw<-ggplot(data=para.dat.hw, aes(x=season.x, y=season.y)) +
  gam_paraplots_params(season=T, hw=T) 



# Non-humpback whale environmental band

gam.plot.seasonal.env<-ggplot(data=sm.seasonal.env, aes(x=x, y=fit)) +
  gam_plots_params(seasonal=T, env=T) 

gam.plot.diel.maui6.env<-ggplot(data=sm.diel.maui6.env, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, env=T, site='Maui6') 

gam.plot.diel.maui7.env<-ggplot(data=sm.diel.maui7.env, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, env=T, site='Maui7')

gam.plot.diel.olo_d.env<-ggplot(data=sm.diel.olo_s.env, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, env=T, site='Olowalu (deep)')

gam.plot.diel.olo_s.env<-ggplot(data=sm.diel.olo_d.env, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, env=T, site='Olowalu (shallow)')

gam.plot.diel.pali.env<-ggplot(data=sm.diel.pali.env, aes(x=x, y=fit)) +
  gam_plots_params(diel=T, env=T, site='Pali') 

gam.paraplot.sites.env<-ggplot(data=para.dat.env, aes(x=site.x, y=site.y)) +
  gam_paraplots_params(site=T, env=T) 

gam.paraplot.seasons.env<-ggplot(data=para.dat.env, aes(x=season.x, y=season.y)) +
  gam_paraplots_params(season=T, env=T)

### ---

#combine subplots into one plot and save to file (Fig. 6 in Kügler et al.)

gam.plots.hw<-(gam.plot.seasonal.hw | gam.plot.diel.olo_s.hw | gam.plot.diel.maui7.hw | gam.plot.diel.olo_d.hw) /
                  (gam.plot.diel.pali.hw | gam.plot.diel.maui6.hw | gam.paraplot.sites.hw | gam.paraplot.seasons.hw)

gam.plots.hw<-(gam.plot.seasonal.hw | gam.paraplot.sites.hw | gam.paraplot.seasons.hw | plot_spacer() | plot_spacer()) /
  (gam.plot.diel.olo_s.hw | gam.plot.diel.maui7.hw | gam.plot.diel.olo_d.hw | gam.plot.diel.pali.hw | gam.plot.diel.maui6.hw )


gam.plots.env<-(gam.plot.seasonal.env | gam.plot.diel.olo_s.env | gam.plot.diel.maui7.env | gam.plot.diel.olo_d.env) /
  (gam.plot.diel.pali.env | gam.plot.diel.maui6.env | gam.paraplot.sites.env | gam.paraplot.seasons.env)


#ggsave(gam.plots.hw, filename='../Plots/Maui_RMSSPL_GAM_hw.png', width=18, height=7, units='in',  bg = "white")    
#ggsave(gam.plots.env, filename='../Plots/Maui_RMSSPL_GAM_env.png', width=18, height=7, units='in',  bg = "white")    


### ---------------------------------------------

#3.2 - early/peak/late for deep vs. shallow recorder

#to illustrate diel patterns during different times of the season corresponding with different whale densities
#using data from one offshore (Maui6) and one nearshore (Olowalu shallow) site as representative examples
#Showing the smooth fits from GAM(M)s fit to a subset of the data for both frequency bands


### ---

#Fit generalized additive (mixed) models for both frequency bands

# Humpback whale band

gam.rms.hw.sub<-gam(SPL ~ s(day.season, k=10) + s(Hour, bs='cc', by=interaction(site, seasontime))  + site + as.factor(season), 
                    data=rms.hw.sub, method='REML', select=T)  

gam.check(gam.rms.hw.sub, old.style=T)

anova(gam.rms.hw.sub)

#check for autocorrelation
acf(residuals(gam.rms.hw.sub))
pacf(residuals(gam.rms.hw.sub))

#There is considerable autocorrelation, so we're refitting as a GAMM:

gamm.rms.hw.sub<-gamm(SPL ~ s(day.season, k=10) + s(Hour, bs='cc', by=interaction(site, seasontime)) + site + as.factor(season), 
                      data=rms.hw.sub, 
                      method='REML',
                      correlation=corAR1(form=~1|Date+site))

#check for autocorrelation 
acf(resid(gamm.rms.hw.sub$lme, type='normalized'))
pacf(resid(gamm.rms.hw.sub$lme, type='normalized'))


# Non-humpback whale environmental band

gam.rms.env.sub<-gam(SPL ~ s(day.season, k=10) + s(Hour, bs='cc', by=interaction(site, seasontime)) + site + as.factor(season), 
                     data=rms.env.sub, method='REML', select=T) 

#check for autocorrelation
acf(residuals(gam.rms.env.sub))
pacf(residuals(gam.rms.env.sub))

#There is considerable autocorrelation, so we're refitting as a GAMM:

gamm.rms.env.sub<-gamm(SPL ~ s(day.season, k=10) + s(Hour, bs='cc', by=interaction(site, seasontime)) + site + as.factor(season), 
                       data=rms.env.sub, 
                       method='REML',
                       correlation=corAR1(form=~1|Date+site))

#check for autocorrelation 
acf(resid(gamm.rms.env.sub$lme, type='normalized'))
pacf(resid(gamm.rms.env.sub$lme, type='normalized'))


#---

#Create effect plots. 

# (a) - pull out all smoothers 

#smoother.dat.hw.sub <- plot(gam.rms.hw.sub, pages=1, scale=0)  
smoother.dat.hw.sub <- plot(gamm.rms.hw.sub$gam, pages=1, scale=0)  
#smoother.dat.env.sub <- plot(gam.rms.env.sub, pages=1, scale=0)
smoother.dat.env.sub <- plot(gamm.rms.env.sub$gam, pages=1, scale=0) 

#get x, fit, and se for each effect plot (only pull out diel smoothers)

# Humpback whale band

sm.diel.maui6.early<-get_GAMsmoothers(smoother.dat.hw.sub,2)
sm.diel.olo_s.early<-get_GAMsmoothers(smoother.dat.hw.sub,4)
sm.diel.maui6.peak<-get_GAMsmoothers(smoother.dat.hw.sub,5)
sm.diel.olo_s.peak<-get_GAMsmoothers(smoother.dat.hw.sub,7)
sm.diel.maui6.late<-get_GAMsmoothers(smoother.dat.hw.sub,8)
sm.diel.olo_s.late<-get_GAMsmoothers(smoother.dat.hw.sub,10)

#add a factor for seasontime
sm.diel.maui6.early$seasontime<-'early'
sm.diel.maui6.peak$seasontime<-'peak'
sm.diel.maui6.late$seasontime<-'late'
#combine into one dataframe
sm.diel.maui6<-rbind(sm.diel.maui6.early, sm.diel.maui6.peak, sm.diel.maui6.late)
#add site factor
sm.diel.maui6$site<-'Maui6'

#add a factor for seasontime
sm.diel.olo_s.early$seasontime<-'early'
sm.diel.olo_s.peak$seasontime<-'peak'
sm.diel.olo_s.late$seasontime<-'late'
#combine into one dataframe
sm.diel.olo_s<-rbind(sm.diel.olo_s.early, sm.diel.olo_s.peak, sm.diel.olo_s.late)
#add site factor
sm.diel.olo_s$site<-'Olowalu (shallow)'

#combine dataframes for both sites
sm.diel.hw.sub<-rbind(sm.diel.maui6, sm.diel.olo_s)
#add factor for frequency band
sm.diel.hw.sub$fband<-'HW'

# Non-humpback whale environmental band

sm.diel.env.maui6.early<-get_GAMsmoothers(smoother.dat.env.sub,2)
sm.diel.env.olo_s.early<-get_GAMsmoothers(smoother.dat.env.sub,4)
sm.diel.env.maui6.peak<-get_GAMsmoothers(smoother.dat.env.sub,5)
sm.diel.env.olo_s.peak<-get_GAMsmoothers(smoother.dat.env.sub,7)
sm.diel.env.maui6.late<-get_GAMsmoothers(smoother.dat.env.sub,8)
sm.diel.env.olo_s.late<-get_GAMsmoothers(smoother.dat.env.sub,10)

#add factor for seasontime
sm.diel.env.maui6.early$seasontime<-'early'
sm.diel.env.maui6.peak$seasontime<-'peak'
sm.diel.env.maui6.late$seasontime<-'late'
#combine into one dataframe
sm.diel.env.maui6<-rbind(sm.diel.env.maui6.early, sm.diel.env.maui6.peak, sm.diel.env.maui6.late)
#add site factor
sm.diel.env.maui6$site<-'Maui6'

#add factor for seasontime
sm.diel.env.olo_s.early$seasontime<-'early'
sm.diel.env.olo_s.peak$seasontime<-'peak'
sm.diel.env.olo_s.late$seasontime<-'late'
#combine into one dataframe
sm.diel.env.olo_s<-rbind(sm.diel.env.olo_s.early, sm.diel.env.olo_s.peak, sm.diel.env.olo_s.late)
#add site factor
sm.diel.env.olo_s$site<-'Olowalu (shallow)'

#combine dataframes for both sites
sm.diel.env.sub<-rbind(sm.diel.env.maui6, sm.diel.env.olo_s)
#add factor for frequency band
sm.diel.env.sub$fband<-'Env'

###

#combine both frequency bands into one dataframe
sm.diel.sub<-rbind(sm.diel.hw.sub, sm.diel.env.sub)

###---

# (b) - create plots and save as objects

#create dummy data for a blank plot because facet_grid doesn't allow to define different axis scales
dummy <- data.frame(fit = c(-5,15, -2.5, 2.5), x = 0, seasontime=c('early', 'peak'),
                    site=c('Maui6', 'Maui6', 'Olowalu (shallow)', 'Olowalu (shallow)'), 
                    stringsAsFactors=FALSE)

#plot diel patterns per seasontime per site for the humpback whale frequency band
gam.diel.plots.sub.hw<-ggplot(sm.diel.sub[sm.diel.sub$fband=='HW',], aes(x=x, y=fit)) + 
  gam_plots_sub.diel_params() +
  ylim(-2.1,2.1) 

#plot diel patterns per seasontime per site for both frequency bands in overlay plots
gam.diel.plots.sub.both<-ggplot(sm.diel.sub, aes(x=x, y=fit)) +
  gam_plots_sub.diel_params(ribbon=F, int=T) +
  geom_blank(data=dummy)

###---

#combine into one plot and save to file (Fig. S3 in Kügler et al.)

gam.diel.plots.sub<- gam.diel.plots.sub.hw / gam.diel.plots.sub.both  +
  plot_annotation(tag_levels = c('A', 'B')) &
  theme(plot.tag = element_text(size = 20),
        plot.tag.position = c(0, 0.95)) 

#ggsave(gam.diel.plots.sub, filename='../Plots/Maui_RMSSPL_GAM_diel_seasontime.png', width=15, height=12, units='in',  bg = "white") 

### ---------------------------------------------

#3.3 - off-season

#To illustrate seasonal and diel patterns during the off-season at the Olowalu (shallow) site for both frequency bands
#using data from the months August and September of 2016. We do not have off-season data for other years.

#---

#Fit generalized additive models for both frequency bands

gam.rms.hw.offseason<-gam(SPL ~ s(day.season) + s(Hour, bs='cc'), data=rms.hw.offseason, method='REML', select=T)
gam.rms.env.offseason<-gam(SPL ~ s(day.season) + s(Hour, bs='cc') , data=rms.env.offseason, method='REML', select=T)

#Check for autocorrelation
acf(residuals(gam.rms.hw.offseason))
acf(residuals(gam.rms.env.offseason))

#---

#Create effect plots. 

# (a) - pull out all smoothers 

smoother.dat.hw.offseason <- plot(gam.rms.hw.offseason, pages=1, scale=0)
smoother.dat.env.offseason <- plot(gam.rms.env.offseason, pages=1, scale=0)

#get x, fit, and se for each effect plot

# Humpback whale band
sm.seasonal.hw.offseason<-get_GAMsmoothers(smoother.dat.hw.offseason,1)
sm.diel.hw.offseason<-get_GAMsmoothers(smoother.dat.hw.offseason,2)

# Non-humpback whale environmental band
sm.seasonal.env.offseason<-get_GAMsmoothers(smoother.dat.env.offseason,1)
sm.diel.env.offseason<-get_GAMsmoothers(smoother.dat.env.offseason,2)

#---

# (b) - create plots and save as objects

# Humpback whale band

gam.plot.seasonal.hw.offseason<-ggplot(data=sm.seasonal.hw.offseason, aes(x=x, y=fit)) +
  gam_plots_params(seasonal=T, hw=T, se=F) +
  scale_y_continuous(limits=c(-4,4)) +
  scale_x_continuous(limits=c(240,303), breaks=c(240,260,280,300))

gam.plot.diel.hw.offseason<-ggplot(data=sm.diel.hw.offseason, aes(x=x, y=fit)) +
  geom_rect(data=sun.times.offseason, aes(xmin = 0, xmax = sunrise, ymin = -Inf, ymax = Inf),alpha = 0.42, fill = "grey85", inherit.aes=F) + #dodgerblue4
  geom_rect(data=sun.times.offseason, aes(xmin = sunset, xmax = Inf, ymin = -Inf, ymax = Inf),alpha = 0.42, fill = "grey85", inherit.aes=F) +
  gam_plots_params(diel=T, hw=T, se=F) +
  scale_y_continuous(limits=c(-4,4))

# Non-humpback whale environmental band

gam.plot.seasonal.env.offseason<-ggplot(data=sm.seasonal.env.offseason, aes(x=x, y=fit)) +
  gam_plots_params(seasonal=T, env=T, se=F) +
  scale_y_continuous(limits=c(-4,4)) +
  scale_x_continuous(limits=c(240,303), breaks=c(240,260,280,300))

gam.plot.diel.env.offseason<-ggplot(data=sm.diel.env.offseason, aes(x=x, y=fit)) +
  geom_rect(data=sun.times.offseason, aes(xmin = 0, xmax = sunrise, ymin = -Inf, ymax = Inf),alpha = 0.42, fill = "grey85", inherit.aes=F) + 
  geom_rect(data=sun.times.offseason, aes(xmin = sunset, xmax = Inf, ymin = -Inf, ymax = Inf),alpha = 0.42, fill = "grey85", inherit.aes=F) +
  gam_plots_params(diel=T, env=T, se=F) +
  scale_y_continuous(limits=c(-4,4))

#---

#combine into one plot and save to file (Fig. S4 in Kügler et al.)

gam.plots.offseason<-((gam.plot.seasonal.hw.offseason | gam.plot.diel.hw.offseason ) /
                        (gam.plot.seasonal.env.offseason | gam.plot.diel.env.offseason ) ) + plot_annotation(title = "Olowalu (shallow) - offseason",
                                                                                                             theme = theme(plot.title = element_text(size = 16)))

#ggsave(gam.plots.offseason, filename='../Plots/Maui_RMSSPL_GAM_offseason.png', width=10, height=8, units='in',  bg = "white")


################################################################


### 4) POOL DATA FOR VISUALIZATIONS/POLAR PLOTS ----

#All hourly median RMS SPL values are normalized by the respective mean per season and site per frequency band
#Then, for each site and frequency band, the median per hour is calculated
#Medians for the humpback whale band are scaled between 0-1 for each site
#Medians for the non-humpback whale band are scaled between 0-1 across sites to illustrate the differing magnitudes of the patterns
#corresponds with Fig. 5 in Kügler et al.

#normalize hourly median SPL per site and year
rms.hw<-as.data.frame(rms.hw %>%
                        group_by(season,site) %>%
                        mutate(SPL.norm = normalize_fun(SPL)
                        ))

rms.hw.season<-as.data.frame(rms.hw[!is.na(rms.hw$seasontime),] %>%
                               group_by(season, site, seasontime) %>%
                               mutate(SPL.scaled = scale_fun(SPL.norm)))
  
rms.env<-as.data.frame(rms.env %>%
                         group_by(season, site) %>%
                         mutate(SPL.norm = normalize_fun(SPL)
                         ))

#calculate median normalized SPL per hour per site
rms.hw.h <- as.data.frame(rms.hw %>%
                            group_by(site,Hour) %>%
                            summarise(SPL.norm.median = median(SPL.norm)
                            ))

rms.env.h <- as.data.frame(rms.env %>%
                             group_by(site,Hour) %>%
                             summarise(SPL.norm.median = median(SPL.norm)
                             ))

#scale medians per site (hw) or across sites (non-hw)
rms.hw.h<-as.data.frame(rms.hw.h %>%
                          group_by(site) %>%
                          mutate(SPL.norm.median.scaled = scale_fun(SPL.norm.median)
                          ))

rms.env.h<-as.data.frame(rms.env.h %>%
                            mutate(SPL.norm.median.scaled = scale_fun(SPL.norm.median)
                            ))

### Diel polar plots --- (Fig. 5 in Kügler et al.)

# humpback whale band

rms.hw.h$Hour_label<-as.POSIXct(as.character(rms.hw.h$Hour), format = "%H:%M:%S", tz='UTC') #force UTC for easier plotting

diel.plots<-ggplot(rms.hw.h, aes(x = Hour_label, y=SPL.norm.median.scaled, fill=SPL.norm.median.scaled)) +
  polarplots_params(hw=T)

#ggsave(diel.plots, filename="../Plots/Maui_RMSSPL_diel_hw.png", width=17.5, height=7.5, units='in', bg = "white")    

# Non-humpback whale environmental band

rms.env.h$Hour_label<-as.POSIXct(as.character(rms.env.h$Hour), format = "%H:%M:%S", tz='UTC') #force UTC for easier plotting

diel.env.plots<-ggplot(rms.env.h, aes(x = Hour_label, y=SPL.norm.median.scaled, fill=SPL.norm.median.scaled)) +
  polarplots_params(env=T)

#ggsave(diel.env.plots, filename="../Plots/Maui_RMSSPL_diel_env.png", width=17.5, height=7.5, units='in', bg = "white")    


################################################################


#### 5) OCTAVE BANDS TIME-SERIES PLOTS FOR VISUALIZATION ----

#to illustrate seasonal and diel patterns for different frequency bands
#using data from one season of one offshore (Maui6) and one nearshore (Olowalu deep) site as representative examples
#daily averages are further smoothed with a running average of n=10
#corresponds with Fig. 2 in Kügler et al.

#get the data
octaves.maui6.2019<-get_rms.EAR(paste(dat.dir, "Maui6_EAR(2018-12-01to2019-04-24)RMS_SPL_50kHz.txt", sep=''),50, 2018)
octaves.olowalu.deep.2019<-get_rms.EAR(paste(dat.dir, "Olowalu_Deep_EAR(2018-12-01to2019-04-25)RMS_SPL_50kHz.txt", sep=''),50, 2018)

#average by day (seasonal) and hour (diel)
octaves.maui6.day<-aggregate(.~Date, data = octaves.maui6.2019[names(octaves.maui6.2019)!='Hour'], median, na.action = na.pass)
octaves.maui6.h<-aggregate(.~Hour, data = octaves.maui6.2019[names(octaves.maui6.2019)!='Date'], median, na.rm=T, na.action=na.pass)

octaves.olowalu.day<-aggregate(.~Date, data = octaves.olowalu.deep.2019[names(octaves.olowalu.deep.2019)!='Hour'], median, na.action = na.pass)
octaves.olowalu.h<-aggregate(.~Hour, data = octaves.olowalu.deep.2019[names(octaves.olowalu.deep.2019)!='Date'], median, na.rm=T, na.action=na.pass)

#convert dataframes into long format
octaves.maui6.day.long<-melt(octaves.maui6.day, id.vars='Date')
octaves.olowalu.day.long<-melt(octaves.olowalu.day, id.vars='Date')

octaves.maui6.h.long<-melt(octaves.maui6.h, id.vars='Hour')
octaves.maui6.h.long$Hour<-as.POSIXct(as.character(octaves.maui6.h$Hour), format = "%H:%M:%S", tz='UTC') #force UTC for easier plotting

octaves.olowalu.h.long<-melt(octaves.olowalu.h, id.vars='Hour')
octaves.olowalu.h.long$Hour<-as.POSIXct(as.character(octaves.olowalu.h$Hour), format = "%H:%M:%S", tz='UTC')

#--- 

#axis labels for the ggplots
lab.y.day<-bquote("Median RMS SPL"~'[dB re 1'~plain("\u03BC")*'Pa]')
lab.y.h<-bquote("Median RMS SPL"~'[dB re 1'~plain("\u03BC")*'Pa]')

#create color palette for different frequency bands, use modified viridis color palette
pkgTest('viridis')
cols<-viridis(6, direction=-1)
cols[1]<-"#FCA50AFF" #remove yellow because it's not very visible
names(cols) <- as.character(names(octaves.maui6.day[-1]))

#---  

#plot the seasonal time-series

octaves.maui6.plot.day<-ggplot(data=octaves.maui6.day.long, aes(x=Date, y=value, shape=variable)) +
  octaves_plots_params(seasonal=T, site='Maui6') +
  theme(plot.title = element_text(hjust=1))

octaves.olowalu.plot.day<-ggplot(data=octaves.olowalu.day.long, aes(x=Date, y=value, shape=variable)) +
  octaves_plots_params(seasonal=T, site='Olowalu (deep)')

#---  

#plot diel averages

octaves.maui6.plot.h<-ggplot(data=octaves.maui6.h.long, aes(x=Hour, y=value, shape=variable)) + 
  octaves_plots_params(diel=T, site='Maui6')  +
  theme(plot.title = element_text(hjust=1))

octaves.olowalu.plot.h<-ggplot(data=octaves.olowalu.h.long, aes(x=Hour, y=value, shape=variable)) + 
  octaves_plots_params(diel=T, site='Olowalu (deep)')

#---

#arrange subplots into single plot and save to file 

Maui.Octaves<-((octaves.olowalu.plot.day | octaves.maui6.plot.day) /
  (octaves.olowalu.plot.h | octaves.maui6.plot.h)) + plot_layout(guides = "collect") & theme(legend.position = "top") 

#ggsave(Maui.Octaves, filename='../Plots/Maui_RMSSPL_Octaves_Comparison.png', width=20, height=12, units='in',  bg = "white")    


################################################################

#### ----

#End of File

# -------

################################################################
