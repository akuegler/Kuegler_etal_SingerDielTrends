#################################################################
#################################################################

### DIEL PATTERNS OF HUMPBACK WHALE CHORUSING OFF MAUI, HAWAIʻI

# Written by Anke Kügler
# Last updated by Anke Kügler, 06/10/2022
# Contact: anke.kuegler@gmail.com

# This script analyses hourly averages of acoustic recordings (RMS SPL) 
# collected from autonomous passive acoustic recorders (Ecological Acoustic Recorders, EARs)
# off Maui, Hawaiʻi between 2017 and 2021 to investigate diel patterns of male humpback whale (HW) chorusing.
# In addition, ambient noise in the adjacent frequency band to HW chorusing is analysed to investigate
# possible impacts of environmental noise. 

# Input data: text files containing RMS SPL values for each 30-second recording per deployment and site

# Levels are analysed in two frequency bands:
# 1: 0-1.5 kHz obtained from subsampled (3 kHz) recordings, associated with humpback whale chorusing
# 2: 1.56-3.23 kHz obtained as octave bands from the fullband recordings, adjacent to HW band

#Generalized Additive Models (GAMs) are fit for hourly averages for both frequency bands and 
#   day of the season
#   hour of the day
#   season (representing respective whale seasons from December through April)
#   site

#Functions:
# pkgTest() is a function that tests if a R package already exists on a system and if not, installs and loads it
# get_rms.EAR() reads in the raw RMS SPL data from file, creates date and hour timestamps, corrects dates for
#   leap years if necessary, and averages (mean or median) the data per hour per day.
# gam_plots_params() and gam_paraplots_params() define settings for the GAM effect plots
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
get_rms.EAR<-function(file=NULL, fs=NULL, year=NULL, site=NULL, round_hour=F, ...){
  
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
  rms$Hour<-round_date(rms$Hour, unit="hour")
  rms$Hour<-times(strftime(rms$Hour, "%H:%M:%S", tz="UTC"))
  
  rms.h<-aggregate(.~Hour+Date, data=rms, median, na.action=na.omit) 
  
  
  ####
  
  
  #RETURN FROM FUNCTION
  
  #Save to object/list for subsequent analysis
  
  return(rms.h[,-c(3,4)]) #take out "File"  and Time columns
  
} 

#---

#GAM plots settings

gam_plots_params<-function(seasonal=F, diel=F,site=NULL, hw=F, env=F){
  list(
    geom_ribbon(aes(ymin=fit-se, ymax=fit+se, group = 1), 
                fill = "grey75", alpha = 0.6),
    geom_line(color='black'),
    theme_bw(),
    theme(axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 15, color = "black"),
          axis.text.y = element_text(size = 15, color = "black", angle = 90),
          plot.title=element_text(size = 16, color = "black")),
    theme(panel.border=element_rect(color="black", size=0.75), 
          panel.grid = element_blank()),
    theme(axis.line.x.top = element_line(color="black", size=0.75)),
    labs(y="Partial effect for RMS SPL"),
    if(diel)scale_x_continuous(breaks=seq(from=0, to=23, by=3)/24, labels=c('0', '3', '6', '9', '12', '15', '18', '21')),
    if(diel) labs(x="Hour", title=site),
    if(seasonal)labs(x="Day of the season", title=""),
    if(env) scale_y_continuous(limits=c(-3.5, 11), breaks=c(-2, 0, 2, 4, 6, 8, 10)), #, labels=c('','0','','4','','8','')
    if(seasonal & hw)scale_y_continuous(limits=c(-14, 7), breaks=c(-10, -5, 0, 5)),
    if(diel & hw)scale_y_continuous(limits=c(-2.5, 2.5), breaks=c(-2, -1, 0, 1,2))
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
    theme(axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 15, color = "black"),
          axis.text.y = element_text(size = 15, color = "black", angle = 90),
          plot.title=element_text(size = 16, color = "black")),
    theme(panel.border=element_rect(color="black", size=0.75), 
          panel.grid = element_blank()),
    theme(axis.line.x.top = element_line(color="black", size=0.75)),
    labs(x='', y="Partial effect for RMS SPL", title=''),
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)),
    if(env) scale_y_continuous(limits=c(-3.5, 11), breaks=c(-2, 0, 2, 4, 6, 8, 10)), #, labels=c('','0','','4','','8','')
    if(hw) scale_y_continuous(limits=c(-2, 9), breaks=c(-2,0,2,4,6,8))
  )
}

#---

#Polarplots and octave band comparison plots settings

polarplots_params<-function(hw=F, env=F){
  list(
    theme_minimal(),
      facet_wrap(~site, nrow=1),
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
            strip.text = element_text (size=14, color='black', margin = margin (t=0, r=0, b=0.5, l=0, 'lines')),
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
      if(seasonal)theme(axis.title.y = element_text(size = 16.5, color="black")),
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

rms.olowalu.deep.2017<-get_rms.EAR(paste(dat.dir, "Olowalu_Deep_EAR(2016-12-22to2017-05-04)RMS_SPL_3kHz.txt", sep=''),3, 2016)
rms.olowalu.deep.2019<-get_rms.EAR(paste(dat.dir, "Olowalu_Deep_EAR(2018-12-01to2019-04-25)RMS_SPL_3kHz.txt", sep=''),3, 2018)
rms.olowalu.deep.2020<-get_rms.EAR(paste(dat.dir, "Olowalu_Deep_EAR(2019-11-06to2020-03-30)RMS_SPL_3kHz.txt", sep=''),3, 2019)
rms.olowalu.deep.2021<-get_rms.EAR(paste(dat.dir, "Olowalu_Deep_EAR(2020-12-01to2021-05-26)RMS_SPL_3kHz.txt", sep=''),3, 2020)

###

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

rms.olowalu.deep.2017.env<-get_rms.EAR(paste(dat.dir, "Olowalu_Deep_EAR(2016-12-22to2017-05-04)RMS_SPL_50kHz.txt", sep=''),50, 2016)
rms.olowalu.deep.2019.env<-get_rms.EAR(paste(dat.dir, "Olowalu_Deep_EAR(2018-12-01to2019-04-25)RMS_SPL_50kHz.txt", sep=''),50, 2018)
rms.olowalu.deep.2020.env<-get_rms.EAR(paste(dat.dir, "Olowalu_Deep_EAR(2019-11-06to2020-03-30)RMS_SPL_25kHz.txt", sep=''),25, 2019)
rms.olowalu.deep.2021.env<-get_rms.EAR(paste(dat.dir, "Olowalu_Deep_EAR(2020-12-01to2021-05-26)RMS_SPL_25kHz.txt", sep=''),25, 2020)

rms.olowalu.deep.2017.env<-rms.olowalu.deep.2017.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.deep.2019.env<-rms.olowalu.deep.2019.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.deep.2020.env<-rms.olowalu.deep.2020.env[,c('Date', 'Hour', '1.56-3.12 kHz')]
rms.olowalu.deep.2021.env<-rms.olowalu.deep.2021.env[,c('Date', 'Hour', '1.56-3.12 kHz')]

###

rms.olowalu.shallow.2017.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2017-01-09to2017-05-17)RMS_SPL_50kHz.txt", sep=''),50, 2017)
rms.olowalu.shallow.2018.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2017-12-07to2018-05-17)RMS_SPL_25kHz.txt", sep=''),25, 2017)
rms.olowalu.shallow.2019.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2018-11-02to2019-03-14)RMS_SPL_25kHz.txt", sep=''),25, 2018)
rms.olowalu.shallow.2020.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2019-11-22to2020-05-09)RMS_SPL_25kHz.txt", sep=''),25, 2019)
rms.olowalu.shallow.2021.env<-get_rms.EAR(paste(dat.dir, "Olowalu_shallow_EAR(2020-12-11to2021-06-18)RMS_SPL_25kHz.txt", sep=''),25, 2020)

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
rms.olowalu.shallow<-rbind(rms.olowalu.shallow.2017, rms.olowalu.shallow.2018, rms.olowalu.shallow.2019, rms.olowalu.shallow.2020, rms.olowalu.shallow.2021) 
rms.pali<-rms.pali.2017
rms.maui7<-rms.maui7.2017

###

rms.maui6.env<-rbind(rms.maui6.2017.env, rms.maui6.2018.env, rms.maui6.2019.env, rms.maui6.2020.env, rms.maui6.2021.env)
rms.olowalu.deep.env<-rbind(rms.olowalu.deep.2017.env, rms.olowalu.deep.2019.env, rms.olowalu.deep.2020.env, rms.olowalu.deep.2021.env)
rms.olowalu.shallow.env<- rbind(rms.olowalu.shallow.2017.env,rms.olowalu.shallow.2018.env, rms.olowalu.shallow.2019.env,rms.olowalu.shallow.2020.env, rms.olowalu.shallow.2021.env)
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
  season.start.diff<-ifelse(leap_year(rms.hw$Date[n]),366-as.numeric(format(as.Date(paste(rms.hw$season[n],'-12-01', sep='')),'%j')),365-as.numeric(format(as.Date(paste(rms.hw$season[n],'-12-01', sep='')),'%j')))
  
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
  season.start.diff<-ifelse(leap_year(rms.env$Date[n]),366-as.numeric(format(as.Date(paste(rms.env$season[n],'-12-01', sep='')),'%j')),365-as.numeric(format(as.Date(paste(rms.env$season[n],'-12-01', sep='')),'%j')))
  
  if (rms.env$julian[n]>=season.start){
    rms.env$day.season[n]<-rms.env$julian[n]-season.start
  }else{
    rms.env$day.season[n]<-rms.env$julian[n]+season.start.diff
  }
  
  n<-n+1
}	

#correct year for December
rms.env$season[which(rms.env$julian>334)]<-rms.env$season[which(rms.env$julian>334)]+1

#subset to whale season= December 1st to April 30
rms.hw<-subset(rms.hw, month(Date)==12 | (month(Date)>=1 & month(Date)<=4))
rms.env<-subset(rms.env, month(Date)==12 | (month(Date)>=1 & month(Date)<=4))


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

#write.csv(rms.hw.summary, 'Results/Maui_RMSSPL_summary_hw.csv', row.names=F)
#write.csv(rms.hw.summary.month, 'Results/Maui_RMSSPL_summary_month_hw.csv', row.names=F)

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
#gam.rms.full.hw<-gam(SPL ~ s(day.season) + te(day.season, Hour, by=site) + s(Hour, bs='cc',by=site) + site + as.factor(season), data=rms.hw, method='REML', select=T)

#check if model is appropriate
gam.check(gam.rms.hw, old.style=T)

#get signifcances
anova(gam.rms.hw)

# Non-humpback whale environmental band

gam.rms.env<-gam(SPL ~ s(day.season, k=14) + s(Hour, bs='cc', by=site) + site + as.factor(season), data=rms.env, method='REML', select=T)
#gam.rms.env.full<-gam(SPL_1 ~ s(day.season) + te(day.season, Hour, by=site, k=20) + s(Hour, bs='cc',by=site) + site + as.factor(season), data=rms.env, method='REML', select=T)

#check if model is appropriate
gam.check(gam.rms.env, old.style=T)

#get signifcances
anova(gam.rms.env)

#---

#Create effect plots. Because plot.gam() doesn't allow for a lot of adjustments of the look of the plots,
#we're going to create them from scratch.

# (a) - pull out all smoothers 
smoother.dat.hw <- plot(gam.rms.hw, pages=1, scale=0)  # plot.gam returns a list of n elements, one per plot
smoother.dat.env <- plot(gam.rms.env, pages=1, scale=0)  # plot.gam returns a list of n elements, one per plot

#get x, fit, and se for each plot and save to DF
get_GAMsmoothers<-function(smoothers,n){
  sm.coef<-smoothers[[n]]
  sm.fit<-as.data.frame(sm.coef$fit)
  names(sm.fit)<-'fit'
  sm.fit$x<-as.vector(sm.coef$x)
  sm.fit$se<-as.vector(sm.coef$se)
  return(sm.fit)
}

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

para.dat.hw <- as.data.frame(termplot(gam.rms.hw, se = TRUE, plot = FALSE))
para.dat.env <- as.data.frame(termplot(gam.rms.env, se = TRUE, plot = FALSE))


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

#combine subplots into one plot and save to file

gam.plots.hw<-(gam.plot.seasonal.hw | gam.plot.diel.maui6.hw | gam.plot.diel.maui7.hw | gam.plot.diel.olo_d.hw) /
                  (gam.plot.diel.olo_s.hw | gam.plot.diel.pali.hw | gam.paraplot.sites.hw | gam.paraplot.seasons.hw)

gam.plots.env<-(gam.plot.seasonal.env | gam.plot.diel.maui6.env | gam.plot.diel.maui7.env | gam.plot.diel.olo_d.env) /
  (gam.plot.diel.olo_s.env | gam.plot.diel.pali.env | gam.paraplot.sites.env | gam.paraplot.seasons.env)


#ggsave(gam.plots.hw, filename='Plots/Maui_RMSSPL_GAM_hw.png', width=20, height=7.5, units='in',  bg = "white")    
#ggsave(gam.plots.env, filename='Plots/Maui_RMSSPL_GAM_env.png', width=20, height=7.5, units='in',  bg = "white")    


################################################################


### 4) POOL DATA FOR VISUALIZATIONS/POLAR PLOTS ----

#All hourly median RMS SPL values are normalized by the respective mean per season and site per frequency band
#Then, for each site and frequency band, the median per hour is calculated
#Medians for the humpback whale band are scaled between 0-1 for each site
#Medians for the non-humpback whale band are scaled between 0-1 across sites to illustrate the differing magnitudes of the patterns

#normalize hourly median SPL per site and year
rms.hw<-as.data.frame(rms.hw %>%
                        group_by(season,site) %>%
                        mutate(SPL.norm = normalize_fun(SPL)
                        ))
  
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

### Diel polar plots ---

# humpback whale band

rms.hw.h$Hour_label<-as.POSIXct(as.character(rms.hw.h$Hour), format = "%H:%M:%S", tz='UTC') #force UTC for easier plotting

diel.plots<-ggplot(rms.hw.h, aes(x = Hour_label, y=SPL.norm.median.scaled, fill=SPL.norm.median.scaled)) +
  polarplots_params(hw=T)

#ggsave(diel.plots, filename="Plots/Maui_RMSSPL_diel_hw.png", width=17.5, height=7.5, units='in', bg = "white")    

# Non-humpback whale environmental band

rms.env.h$Hour_label<-as.POSIXct(as.character(rms.env.h$Hour), format = "%H:%M:%S", tz='UTC') #force UTC for easier plotting

diel.env.plots<-ggplot(rms.env.h, aes(x = Hour_label, y=SPL.norm.median.scaled, fill=SPL.norm.median.scaled)) +
  polarplots_params(env=T)

#ggsave(diel.env.plots, filename="Plots/Maui_RMSSPL_diel_env.png", width=15, height=5, units='in', bg = "white")    


################################################################


#### 5) OCTAVE BANDS TIME-SERIES PLOTS FOR VISUALIZATION ----

#to illustrate seasonal and diel patterns for different frequency bands
#using data from one season of one offshore (Maui6) and one nearshore (Olowalu deep) site as representative examples
#daily averages are further smoothed with a running average of n=10

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
lab.y.day<-bquote("Median RMS SPL"~'[dB re 1'~plain("\u03BC")*'Pa], running average n=10')
lab.y.h<-bquote("Median RMS SPL"~'[dB re 1'~plain("\u03BC")*'Pa]')

#create color palette for different frequency bands, use modified viridis color palette
pkgTest('viridis')
cols<-viridis(6, direction=-1)
cols[1]<-"#FCA50AFF" #remove yellow because it's not very visible
names(cols) <- as.character(names(octaves.maui6.day[-1]))

#---  

#plot the seasonal time-series

octaves.maui6.plot.day<-ggplot(data=octaves.maui6.day.long, aes(x=Date, y=value, shape=variable)) +
  octaves_plots_params(seasonal=T, site='Maui6')
  
octaves.olowalu.plot.day<-ggplot(data=octaves.olowalu.day.long, aes(x=Date, y=value, shape=variable)) +
  octaves_plots_params(seasonal=T, site='Olowalu (deep)') +
  theme(plot.title = element_text(hjust=1))

#---  

#plot diel averages

octaves.maui6.plot.h<-ggplot(data=octaves.maui6.h.long, aes(x=Hour, y=value, shape=variable)) + 
  octaves_plots_params(diel=T, site='Maui6')

octaves.olowalu.plot.h<-ggplot(data=octaves.olowalu.h.long, aes(x=Hour, y=value, shape=variable)) + 
  octaves_plots_params(diel=T, site='Olowalu (deep)') +
  theme(plot.title = element_text(hjust=1))

#---

#arrange subplots into single plot and save to file

Maui.Octaves<-ggarrange(octaves.maui6.plot.day, octaves.olowalu.plot.day, octaves.maui6.plot.h,  octaves.olowalu.plot.h, ncol = 2, nrow = 2, common.legend=T)

#ggsave(Maui.Octaves, filename='Plots/Maui_RMSSPL_Octaves_Comparison.png', width=20, height=12, units='in',  bg = "white")    


################################################################

#### ----

#End of File

# -------

################################################################
