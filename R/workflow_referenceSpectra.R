## Workflow for creating reference spectra for event reports
#
# Script to process groups of clicks that have been manually selected in 
# Pamguard into a reference spectrum that can be used for assessing clicks 
# within an event summary report
#
# Process Overview
# 1

# ------ install/update necessary packages --------------------------------

# make sure you have Rtools installed
# if(!require('devtools')) install.packages('devtools')

# install latest ver from github (won't do anything if up to date!)
# devtools::install_github('taikisan21/PAMpal')

# # Installation of CRAN versions - only have to run once
# install.packages('PAMpal')

library(here)
library(PAMpal)
library(tidyverse)
library(RSQLite)


# ------ USER DEFINED INPUTS ----------------------------------------------

# calFile = '//picqueenfish/psd2/crp/LLHARP/calibration/transfer_functions/RevN_PAMpal_invSensit.csv'
sampleRate <- 288000 #576000 or 288000 for MACS2018  #384000 for HICEAS 2023

# define paths
# path_database <- '//piccrpnas/CRP2/McCullough/KOGIA_DASBR_POST/2023_HICEAS/Database'
# databaseFile <- file.path(path_database, 'HICEAS_2023_DS12a_Kogia_PAM20210a.sqlite3')
# path_binaries <- '//piccrpnas/CRP2/McCullough/KOGIA_DASBR_POST/2023_HICEAS/Binaries/DS12a'

path_database <- '//piccrpnas/CRP/MACS_2018_DASBR/Post_Processing/MACS_DASBR_PostProcessing_JLK/Databases'
databaseFile <- file.path(path_database, 'DS3_ST-3_PAM20014c.sqlite3')
path_binaries <- '//piccrpnas/CRP/MACS_2018_DASBR/Post_Processing/MACS_DASBR_PostProcessing_JLK/Binaries/DS3_ST-3'

# specify path to save to
path_save <- 'C:/Users/Selene.Fregosi/Documents/GitHub/SpermWhale_DASBR/refSpec'
# specify some output filenames
paramFile <- file.path(path_save, 'Pm_refSpec_MACS2018_DS3_params.rda')
detFile <- file.path(path_save, 'Pm_refSpec_MACS2018_DS3_dets.rda')
clickFile <- file.path(path_save, 'Pm_refSpec_MACS2018_DS3_eventClicks.rda')
refSpecFile <- file.path(path_save, 'Pm_refSpec_MACS2018_DS3.csv')


# ------ PAMpal steps -----------------------------------------------------

#Create a settings object (loading databases and binaries) 
#below example filter setting is for Kogia data
# pps <- PAMpalSettings(db = databaseFile, binaries = path_binaries, 
#                       sr_hz = 'auto', winLen_sec = 0.0025, 
#                       filterfrom_khz = 80, filterto_khz = NULL)
# typical settings for other odonts are:
pps <- PAMpalSettings(db = databaseFile, binaries = path_binaries,
                      sr_hz = 'auto', winLen_sec = 0.0025,
                      filterfrom_khz = 2, filterto_khz = NULL)
# save it because it can take a while
saveRDS(pps, file = paramFile)

# create the acoustic study based on the single event
dets <- processPgDetections(pps, mode = 'db', id = 'Id-1', 
                            progress = TRUE)
saveRDS(dets, file = detFile)


# set species
dets <- setSpecies(dets, method='pamguard')
# dets = setSpecies(dets, method = 'manual', value = 'Pc')
# # filter out only our species of interest
# dets <- filter(dets, species %in% c('Pc'))


# extract clicks and events from database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, databaseFile)
events <- dbReadTable(conn, "Click_Detector_OfflineEvents")      
eventClicks <- dbReadTable(conn, "Click_Detector_OfflineClicks")    
dbDisconnect(conn)

# clicks can be duplicated if they are classified more than once so remove them
eventClicks <- eventClicks[!duplicated(eventClicks$UID),]

avSpec <- calculateAverageSpectra(dets, evNum = 1, wl = 256, channel = 1, 
                                      sr = sampleRate, norm = FALSE, plot = TRUE, 
                                      noise = FALSE)

# save event click spectrum data
saveRDS(avSpec, file = clickFile)


# ------ Plot to check ----------------------------------------------------

# calc average (not normalized) for all clicks
# this is what PAMpal is using...which makes sense - see function 'doLogAvg'
pres <- 10^(avSpec$allSpec/20) # back convert to pressure
presAvg <- rowMeans(pres) # average pressures
dBAvg <- 20*log10(presAvg) # convert pressures to dB

# try plotting each individual click on one plot with avg overlaid
yMin <- round(min(avSpec$allSpec), -1) - 10
yMax <- round(max(avSpec$allSpec), -1) + 10
plot(1, type = 'n', xlim = c(0, sampleRate/2/1000), ylim = c(yMin, yMax), 
     xlab = 'Frequency (kHz)', ylab = 'Magnitude (dB)', 
     main = 'Average Spectrum')
# add grid lines for frequency at 10 kHz intervals
for (iline in ((0:floor(sampleRate/20000))*10)) {
  lines(c(iline,iline), c(yMin-10, yMax+10),col="darkgray")}
for (f in 1:ncol(avSpec$allSpec)){
  lines(avSpec$freq/1000, avSpec$allSpec[,f], 
        col = rgb(red = 0.8, green = 0.8, blue = 0.8, alpha = 0.2))
}

lines(avSpec$freq/1000, dBAvg_20log, col = 'orange', lwd = 2) # CORRECT
# lines(Spec_LL061$freq/1000, Spec_LL061$avgSpec-normAdj_10log, col = 'black', lty = 2)
lines(avSpec$freq/1000, avSpec$avgSpec, col = 'black', lty = 2) # from PAMpal
legend('bottomleft', legend = c('individual clicks', '20*log10', 'PAMpal output'), 
       col = c('grey', 'orange', 'black'), lty = c(1, 1, 2), 
       lwd = c(0.5, 2, 1))


# ------ Final clean up and save to CSV -----------------------------------

# coerce into a cleaner data frame and save
refSpec <- data.frame(dB = dBAvg, frq = avSpec$freq/1000) 
write.csv(refSpec, file = refSpecFile, row.names = FALSE)


# ------ WORKING WITH MULTIPLE EVENTS -------------------------------------

# compile pressure measures for each click
pres <- cbind(pres_20log_ev1, pres_20log_ev2)
# average the pressures and convert to dB
presAvg <- rowMeans(pres) #average pressures
dBAvg <- 20*log10(presAvg) # convert to dB

# plot them all together
plot(1, type = 'n', xlim = c(0, sampleRate/2/1000), ylim = c(yMin, yMax), 
     xlab = 'Frequency (kHz)', ylab = 'Magnitude (dB)', 
     main = 'Average Spectrum - Multiple Events')
# add grid lines for frequency at 10 kHz intervals
for (iline in ((0:floor(sampleRate/20000))*10)) {
  lines(c(iline,iline), c(yMin-10, yMax+10),col="darkgray")}
# add the individual clicks of each event
for (f in 1:ncol(avSpec_ev1$allSpec)){
  lines(avSpec_ev1$freq/1000,avSpec_ev1$allSpec[,f], 
        col = rgb(red = 0.8, green = 0.8, blue = 0.8, alpha = 0.2))
}
for (f in 1:ncol(avSpec_ev2$allSpec)){
  lines(avSpec_ev2$freq/1000, avSpec_ev2$allSpec[,f], 
        col = rgb(red = 0.8, green = 0.8, blue = 0.8, alpha = 0.2))
}
# overlay the individual averages and combined average
lines(avSpec_ev1$freq/1000, dBAvg_20log_ev1, col = 'navyblue', lwd = 1) 
lines(avSpec_ev2$freq/1000, dBAvg_20log_ev2, col = 'darkred', lwd = 1) 
lines(avSpec_ev1$freq/1000, dBAvg, col = 'mediumpurple4', lwd = 2) 

# coerce into a cleaner data frame and save
refSpec <- data.frame(dB = dBAvg, frq = avSpec$freq/1000) 
write.csv(refSpec, file = refSpecFile, row.names = FALSE)


