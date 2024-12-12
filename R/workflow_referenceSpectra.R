## LLHARP Workflow for creating reference spectra
#
# Script to process clicks manually selected in Pamguard into a reference 
# spectrum that can be used for assessing clicks to species level (ideally!)
#
# Process Overview
# 1

# ------ install/update necessary packages --------------------------------

# make sure you have Rtools installed
# if(!require('devtools')) install.packages('devtools')

# install latest ver from github (won't do anything if up to date!)
# devtools::install_github('taikisan21/PAMpal')


# # Installation of CRAN versions - only have to run once
# install.packages('banter')
# install.packages('PAMpal')

library(here)
library(PAMpal)
# library(banter)
# library(rfPermute)
library(tidyverse)
library(RSQLite)


# --- LL061 - FIRST EVENT -------------------------------------------------

# ------ USER DEFINED INPUTS ----------------------------------------------

tripStr = 'LL061'
recStr = 'FR45_DL84'
# Pamguard version
pgVer = '20207b'


# define transfer function
# calFile = 0         # for no calibration
# for LL050 onward...use this one
calFile = '//picqueenfish/psd2/crp/LLHARP/calibration/transfer_functions/RevN_PAMpal_invSensit.csv'

# define paths
# path_analysis should point to the folder above the 'pamguard' folder (that 
#               contains the 'databases' and 'binaries' subdirectories), must
#               contain a 'manualPicks' folder of event .csvs, a 'banterModels' 
#               folder with models for prediction, and an 'eventData folder
#               for outputs
path_analysis = '//piccrp4nas/psuedorcanassidens/analysis_phase3/'
# path_crptools = paste0(dirname(here()), '/crptools/') # assuming llamp and crptools are both in the users github directory
path_save = 'C:/Users/Selene.Fregosi/Documents/GitHub/llamp/files/reference_spectra'

# below objects are built from the above user defined inputs
path_pg = paste0(path_analysis, 'pamguard')

pgVerPrfx = paste0('pam', pgVer)
trStr = paste0(tripStr, '_', recStr)


# ------ PAMpal steps -----------------------------------------------------

# define parameters (database, binaries, default settings)
paramFile = file.path(path_analysis, 'eventData', tripStr, 
                      paste0(pgVerPrfx, '_', trStr, '_params.rda'))
if (!file.exists(paramFile)){
  source(here('R', 'functions', 'defineParams_phase3.R'))
  fkwPps = defineParams_phase3(pgVerPrfx, trStr, calFile, path_pg)
  # calibration prompts = 2 - dB re uPa/Counts, 16 - bit rate, 1 - apply
  save(fkwPps, file = paramFile)
}else if(file.exists(paramFile)){
  load(paramFile)
}


# Extracting data to build the AcousticStudy
dets <- processPgDetections(fkwPps, mode = 'db', id = paste0(trStr), 
                            progress = TRUE)
# event-based detection output
detFile = file.path(path_save, 'creatingNewReferenceSpectra_20230924',
                    paste0(pgVerPrfx, '_', trStr, '_dets_refSpecEvent.rda'))
save(dets, file = detFile)
# dets_og = dets # for testing

# set species
dets <- setSpecies(dets, method='pamguard')
# dets = setSpecies(dets, method = 'manual', value = 'Pc')
# # filter out only our species of interest
# dets <- filter(dets, species %in% c('Pc'))




# extract click data
ClickData <- getClickData(dets)
# some clicks are duplicated if they are classified under multiple Click Detectors
Clicks <- ClickData[!duplicated(ClickData$UID),]
# not sure why ClickData is bigger than just the event clicks (2112 vs 1743)

# extract clicks and events from database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, file.path(path_analysis, 'pamguard', 'databases', 
                                    tripStr, paste0(pgVerPrfx, '_llharp_banter_', 
                                                    trStr, '.sqlite3')))
Events <- dbReadTable(conn, "Click_Detector_OfflineEvents")      
EventClicks <- dbReadTable(conn, "Click_Detector_OfflineClicks")    
dbDisconnect(conn)

# EventClicks also still contains some duplicate clicks
# clear them
EventClicks = EventClicks[!duplicated(EventClicks$UID),]
# filter out specific events and clicks of our species of interest
# Events$eventType <- trimws(Events$eventType)
# Events <- subset(Events, eventType=="BWC" | eventType=="IP" | eventType=="MD" | eventType=="ZC" | eventType=="BW")
# EventClicks <- subset(EventClicks, parentID %in% c(Events$Id))

Spec_LL061 <- calculateAverageSpectra(dets, evNum = 1, wl = 256, channel = 1, 
                                      sr = 200000, norm = FALSE, plot = TRUE, 
                                      noise = FALSE)

# event-based detection output
clickFile = file.path(path_save, 'creatingNewReferenceSpectra_20230924', 
                      paste0(pgVerPrfx, '_', trStr, '_refSpecEventClicksOnly.rda'))
save(Spec_LL061, file = clickFile)




# --- LL067 - SECOND EVENT -------------------------------------------------

# ------ USER DEFINED INPUTS ----------------------------------------------

tripStr = 'LL067'
recStr = 'FR64_DL88'
# Pamguard version
pgVer = '20207b'


# define transfer function
# calFile = 0         # for no calibration
# for LL050 onward...use this one
calFile = '//picqueenfish/psd2/crp/LLHARP/calibration/transfer_functions/RevN_PAMpal_invSensit.csv'

# define paths
# path_analysis should point to the folder above the 'pamguard' folder (that 
#               contains the 'databases' and 'binaries' subdirectories), must
#               contain a 'manualPicks' folder of event .csvs, a 'banterModels' 
#               folder with models for prediction, and an 'eventData folder
#               for outputs
path_analysis = '//piccrp4nas/psuedorcanassidens/analysis_phase3/'
# path_crptools = paste0(dirname(here()), '/crptools/') # assuming llamp and crptools are both in the users github directory
path_save = 'C:/Users/Selene.Fregosi/Documents/GitHub/llamp/files/reference_spectra'

# below objects are built from the above user defined inputs
path_pg = paste0(path_analysis, 'pamguard')

pgVerPrfx = paste0('pam', pgVer)
trStr = paste0(tripStr, '_', recStr)


# ------ PAMpal steps -----------------------------------------------------

# define parameters (database, binaries, default settings)
paramFile = file.path(path_analysis, 'eventData', tripStr, 
                      paste0(pgVerPrfx, '_', trStr, '_params.rda'))
if (!file.exists(paramFile)){
  source(here('R', 'functions', 'defineParams_phase3.R'))
  fkwPps = defineParams_phase3(pgVerPrfx, trStr, calFile, path_pg)
  # calibration prompts = 2 - dB re uPa/Counts, 16 - bit rate, 1 - apply
  save(fkwPps, file = paramFile)
}else if(file.exists(paramFile)){
  load(paramFile)
}


# # point to Encounter Times file
# # phase 3 naming scheme
# encTimesFile = paste0(path_analysis, 'manualPicks/', tripStr, '/', trStr,
#                       '_encounterTimes.csv')

# # process the pamguard detections into an acoustics study object
# # (slow so only run if doesn't exist already)
# detFile = file.path(path_analysis, 'eventData', tripStr, 
#                     paste0(pgVerPrfx, '_', trStr, '_dets.rda'))
# if (!file.exists(detFile)){
#   dets = processPgDetections(fkwPps, mode = 'time', id = paste0(trStr),
#                              grouping = encTimesFile, 
#                              format = '%m/%d/%Y %H:%M')
#   # format = '%m/%d/%Y %H:%M:%S')
#   # the time format may need to be modified depending on how the csv was made
#   save(dets, file = detFile)
# }else if(file.exists(detFile)){
#   load(detFile)
# }


# Extracting data to build the AcousticStudy
dets <- processPgDetections(fkwPps, mode = 'db', id = paste0(trStr), 
                            progress = TRUE)
# event-based detection output
detFile = file.path(path_save, 'creatingNewReferenceSpectra_20230924', 
                    paste0(pgVerPrfx, '_', trStr, '_dets_refSpecEvent.rda'))
save(dets, file = detFile)
# dets_og = dets

# set species
dets <- setSpecies(dets, method='pamguard')
# dets = setSpecies(dets, method = 'manual', value = 'Pc')

# filter out only our species of interest
dets <- filter(dets, species %in% c('Pc'))




# extract click data
ClickData <- getClickData(dets)
# some clicks are duplicated if they are classified under multiple Click Detectors
Clicks <- ClickData[!duplicated(ClickData$UID),]
# not sure why ClickData is bigger than just the event clicks (2112 vs 1743)

# extract clicks and events from database
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, file.path(path_analysis, 'pamguard', 'databases', 
                                    tripStr, paste0(pgVerPrfx, '_llharp_banter_', 
                                                    trStr, '.sqlite3')))
Events <- dbReadTable(conn, "Click_Detector_OfflineEvents")      
EventClicks <- dbReadTable(conn, "Click_Detector_OfflineClicks")         
dbDisconnect(conn)

# EventClicks also still contains some duplicate clicks
# clear them
EventClicks = EventClicks[!duplicated(EventClicks$UID),]
# filter out specific events and clicks of our species of interest
# Events$eventType <- trimws(Events$eventType)
# Events <- subset(Events, eventType=="BWC" | eventType=="IP" | eventType=="MD" | eventType=="ZC" | eventType=="BW")
# EventClicks <- subset(EventClicks, parentID %in% c(Events$Id))

Spec_LL067 <- calculateAverageSpectra(dets, evNum = 1, wl = 256, channel = 1, 
                                      sr = 200000, norm = FALSE, plot = TRUE, 
                                      noise = FALSE)

# event-based detection output
clickFile = file.path(path_save, 'creatingNewReferenceSpectra_20230924', 
                      paste0(pgVerPrfx, '_', trStr, '_refSpecEventClicksOnly.rda'))
save(Spec_LL067, file = clickFile)


# --- COMBINE THEM --------------------------------------------------------

# reload both
load(file.path(path_save, 'creatingNewReferenceSpectra_20230924', 
               'PAM20207b_LL061_FR45_DL84_refSpecEventClicksOnly.rda'))
load(file.path(path_save, 'creatingNewReferenceSpectra_20230924', 
               'PAM20207b_LL067_FR64_DL88_refSpecEventClicksOnly.rda'))


# ------ Plot LL061 -------------------------------------------------------

# calc average (not normalized) for all clicks
# this is what PAMpal is using...which makes sense - see function 'doLogAvg
pres61_20log = 10^(Spec_LL061$allSpec/20) # back convert to pressure
presAvg61_20log = rowMeans(pres61_20log) # average pressures
dBAvg61_20log = 20*log10(presAvg61_20log) # convert pressures to dB

#WRONG
pow61_10log = 10^(Spec_LL061$allSpec/10) # back convert to power
powAvg61_10log = rowMeans(pow61_10log) # average powers
dBAvg61_10log = 10*log10(powAvg61_10log) # convert powers to dB

simpAvg61 = rowMeans(Spec_LL061$allSpec) # avg of dB without back conv (WRONG!)

# normAdj_20log = Spec_LL061$avgSpec[1] - dBAvg61_20log[1] # get offset from normalized
# normAdj_10log = Spec_LL061$avgSpec[1] - dBAvg61_10log[1] # get offset from normalized
# try plotting each individual click on one plot with avg overlaid
plot(1, type = 'n', xlim = c(0, 100), ylim = c(50, 190), 
     xlab = 'Frequency (kHz)', ylab = 'Magnitude (dB)', 
     main = 'Average Spectrum - LL061')
# add grid lines for frequency at 10 kHz intervals
for (iline in ((0:15)*10)) {lines(c(iline,iline),c(40, 200),col="darkgray")}
for (f in 1:ncol(Spec_LL061$allSpec)){
  lines(Spec_LL061$freq/1000, Spec_LL061$allSpec[,f], 
        col = rgb(red = 0.8, green = 0.8, blue = 0.8, alpha = 0.2))
}

lines(Spec_LL061$freq/1000, dBAvg61_20log, col = 'orange', lwd = 2) # CORRECT
lines(Spec_LL061$freq/1000, dBAvg61_10log, col = 'blue', lwd = 2) #WRONG
# lines(Spec_LL061$freq/1000, Spec_LL061$avgSpec-normAdj_10log, col = 'black', lty = 2)
lines(Spec_LL061$freq/1000, Spec_LL061$avgSpec, col = 'black', lty = 2) # from PAMpal
lines(Spec_LL061$freq/1000, simpAvg61, col = 'black') # WRONG
legend('bottomleft', legend = c('10*log10', '20*log10', 'PAMpal output', 'dB avg'), 
       col = c('blue', 'orange', 'black', 'black'), lty = c(1, 1, 2, 1), 
       lwd = c(2, 2, 1, 1))


# ------ Plot LL067 -------------------------------------------------------

# calc average (not normalized) for all clicks
# this is what PAMpal is using...which makes sense - see function 'doLogAvg
pres67_20log = 10^(Spec_LL067$allSpec/20) # back convert to pressure
presAvg67_20log = rowMeans(pres67_20log) # average pressures
dBAvg67_20log = 20*log10(presAvg67_20log) # convert pressures to dB

#WRONG
pow67_10log = 10^(Spec_LL067$allSpec/10) # back convert to power
powAvg67_10log = rowMeans(pow67_10log) # average powers
dBAvg67_10log = 10*log10(powAvg67_10log) # convert powers to dB

simpAvg67 = rowMeans(Spec_LL067$allSpec) # avg of dB without back conv (WRONG!)

# try plotting each individual click on one plot with avg overlaid
plot(1, type = 'n', xlim = c(0, 100), ylim = c(50, 190), 
     xlab = 'Frequency (kHz)', ylab = 'Magnitude (dB)', 
     main = 'Average Spectrum - LL067')
# add grid lines for frequency at 10 kHz intervals
for (iline in ((0:15)*10)) {lines(c(iline,iline),c(40, 200),col="darkgray")}
for (f in 1:ncol(Spec_LL067$allSpec)){
  lines(Spec_LL067$freq/1000, Spec_LL067$allSpec[,f], 
        col = rgb(red = 0.8, green = 0.8, blue = 0.8, alpha = 0.2))
}

lines(Spec_LL067$freq/1000, dBAvg67_20log, col = 'orange', lwd = 2) # CORRECT
lines(Spec_LL067$freq/1000, dBAvg67_10log, col = 'blue', lwd = 2) #WRONG
# lines(Spec_LL067$freq/1000, Spec_LL067$avgSpec-normAdj_10log, col = 'black', lty = 2)
lines(Spec_LL067$freq/1000, Spec_LL067$avgSpec, col = 'black', lty = 2) # from PAMpal
lines(Spec_LL067$freq/1000, simpAvg67, col = 'black') # WRONG
legend('bottomleft', legend = c('10*log10', '20*log10', 'PAMpal output', 'dB avg'), 
       col = c('blue', 'orange', 'black', 'black'), lty = c(1, 1, 2, 1), 
       lwd = c(2, 2, 1, 1))


# ------ Average the two events -------------------------------------------
# use 20 log because that is for pressure!
ncol(pres61_20log) #1708
ncol(pres67_20log) #1668
pres = cbind(pres61_20log, pres67_20log)
ncol(pres) #3376

presAvg = rowMeans(pres) #average pressures
dBAvg = 20*log10(presAvg) # convert to dB

# ------ Plot both together -----------------------------------------------
plot(1, type = 'n', xlim = c(0, 100), ylim = c(50, 190), 
     xlab = 'Frequency (kHz)', ylab = 'Magnitude (dB)', 
     main = 'Average Spectrum - LL061 and LL067')
# add grid lines for frequency at 10 kHz intervals
for (iline in ((0:15)*10)) {lines(c(iline,iline),c(40, 200),col="darkgray")}
for (f in 1:ncol(Spec_LL061$allSpec)){
  lines(Spec_LL061$freq/1000, Spec_LL061$allSpec[,f], 
        col = rgb(red = 0.8, green = 0.8, blue = 0.8, alpha = 0.2))
}
for (f in 1:ncol(Spec_LL067$allSpec)){
  lines(Spec_LL067$freq/1000, Spec_LL067$allSpec[,f], 
        col = rgb(red = 0.8, green = 0.8, blue = 0.8, alpha = 0.2))
}

lines(Spec_LL061$freq/1000, dBAvg61_20log, col = 'navyblue', lwd = 1) # CORRECT
lines(Spec_LL067$freq/1000, dBAvg67_20log, col = 'darkred', lwd = 1) # CORRECT
lines(Spec_LL067$freq/1000, dBAvg, col = 'mediumpurple4', lwd = 2) # CORRECT

# and compare to old one.
path_to_refSpec = '//piccrp4nas/psuedorcanassidens/analysis_phase3/TritonLogs/GitHub/llamp/files/reference_spectra'
refSpec = read.csv(file.path(path_to_refSpec, paste0('meanSpecPc.csv')))

# normalize to line up with our spec
rsAdj = max(dBAvg) - max(refSpec$dB)
lines(refSpec$frq, refSpec$dB + rsAdj, col = 'grey40', lwd = 1.5, lty = 2) 

# replot newest one so on top
lines(Spec_LL067$freq/1000, dBAvg, col = 'mediumpurple4', lwd = 2)
legend('bottomleft', legend = c('LL061', 'LL067', 'LL combined', 'old'), 
       col = c('navyblue', 'darkred', 'mediumpurple4', 'grey40'), 
       lty = c(1, 1, 1, 2), lwd = c(1, 1, 2, 1.5))


# --- SAVE IT -------------------------------------------------------------
# coerce into data.frame
refSpecNew = data.frame(dB = dBAvg, frq = Spec_LL067$freq/1000) 
write.csv(refSpecNew, file = file.path(path_save, 'refSpecPc_LLHARP.csv'), 
          row.names = FALSE)

