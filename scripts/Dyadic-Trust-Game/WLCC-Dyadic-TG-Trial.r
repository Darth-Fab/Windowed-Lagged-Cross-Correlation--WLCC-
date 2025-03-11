## 06-04-2024 - Fabiola Diana
# Load libraries & data files -----

library(lmerTest)
library(ggplot2)
library(sjPlot)
library(pbkrtest)
library(effects)
library(piecewiseSEM)
library(effectsize)
library(dplyr)
library(magrittr)
library(viridis) #NOTE: added for colour palette
options(scipen = 999) # forces R to convert scientific notations to decimals

# set working directory
setwd('C:/Users/fabio/Documents/GitHub/Windowed-Lagged-Cross-Correlation--WLCC-/scripts/Dyadic-Trust-Game/')

d_raw <- read.csv(file= 'data/Merged_SC_ECG_12Trial_Level_clean.csv')
#d_raw <- read.csv(file= 'data/WCC_DIA_SC_ECG-trial.csv', sep=",")
unique(d_raw$Condition)
unique(d_raw$trial)
unique(d_raw$Dyad)
colnames(d_raw)

### Create numeric trial 
# Numbers 1 - 6 = FaceBlocked Trials
# Numbers 7 - 12 = FaceToFace Trials
d_raw$Trial_numeric <- as.integer(d_raw$trial)

unique(d_raw$Trial_numeric)

######### LOAD FUNCTIONS ########

source("functions/wcc-peak.R")
source("functions/mean-se.R")
source("functions/vif.R")
source("functions/rsquared.R")


######### PERFORM WCC ANALYSIS #########
d_raw_all<-d_raw
num = seq(1,5,by=1)
num2 = seq(11,15, by=1)
num3= seq(12,15, by=1)

#Only selected dyad
d_raw <- d_raw[is.element(d_raw$Dyad, "13"),]%>%droplevels
# to run the analysis, we first make lists containing matrices with the raw HR and SCL data per face condition for each dyad


# Get unique dyads
dyads <- unique(d_raw$Dyad)
numdyad <- length(dyads)

# Initialize lists dynamically for all trials
nofaceHRlist <- vector("list", 6)
faceHRlist <- vector("list", 6)
nofaceSClist <- vector("list", 6)
faceSClist <- vector("list", 6)

for (t in 1:6) {
  nofaceHRlist[[t]] <- vector("list", numdyad)
  faceHRlist[[t]] <- vector("list", numdyad)
  nofaceSClist[[t]] <- vector("list", numdyad)
  faceSClist[[t]] <- vector("list", numdyad)
}

# Populate lists for all dyads
for (l in seq_along(dyads)) {
  dyad_data <- d_raw[d_raw$Dyad == dyads[l], ]
  
  for (t in 1:6) {
    nofaceHRlist[[t]][[l]] <- dyad_data[dyad_data$Trial_numeric == t, ]
    faceHRlist[[t]][[l]] <- dyad_data[dyad_data$Trial_numeric == t + 6, ]
    nofaceSClist[[t]][[l]] <- dyad_data[dyad_data$Trial_numeric == t, ]
    faceSClist[[t]][[l]] <- dyad_data[dyad_data$Trial_numeric == t + 6, ]
  }
}

# Skin conductance ----
numdyadWCC <- length(dyads)

# Define the number of trials
numTrials <- 6

# Initialize lists dynamically
FaceWCCList_SC <- vector("list", numTrials)
FacePeakList_SC <- vector("list", numTrials)
NoFaceWCCList_SC <- vector("list", numTrials)
NoFacePeakList_SC <- vector("list", numTrials)

for (t in 1:numTrials) {
  FaceWCCList_SC[[t]] <- vector("list", numdyadWCC)
  FacePeakList_SC[[t]] <- vector("list", numdyadWCC)
  NoFaceWCCList_SC[[t]] <- vector("list", numdyadWCC)
  NoFacePeakList_SC[[t]] <- vector("list", numdyadWCC)
}
# Number of samples (sampling rate of 20 Hz)
wSize <- 20 * 8   # 8 sec
tMax <- 20 * 4    # 4 sec
wInc <- 20 * 2    # 2 sec
tInc <- 2         # 100 ms
LoessSpan <- 0.25
L <- 25
SamplingRate <- 20  # Hz

# Color palette for plots
col1 <- colorRampPalette(viridis(100))

# Progress bar
pb <- txtProgressBar(min = 0, max = numdyadWCC, style = 3)

# Perform analysis for all trials and all dyads
for (t in 1:numTrials) {  
    for (i in 1:numdyadWCC) {  
        try(if (nrow(faceSClist[[t]][[i]]) != 0 & nrow(nofaceSClist[[t]][[i]]) != 0) {
            
            print(Sys.time())
            
            # WCC Analysis
            FaceWCCList_SC[[t]][[i]] <- WCC(faceSClist[[t]][[i]][,"SC_PPN1_SkinConductance"], 
                                            faceSClist[[t]][[i]][,"SC_PPN2_SkinConductance"])
            FaceWCCList_SC[[t]][[i]] <- na.omit(t(FaceWCCList_SC[[t]][[i]]))
            FaceWCCList_SC[[t]][[i]] <- t(FaceWCCList_SC[[t]][[i]]) 
            write.csv(FaceWCCList_SC[[t]][[i]], 
                      file = paste("SC_D", i, "_Trial-", t, "_Face_W", wSize, "_Incr", wInc, ".csv", sep = ""), 
                      row.names = FALSE)
            
            NoFaceWCCList_SC[[t]][[i]] <- WCC(nofaceSClist[[t]][[i]][,"SC_PPN1_SkinConductance"], 
                                              nofaceSClist[[t]][[i]][,"SC_PPN2_SkinConductance"])
            NoFaceWCCList_SC[[t]][[i]] <- na.omit(t(NoFaceWCCList_SC[[t]][[i]]))
            NoFaceWCCList_SC[[t]][[i]] <- t(NoFaceWCCList_SC[[t]][[i]])
            write.csv(NoFaceWCCList_SC[[t]][[i]], 
                      file = paste("SC_D", i, "_Trial-", t, "_NoFace_W", wSize, "_Incr", wInc, ".csv", sep = ""), 
                      row.names = FALSE)
            
            # Peak Picking
            FacePeakList_SC[[t]][[i]] <- peakpick(FaceWCCList_SC[[t]][[i]], pspan = LoessSpan, graphs = 1)
            write.csv(FacePeakList_SC[[t]][[i]], 
                      file = paste("D", i, "_FacePeakList_SC_", t, ".csv", sep = ""), 
                      row.names = FALSE)
            
            NoFacePeakList_SC[[t]][[i]] <- peakpick(NoFaceWCCList_SC[[t]][[i]], pspan = LoessSpan)
            write.csv(NoFacePeakList_SC[[t]][[i]], 
                      file = paste("D", i, "_NoFacePeakList_SC_", t,  ".csv", sep = ""), 
                      row.names = FALSE)

            # Update progress bar
            setTxtProgressBar(pb, i)
        })
    }
}

close(pb)

# Save results
for (t in 1:numTrials) {
  # Save the individual sublist for Face and NoFace conditions for each trial
  save(list = c("FacePeakList_SC"), file = paste("Trial-", t, "-SC_PeakList_Face_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))
  save(list = c("NoFacePeakList_SC"), file = paste("Trial-", t, "-SC_PeakList_NoFace_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))
}


# Heart Rate -------
# For HR, create lists for Face and NoFace conditions for all trials
FaceWCCList_HR <- vector("list", numdyadWCC)
FacePeakList_HR <- vector("list", numdyadWCC)
NoFaceWCCList_HR <- vector("list", numdyadWCC)
NoFacePeakList_HR <- vector("list", numdyadWCC)

for (t in 1:numTrials) {
  FaceWCCList_HR[[t]] <- vector("list", numdyadWCC)
  FacePeakList_HR[[t]] <- vector("list", numdyadWCC)
  NoFaceWCCList_HR[[t]] <- vector("list", numdyadWCC)
  NoFacePeakList_HR[[t]] <- vector("list", numdyadWCC)
}

# Number of samples (using the same as SC)
wSize <- 20 * 8   # 8 sec
tMax <- 20 * 4    # 4 sec
wInc <- 20 * 2    # 2 sec
tInc <- 2         # 100 ms
LoessSpan <- 0.25
L <- 25
SamplingRate <- 20  # Hz

# Set up the color palette
col1 <- colorRampPalette(viridis(100))

# Progress bar
pb <- txtProgressBar(min = 0, max = numdyadWCC, style = 3)

### Loop for all trials
for (t in 1:numTrials) {
  for (i in 1:numdyadWCC) {
    try({
      # Make sure there is data for both conditions
      if (nrow(faceHRlist[[t]][[i]]) != 0 & nrow(nofaceHRlist[[t]][[i]]) != 0) {
        
        print(Sys.time())
        
        # WCC Analysis for Face condition (Heart Rate)
        FaceWCCList_HR[[t]][[i]] <- WCC(faceHRlist[[t]][[i]][, "ECG_PPN1_HeartRate"], 
                                        faceHRlist[[t]][[i]][, "ECG_PPN2_HeartRate"])
        FaceWCCList_HR[[t]][[i]] <- na.omit(t(FaceWCCList_HR[[t]][[i]]))
        FaceWCCList_HR[[t]][[i]] <- t(FaceWCCList_HR[[t]][[i]])  # Adjust dimensions for WCC
        write.csv(FaceWCCList_HR[[t]][[i]], 
                  file = paste("HR_D", i, "_Trial-", t, "_Face_W", wSize, "_Incr", wInc, ".csv", sep = ""), 
                  row.names = FALSE)
        
        # WCC Analysis for NoFace condition (Heart Rate)
        NoFaceWCCList_HR[[t]][[i]] <- WCC(nofaceHRlist[[t]][[i]][, "ECG_PPN1_HeartRate"], 
                                          nofaceHRlist[[t]][[i]][, "ECG_PPN2_HeartRate"])
        NoFaceWCCList_HR[[t]][[i]] <- na.omit(t(NoFaceWCCList_HR[[t]][[i]]))
        NoFaceWCCList_HR[[t]][[i]] <- t(NoFaceWCCList_HR[[t]][[i]])  # Adjust dimensions for WCC
        write.csv(NoFaceWCCList_HR[[t]][[i]], 
                  file = paste("HR_D", i, "_Trial-", t, "_NoFace_W", wSize, "_Incr", wInc, ".csv", sep = ""), 
                  row.names = FALSE)
        
        # Peak picking for Face condition
        FacePeakList_HR[[t]][[i]] <- peakpick(FaceWCCList_HR[[t]][[i]], pspan = LoessSpan, graphs = 1)
        write.csv(FacePeakList_HR[[t]][[i]], 
                  file = paste("D", i, "_FacePeakList_HR_", t, ".csv", sep = ""), 
                  row.names = FALSE)
        
        # Peak picking for NoFace condition
        NoFacePeakList_HR[[t]][[i]] <- peakpick(NoFaceWCCList_HR[[t]][[i]], pspan = LoessSpan)
        write.csv(NoFacePeakList_HR[[t]][[i]], 
                  file = paste("D", i, "_NoFacePeakList_HR_", t, ".csv", sep = ""), 
                  row.names = FALSE)
        
        # Update progress bar
        setTxtProgressBar(pb, i)
      }
    })
  }
}

close(pb)

# Save results for HR
for (t in 1:numTrials) {
  # Save the individual sublist for Face and NoFace conditions for each trial (Heart Rate)
  save(list = c("FacePeakList_HR"), file = paste("Trial-", t, "-HR_PeakList_Face_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))
  save(list = c("NoFacePeakList_HR"), file = paste("Trial-", t, "-HR_PeakList_NoFace_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))
}


######## SUMMARY STATISTICS ########
summary_data <- data.frame()

# Define the number of trials and dyads (from your code context)
numTrials <- 6
numdyadWCC <- length(dyads)
SamplingRate = 20
tInc_HR = 20*0.5
tInc_SC = 20*0.5

# Loop over each dyad and trial
for (i in 1:numdyadWCC) {
  for (t in 1:numTrials) {
    # Construct the file paths for the relevant CSV files
    FaceHRFile <- paste("D", i, "_FacePeakList_HR_", t, ".csv", sep = "")
    NoFaceHRFile <- paste("D", i, "_NoFacePeakList_HR_", t, ".csv", sep = "")
    FaceSCFile <- paste("D", i, "_FacePeakList_SC_", t, ".csv", sep = "")
    NoFaceSCFile <- paste("D", i, "_NoFacePeakList_SC_", t, ".csv", sep = "")
    
    # Read the files (make sure they exist)
    if (file.exists(FaceHRFile) & file.exists(NoFaceHRFile) & 
        file.exists(FaceSCFile) & file.exists(NoFaceSCFile)) {
      
      # Load the data
      faceHR_data <- read.csv(FaceHRFile)
      nofaceHR_data <- read.csv(NoFaceHRFile)
      faceSC_data <- read.csv(FaceSCFile)
      nofaceSC_data <- read.csv(NoFaceSCFile)
      
      # Calculate the summary statistics for HR Face condition
      WCCmean_HR_Face <- mean(faceHR_data$maxValue, na.rm = TRUE)
      WCCsd_HR_Face <- sd(faceHR_data$maxValue, na.rm = TRUE)
      TLAmean_HR_Face <- mean(abs(faceHR_data$maxIndex / (2 * SamplingRate / tInc_HR)), na.rm = TRUE)
      TLAsd_HR_Face <- sd(abs(faceHR_data$maxIndex / (2 * SamplingRate / tInc_HR)), na.rm = TRUE)
      
      # Calculate the summary statistics for SC Face condition
      WCCmean_SC_Face <- mean(faceSC_data$maxValue, na.rm = TRUE)
      WCCsd_SC_Face <- sd(faceSC_data$maxValue, na.rm = TRUE)
      TLAmean_SC_Face <- mean(abs(faceSC_data$maxIndex / (2 * SamplingRate / tInc_HR)), na.rm = TRUE)
      TLAsd_SC_Face <- sd(abs(faceSC_data$maxIndex / (2 * SamplingRate / tInc_HR)), na.rm = TRUE)
      
      # Calculate the summary statistics for HR NoFace condition
      WCCmean_HR_NoFace <- mean(nofaceHR_data$maxValue, na.rm = TRUE)
      WCCsd_HR_NoFace <- sd(nofaceHR_data$maxValue, na.rm = TRUE)
      TLAmean_HR_NoFace <- mean(abs(nofaceHR_data$maxIndex / (2 * SamplingRate / tInc_HR)), na.rm = TRUE)
      TLAsd_HR_NoFace <- sd(abs(nofaceHR_data$maxIndex / (2 * SamplingRate / tInc_HR)), na.rm = TRUE)
      
      # Calculate the summary statistics for SC NoFace condition
      WCCmean_SC_NoFace <- mean(nofaceSC_data$maxValue, na.rm = TRUE)
      WCCsd_SC_NoFace <- sd(nofaceSC_data$maxValue, na.rm = TRUE)
      TLAmean_SC_NoFace <- mean(abs(nofaceSC_data$maxIndex / (2 * SamplingRate / tInc_HR)), na.rm = TRUE)
      TLAsd_SC_NoFace <- sd(abs(nofaceSC_data$maxIndex / (2 * SamplingRate / tInc_HR)), na.rm = TRUE)
      
      # Add the row to the summary data
      summary_data <- rbind(summary_data, data.frame(
        Dyad = i,
        Trial = t,
        WCCmean_HR_Face = WCCmean_HR_Face,
        WCCsd_HR_Face = WCCsd_HR_Face,
        TLAmean_HR_Face = TLAmean_HR_Face,
        TLAsd_HR_Face = TLAsd_HR_Face,
        WCCmean_SC_Face = WCCmean_SC_Face,
        WCCsd_SC_Face = WCCsd_SC_Face,
        TLAmean_SC_Face = TLAmean_SC_Face,
        TLAsd_SC_Face = TLAsd_SC_Face,
        WCCmean_HR_NoFace = WCCmean_HR_NoFace,
        WCCsd_HR_NoFace = WCCsd_HR_NoFace,
        TLAmean_HR_NoFace = TLAmean_HR_NoFace,
        TLAsd_HR_NoFace = TLAsd_HR_NoFace,
        WCCmean_SC_NoFace = WCCmean_SC_NoFace,
        WCCsd_SC_NoFace = WCCsd_SC_NoFace,
        TLAmean_SC_NoFace = TLAmean_SC_NoFace,
        TLAsd_SC_NoFace = TLAsd_SC_NoFace
      ))
    }
  }
}

# Write the final summary dataset to a CSV file
write.csv(summary_data, "summary_statistics.csv", row.names = FALSE)
