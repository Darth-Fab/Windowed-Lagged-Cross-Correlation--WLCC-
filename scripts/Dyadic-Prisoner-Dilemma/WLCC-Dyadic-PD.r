## 17-03-2025 - Fabiola Diana

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
library(tidyr)
library(viridis) #NOTE: added for colour palette
options(scipen = 999) # forces R to convert scientific notations to decimals

# set working directory
setwd('C:/Users/fabio/Documents/GitHub/Windowed-Lagged-Cross-Correlation--WLCC-/scripts/Dyadic-Prisoner-Dilemma/')

d_raw <- read.csv(file= 'data/ECG-SC-20Hz-WideFormat-AllEpochs.csv')
d_raw <- d_raw %>% drop_na()

######### LOAD FUNCTIONS ########

source("functions/wcc-peak.R")
source("functions/mean-se.R")
source("functions/vif.R")
source("functions/rsquared.R")

################ PERFORM WCC ANALYSIS ###############
d_raw_all<-d_raw
#d_raw <- d_raw[is.element(d_raw$Dyad, "13"),]%>%droplevels
d_raw$Dyad <- as.numeric(d_raw$Dyad)
d_raw <- d_raw[order(d_raw$Dyad), ]
d_raw$Dyad_new <- match(d_raw$Dyad, unique(d_raw$Dyad))
d_raw$Dyad_new <- as.numeric(d_raw$Dyad_new)

# Create a new dataframe with only the relevant columns
dyad_mapping <- unique(d_raw[, c("Dyad", "Dyad_new")])
names(dyad_mapping)[names(dyad_mapping) == "Dyad"] <- "Real_Dyad"
names(dyad_mapping)[names(dyad_mapping) == "Dyad_new"] <- "WCC_Dyad"

# Save it to a CSV file
write.csv(dyad_mapping, "DyadMapping-CondLevel-PrisonerDilemma.csv", row.names = FALSE)

dyads <- unique(d_raw$Dyad_new)
numdyad <- length(dyads)
# Create list for SC data per dyad
SClist <- vector("list", numdyad)
HRlist <- vector("list", numdyad)

SCmatrix <- matrix(NA, nrow = 160000, ncol = ncol(d_raw))
HRmatrix <- matrix(NA, nrow = 160000, ncol = ncol(d_raw))

# this is not working - take from d_raw instead of d (at least values in table then)

l = 1
for (i in 1:numdyad){
  HRmatrix <- d_raw[which(d_raw$Dyad_new == i),]
  SCmatrix <- d_raw[which(d_raw$Dyad_new == i),]
  HRlist[[l]] <- HRmatrix
  SClist[[l]] <- SCmatrix
  l = l + 1
}

# Skin Conductance ------ 
numdyadWCC <- length(dyads)

WCCList_SC<-vector("list", numdyadWCC)
PeakList_SC<-vector("list", numdyadWCC)

# number of samples (we have a sample rate of 20 Hz, so 20 samples = 1 sec; 1 sample = 50 ms)
wSize <- 20*8 # 8 sec
tMax <- 20*4 # 4 sec
wInc <- 20*2 # 2 sec
tInc <- 2 # 100 ms
LoessSpan <- 0.25
L <- 25
SamplingRate <- 20 # Hz

col1<-colorRampPalette(viridis(100)) # blue-yellow colors for WCC plots - NOTE: was missing!! viridis package also has to be installed + loaded

pb <- txtProgressBar(min = 0, max = numdyadWCC, style = 3)


# Loop through dyads
for (i in 1:numdyadWCC) { 
  
  try(if (nrow(SClist[[i]]) != 0) {  # Check if SC data is not empty
    
    print(Sys.time())
    
    # WCC analysis for SC data
    WCCList_SC <- WCC(SClist[[i]][,"Resampled_Filtered_SCL_1"], SClist[[i]][,"Resampled_Filtered_SCL_2"])
    WCCList_SC <- na.omit(t(WCCList_SC))
    WCCList_SC <- t(WCCList_SC) # columns: lags; rows: windows
    
    # Save WCC analysis result
    write.csv(WCCList_SC, file = paste("SC_D", i, "_W", wSize, "_Incr", wInc, ".csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # Peak picking
    PeakList_SC <- peakpick(WCCList_SC, pspan = LoessSpan)
    write.csv(PeakList_SC, file = paste("D", i, "_PeakList_SC.csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # Plotting
    labelsPlot_SC <- seq(6, (((nrow(WCCList_SC)-1)/(tInc/SamplingRate))+(tMax/SamplingRate)), by = ((98-1)/(tInc/SamplingRate))+(tMax/SamplingRate))
    atPlot_SC <- seq(1, nrow(WCCList_SC), by = 98)  # 98 is used as a time step for the x-axis
    if (length(atPlot_SC) < length(labelsPlot_SC)) { 
      labelsPlot_SC <- labelsPlot_SC[-length(labelsPlot_SC)]
    }
    if (length(atPlot_SC) > length(labelsPlot_SC)) {
      atPlot_SC <- atPlot_SC[-length(atPlot_SC)]
    }
    colLen <- ncol(WCCList_SC)
    
    pdf(paste("WCC_SC_D", i, "_W", wSize, "_Incr", wInc, "_L", L, ".pdf", sep = ""), width = 10, height = 5) 
    image(1:nrow(WCCList_SC), 1:ncol(WCCList_SC), WCCList_SC, col = col1(100), xlab = "Elapsed Time (seconds)", ylab = "Time Lag (seconds)", axes = F, main = paste("Dyad ", i, " - Skin Conductance"), cex.lab = 1.5, cex.main = 1.5)
    box()
    axis(2, at = seq(1, ncol(WCCList_SC), by = (ncol(WCCList_SC)-1)/(tMax/SamplingRate)), labels = seq( (-tMax/SamplingRate), (tMax/SamplingRate), by = 2), las=1)
    axis(1, at = atPlot_SC, labels = labelsPlot_SC)
    lines((PeakList_SC$maxIndex/(2*SamplingRate/tInc))+(colLen + 1)/2,lwd=3)
    abline(h = (colLen + 1)/2)
    dev.off()
    
    print(Sys.time())
    
    setTxtProgressBar(pb, i)
  })
}

close(pb)

# Save results
save("PeakList_SC", file = paste("SC_PeakList_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))

# Heart Rate -------

WCCList_HR<-vector("list", numdyadWCC)
PeakList_HR<-vector("list", numdyadWCC)

# number of samples (we have a sample rate of 20 Hz, so 20 samples = 1 sec; 1 sample = 50 ms)
wSize <- 20*8 # 8 sec
tMax <- 20*4 # 4 sec
wInc <- 20*2 # 2 sec
tInc <- 2 # 100 ms
LoessSpan <- 0.25
L <- 25
SamplingRate <- 20 # Hz

col1<-colorRampPalette(viridis(100)) 
pb <- txtProgressBar(min = 0, max = numdyadWCC, style = 3)

for (i in 1:numdyadWCC) { 
  
  try(if (nrow(HRlist[[i]]) != 0) {  # Check if HR data is not empty
    
    print(Sys.time())
    
    # WCC analysis for HR data
    WCCList_HR <- WCC(HRlist[[i]][,"IHR_1"], HRlist[[i]][,"IHR_2"])
    WCCList_HR <- na.omit(t(WCCList_HR))
    WCCList_HR <- t(WCCList_HR) # columns: lags; rows: windows
    
    # Save WCC analysis result
    write.csv(WCCList_HR, file = paste("HR_D", i, "_W", wSize, "_Incr", wInc, ".csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # Peak picking
    PeakList_HR <- peakpick(WCCList_HR, pspan = LoessSpan)
    write.csv(PeakList_HR, file = paste("D", i, "_PeakList_HR.csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # Plotting
    labelsPlot_HR <- seq(6, (((nrow(WCCList_HR)-1)/(tInc/SamplingRate))+(tMax/SamplingRate)), by = ((98-1)/(tInc/SamplingRate))+(tMax/SamplingRate))
    atPlot_HR <- seq(1, nrow(WCCList_HR), by = 98)  # 98 is used as a time step for the x-axis
    if (length(atPlot_HR) < length(labelsPlot_HR)) { 
      labelsPlot_HR <- labelsPlot_HR[-length(labelsPlot_HR)]
    }
    if (length(atPlot_HR) > length(labelsPlot_HR)) {
      atPlot_HR <- atPlot_HR[-length(atPlot_HR)]
    }
    colLen <- ncol(WCCList_HR)
    
    pdf(paste("WCC_HR_D", i, "_W", wSize, "_Incr", wInc, "_L", L, ".pdf", sep = ""), width = 10, height = 5) 
    image(1:nrow(WCCList_HR), 1:ncol(WCCList_HR), WCCList_HR, col = col1(100), xlab = "Elapsed Time (seconds)", ylab = "Time Lag (seconds)", axes = F, main = paste("Dyad ", i, " - Heart Rate"), cex.lab = 1.5, cex.main = 1.5)
    box()
    axis(2, at = seq(1, ncol(WCCList_HR), by = (ncol(WCCList_HR)-1)/(tMax/SamplingRate)), labels = seq( (-tMax/SamplingRate), (tMax/SamplingRate), by = 2), las=1)
    axis(1, at = atPlot_HR, labels = labelsPlot_HR)
    lines((PeakList_HR$maxIndex/(2*SamplingRate/tInc))+(colLen + 1)/2,lwd=3)
    abline(h = (colLen + 1)/2)
    dev.off()
    
    print(Sys.time())
    
    setTxtProgressBar(pb, i)
  })
}

#####  COMPUTE SUMMARY STATISTICS ####
SynchronyMeasures_matrix <- matrix(NA, nrow = numdyadWCC, ncol = 9)

# Set the column names for the matrix
colnames(SynchronyMeasures_matrix) <- c("dyad_i", 
                                        "WCCmean_HR", "WCCsd_HR", "TLAmean_HR", "TLAsd_HR",
                                        "WCCmean_SC", "WCCsd_SC", "TLAmean_SC", "TLAsd_SC")

# Set the sampling rate
SamplingRate = 20

# Loop through each dyad
for (i in 1:numdyadWCC){
  SynchronyMeasures_matrix[i,1] <- i
  
  # Construct the file paths for the relevant CSV files for HR and SC (no need for Face and NoFace distinction)
  HRFile <- paste("D", i, "_PeakList_HR.csv", sep = "")
  SCFile <- paste("D", i, "_PeakList_SC.csv", sep = "")
  
  # Check if both files exist
  if (file.exists(HRFile) & file.exists(SCFile)) {
    
    # Load the data
    HR_data <- read.csv(HRFile)
    SC_data <- read.csv(SCFile)
    
    # Calculate the summary statistics for HR data
    SynchronyMeasures_matrix[i,2] <- mean(HR_data$maxValue, na.rm = TRUE)
    SynchronyMeasures_matrix[i,3] <- sd(HR_data$maxValue, na.rm = TRUE)
    SynchronyMeasures_matrix[i,4] <- mean(abs(HR_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
    SynchronyMeasures_matrix[i,5] <- sd(abs(HR_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
    
    # Calculate the summary statistics for SC data
    SynchronyMeasures_matrix[i,6] <- mean(SC_data$maxValue, na.rm = TRUE)
    SynchronyMeasures_matrix[i,7] <- sd(SC_data$maxValue, na.rm = TRUE)
    SynchronyMeasures_matrix[i,8] <- mean(abs(SC_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
    SynchronyMeasures_matrix[i,9] <- sd(abs(SC_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
  }
}

# Convert the matrix to a data frame
summary_data <- as.data.frame(SynchronyMeasures_matrix)
colnames(summary_data) <- c("WCC_Dyad", 
                 "WCCmean_HR", "WCCsd_HR", "TLAmean_HR", "TLAsd_HR", 
                 "WCCmean_SC", "WCCsd_SC", "TLAmean_SC", "TLAsd_SC")

# Merge the two datasets by the 'Dyad' column
summary_data <- merge(summary_data, dyad_mapping, by = "WCC_Dyad", all.x = TRUE)

# Optional: move Real_Dyad to the front of the data frame
summary_data <- summary_data[, c("Real_Dyad", setdiff(names(summary_data), "Real_Dyad"))]

# Save the summary data frame to a CSV file
write.csv(d, file = "Synchrony-20Hz-ConditionLevel-AllDyads(56).csv", row.names = FALSE)
