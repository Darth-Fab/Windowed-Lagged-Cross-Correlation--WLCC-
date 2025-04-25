## 01-03-2025 - Fabiola Diana
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
library(viridis) 
library(psych)#NOTE: added for colour palette
options(scipen = 999) # forces R to convert scientific notations to decimals

setwd('C:/Users/fabio/Documents/GitHub/Windowed-Lagged-Cross-Correlation--WLCC-/scripts/Dyadic-Trust-Game/')

d_raw <- read.csv(file= 'data/TrustGame-ECG-SC-20Hz-EpochLevel-AllEpochs.csv')
unique(d_raw$Condition)
colnames(d_raw)
d_raw <- d_raw[is.element(d_raw$Condition, c('Face-blocked', 'Face-to-face')),]%>% droplevels
d_raw$FaceCondition <- ifelse(d_raw$Condition == "Face-to-face", 1, 0)


######### LOAD FUNCTIONS ########

source("functions/wcc-peak.R")
source("functions/mean-se.R")
source("functions/vif.R")
source("functions/rsquared.R")


######### PERFORM WCC ANALYSIS #########
d_raw_all<-d_raw
d_raw$Dyad <- as.numeric(d_raw$Dyad)
d_raw <- d_raw[order(d_raw$Dyad), ]
d_raw$Dyad_new <- match(d_raw$Dyad, unique(d_raw$Dyad))
d_raw$Dyad_new <- as.numeric(d_raw$Dyad_new)

# Create a new dataframe with only the relevant columns
dyad_mapping <- unique(d_raw[, c("Dyad", "Dyad_new")])
names(dyad_mapping)[names(dyad_mapping) == "Dyad"] <- "Real_Dyad"
names(dyad_mapping)[names(dyad_mapping) == "Dyad_new"] <- "WCC_Dyad"

# Save it to a CSV file
write.csv(dyad_mapping, "DyadMapping-TrustGame.csv", row.names = FALSE)

dyads <- unique(d_raw$Dyad_new)
numdyad <- length(dyads)

unique(d_raw$trial) # 1-2-3-4-5-6
d_raw$Trial_numeric <- as.integer(d_raw$trial)
unique(d_raw$Trial_numeric)


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
  dyad_data <- d_raw[d_raw$Dyad_new== dyads[l], ]
  
  for (t in 1:6) {
    nofaceHRlist[[t]][[l]] <- dyad_data[which(dyad_data$FaceCondition == 0 & dyad_data$Trial_numeric == t), ]
    faceHRlist[[t]][[l]] <- dyad_data[which(dyad_data$FaceCondition == 1 & dyad_data$Trial_numeric == t), ]
    nofaceSClist[[t]][[l]] <- dyad_data[which(dyad_data$FaceCondition == 0 & dyad_data$Trial_numeric == t), ]
    faceSClist[[t]][[l]] <- dyad_data[which(dyad_data$FaceCondition == 1 & dyad_data$Trial_numeric == t), ]
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
tInc = 2

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
      TLAmean_HR_Face <- mean(abs(faceHR_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
      TLAsd_HR_Face <- sd(abs(faceHR_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
      
      # Calculate the summary statistics for SC Face condition
      WCCmean_SC_Face <- mean(faceSC_data$maxValue, na.rm = TRUE)
      WCCsd_SC_Face <- sd(faceSC_data$maxValue, na.rm = TRUE)
      TLAmean_SC_Face <- mean(abs(faceSC_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
      TLAsd_SC_Face <- sd(abs(faceSC_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
      
      # Calculate the summary statistics for HR NoFace condition
      WCCmean_HR_NoFace <- mean(nofaceHR_data$maxValue, na.rm = TRUE)
      WCCsd_HR_NoFace <- sd(nofaceHR_data$maxValue, na.rm = TRUE)
      TLAmean_HR_NoFace <- mean(abs(nofaceHR_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
      TLAsd_HR_NoFace <- sd(abs(nofaceHR_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
      
      # Calculate the summary statistics for SC NoFace condition
      WCCmean_SC_NoFace <- mean(nofaceSC_data$maxValue, na.rm = TRUE)
      WCCsd_SC_NoFace <- sd(nofaceSC_data$maxValue, na.rm = TRUE)
      TLAmean_SC_NoFace <- mean(abs(nofaceSC_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
      TLAsd_SC_NoFace <- sd(abs(nofaceSC_data$maxIndex / (2 * SamplingRate / tInc)), na.rm = TRUE)
      
      # Calculate the mean of SC PPT1 and PPT2 for each condition
      SCmean_PPT1_Face <- mean(faceSC_data$PPT1, na.rm = TRUE)
      SCmean_PPT2_Face <- mean(faceSC_data$PPT2, na.rm = TRUE)
      SCmean_PPT1_NoFace <- mean(nofaceSC_data$PPT1, na.rm = TRUE)
      SCmean_PPT2_NoFace <- mean(nofaceSC_data$PPT2, na.rm = TRUE)
      
      # Calculate the mean of ECG PPT1 and PPT2 for each condition
      ECGmean_PPT1_Face <- mean(faceECG_data$PPT1, na.rm = TRUE)
      ECGmean_PPT2_Face <- mean(faceECG_data$PPT2, na.rm = TRUE)
      ECGmean_PPT1_NoFace <- mean(nofaceECG_data$PPT1, na.rm = TRUE)
      ECGmean_PPT2_NoFace <- mean(nofaceECG_data$PPT2, na.rm = TRUE)
    
      
      # Add the row to the summary data
      summary_data <- rbind(summary_data, data.frame(
        WCC_Dyad = i,
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

# Merge the two datasets by the 'Dyad' column
summary_data <- merge(summary_data, dyad_mapping, by = "WCC_Dyad", all.x = TRUE)

# Optional: move Real_Dyad to the front of the data frame
summary_data <- summary_data[, c("Real_Dyad", setdiff(names(summary_data), "Real_Dyad"))]

# Save the summary data frame to a CSV file
write.csv(summary_data, file = "SynchronyTrustGame-20Hz-TrialLevel-AllDyads(58)-OLDDATASET.csv", row.names = FALSE)

#### PSEUDOSYNCHRONY ----
# HR ------
### Face-to-face and Face-blocked together
dHR <- d_raw
fake_cor <- matrix(NA, nrow = length(unique(dHR$Dyad_new)), ncol = length(unique(dHR$Trial_numeric)))
fake_Z <- matrix(NA, nrow = length(unique(dHR$Dyad_new)), ncol = length(unique(dHR$Trial_numeric)))
minlengthcount <- matrix(NA, nrow = length(unique(dHR$Dyad_new)), ncol = length(unique(dHR$Trial_numeric)))

pb <- txtProgressBar(min=0, max = length(unique(dHR$Dyad_new)), style=3)

for (i in 1: length(unique(dHR$Dyad_new))){
  for (j in 1:length(unique(dHR$Trial_numeric))) {
    t1 <- na.omit(dHR$ECG_PPN1_HeartRate[which(dHR$Dyad_new == unique(dHR$Dyad_new)[i] & dHR$Trial_numeric == j)])
    t2 <- na.omit(dHR$ECG_PPN2_HeartRate[which(dHR$Dyad_new == unique(dHR$Dyad_new)[i + 1] & dHR$Trial_numeric == j)])
    minlength <- min(length(t1), length(t2))
    
    t1cut <- t1[1:minlength]
    t2cut <- t2[1:minlength]
    
    fake_cor[i,j] <- cor(t1cut, t2cut)
    fake_Z[i,j] <- fisherz(fake_cor[i,j])
    minlengthcount[i,j] <- minlength
    
  }
  setTxtProgressBar(pb,i)
}
close(pb)


real_cor <- matrix(NA, nrow = length(unique(dHR$Dyad_new)), ncol = length(unique(dHR$Trial_numeric)))
real_Z <- matrix(NA, nrow = length(unique(dHR$Dyad_new)), ncol = length(unique(dHR$Trial_numeric)))

pb <- txtProgressBar(min=0, max = length(unique(dHR$Dyad_new)), style=3)

for (i in 1: length(unique(dHR$Dyad_new))){
  for (j in 1:length(unique(dHR$Trial_numeric))) {
    t1 <- na.omit(dHR$ECG_PPN1_HeartRate[which(dHR$Dyad_new == unique(dHR$Dyad_new)[i] & dHR$Trial_numeric == j)])
    t2 <- na.omit(dHR$ECG_PPN2_HeartRate[which(dHR$Dyad_new == unique(dHR$Dyad_new)[i] & dHR$Trial_numeric == j)])
    if (length(t1) > 0 & length(t2) > 0){
      real_cor[i,j] <- cor(t1, t2)
      real_Z[i,j] <- fisherz(real_cor[i,j])
    } else{
      real_cor[i,j] <- NA
      real_Z[i,j] <- NA
    }
  }
  setTxtProgressBar(pb,i)
}
close(pb)


fake_Z_v <- as.vector(t(fake_Z))
fake_Z_v1 <- fake_Z_v
fake_Z_v1[is.infinite(fake_Z_v1)] <- NA 

real_Z_v <- as.vector(t(real_Z))
real_Z_v1 <- real_Z_v
real_Z_v1[is.infinite(real_Z_v1)] <- NA

t.test(na.omit(real_Z_v1), na.omit(fake_Z_v1)) # t(2.468) = 8.06, p  = .013 ; real Z mean estimate = .104 ; fake Z mean estimate = .026

plot(density(na.omit(real_Z_v1)))
lines(density(na.omit(fake_Z_v1)))

# Create a density plot comparing real vs. pseudo synchrony
plot(
  density(na.omit(real_Z_v1)), 
  col = "blue", 
  lwd = 2, 
  main = "Distribution of Synchrony (Fisher-Z Transformed)", 
  xlab = "Fisher Z-Correlation", 
  ylim = c(0, max(density(na.omit(real_Z_v1))$y, density(na.omit(fake_Z_v1))$y))
)

lines(
  density(na.omit(fake_Z_v1)), 
  col = "red", 
  lwd = 2, 
  lty = 2
)

# Add legend
legend(
  "topright", 
  legend = c("Real Dyads", "Pseudo Dyads"), 
  col = c("blue", "red"), 
  lwd = 2, 
  lty = c(1, 2),
  bty = "n"
)

# Optionally, add vertical lines for means
abline(v = mean(na.omit(real_Z_v1)), col = "blue", lwd = 2, lty = 3)
abline(v = mean(na.omit(fake_Z_v1)), col = "red", lwd = 2, lty = 3)

# SC ----
### Condition collapsed together
dSC <- d_raw
fake_cor <- matrix(NA, nrow = length(unique(dSC$Dyad_new)), ncol = length(unique(dSC$Trial_numeric)))
fake_Z <- matrix(NA, nrow = length(unique(dSC$Dyad_new)), ncol = length(unique(dSC$Trial_numeric)))
minlengthcount <- matrix(NA, nrow = length(unique(dSC$Dyad_new)), ncol = length(unique(dSC$Trial_numeric)))

pb <- txtProgressBar(min=0, max = length(unique(dSC$Dyad_new)), style=3)

for (i in 1: length(unique(dSC$Dyad_new))){
  for (j in 1:length(unique(dSC$Trial_numeric))) {
    t1 <- na.omit(dSC$SC_PPN1_SkinConductance[which(dSC$Dyad_new == unique(dSC$Dyad_new)[i] & dSC$Trial_numeric == j)])
    t2 <- na.omit(dSC$SC_PPN2_SkinConductance[which(dSC$Dyad_new == unique(dSC$Dyad_new)[i + 1] & dSC$Trial_numeric == j)])
    minlength <- min(length(t1), length(t2))
    
    t1cut <- t1[1:minlength]
    t2cut <- t2[1:minlength]
    
    fake_cor[i,j] <- cor(t1cut, t2cut)
    fake_Z[i,j] <- fisherz(fake_cor[i,j])
    minlengthcount[i,j] <- minlength
    
  }
  setTxtProgressBar(pb,i)
}
close(pb)


real_cor <- matrix(NA, nrow = length(unique(dSC$Dyad_new)), ncol = length(unique(dSC$Trial_numeric)))
real_Z <- matrix(NA, nrow = length(unique(dSC$Dyad_new)), ncol = length(unique(dSC$Trial_numeric)))

pb <- txtProgressBar(min=0, max = length(unique(dSC$Dyad_new)), style=3)

for (i in 1: length(unique(dSC$Dyad_new))){
  for (j in 1:length(unique(dSC$Trial_numeric))) {
    t1 <- na.omit(dSC$SC_PPN1_SkinConductance[which(dSC$Dyad_new == unique(dSC$Dyad_new)[i] & dSC$Trial_numeric == j)])
    t2 <- na.omit(dSC$SC_PPN2_SkinConductance[which(dSC$Dyad_new == unique(dSC$Dyad_new)[i] & dSC$Trial_numeric == j)])
    if (length(t1) > 0 & length(t2) > 0){
      real_cor[i,j] <- cor(t1, t2)
      real_Z[i,j] <- fisherz(real_cor[i,j])
    } else{
      real_cor[i,j] <- NA
      real_Z[i,j] <- NA
    }
  }
  setTxtProgressBar(pb,i)
}
close(pb)


fake_Z_v <- as.vector(t(fake_Z))
fake_Z_v1 <- fake_Z_v
fake_Z_v1[is.infinite(fake_Z_v1)] <- NA 

real_Z_v <- as.vector(t(real_Z))
real_Z_v1 <- real_Z_v
real_Z_v1[is.infinite(real_Z_v1)] <- NA

t.test(na.omit(real_Z_v1), na.omit(fake_Z_v1)) # t(2.468) = 8.06, p  = .013 ; real Z mean estimate = .104 ; fake Z mean estimate = .026

# Create a density plot comparing real vs. pseudo synchrony
plot(
  density(na.omit(real_Z_v1)), 
  col = "blue", 
  lwd = 2, 
  main = "Distribution of SCL Synchrony (Fisher-Z Transformed)", 
  xlab = "Fisher Z-Correlation", 
  ylim = c(0, max(density(na.omit(real_Z_v1))$y, density(na.omit(fake_Z_v1))$y))
)

lines(
  density(na.omit(fake_Z_v1)), 
  col = "red", 
  lwd = 2, 
  lty = 2
)

# Add legend
legend(
  "topright", 
  legend = c("Real Dyads", "Pseudo Dyads"), 
  col = c("blue", "red"), 
  lwd = 2, 
  lty = c(1, 2),
  bty = "n"
)

# Optionally, add vertical lines for means
abline(v = mean(na.omit(real_Z_v1)), col = "blue", lwd = 2, lty = 3)
abline(v = mean(na.omit(fake_Z_v1)), col = "red", lwd = 2, lty = 3)

### SC face-to-face
dSC_Face <- d_raw[d_raw$FaceCondition == 1,]
fake_cor <- matrix(NA, nrow = length(unique(dSC_Face$Dyad)), ncol = length(unique(dSC_Face$Trial_numeric)))
fake_Z <- matrix(NA, nrow = length(unique(dSC_Face$Dyad)), ncol = length(unique(dSC_Face$Trial_numeric)))
minlengthcount <- matrix(NA, nrow = length(unique(dSC_Face$Dyad)), ncol = length(unique(dSC_Face$Trial_numeric)))

pb <- txtProgressBar(min=0, max = length(unique(dSC_Face$Dyad)), style=3)

for (i in 1: length(unique(dSC_Face$Dyad))){
  for (j in 1:length(unique(dSC_Face$Trial_numeric))) {
    t1 <- na.omit(dSC_Face$SC_PPN1_SkinConductance[which(dSC_Face$Dyad == unique(dSC_Face$Dyad)[i] & dSC_Face$Trial_numeric == j)])
    t2 <- na.omit(dSC_Face$SC_PPN2_SkinConductance[which(dSC_Face$Dyad == unique(dSC_Face$Dyad)[i + 1] & dSC_Face$Trial_numeric == j)])
    minlength <- min(length(t1), length(t2))
    
    t1cut <- t1[1:minlength]
    t2cut <- t2[1:minlength]
    
    fake_cor[i,j] <- cor(t1cut, t2cut)
    fake_Z[i,j] <- fisherz(fake_cor[i,j])
    minlengthcount[i,j] <- minlength
    
  }
  setTxtProgressBar(pb,i)
}
close(pb)


real_cor <- matrix(NA, nrow = length(unique(dSC_Face$Dyad)), ncol = length(unique(dSC_Face$Trial_numeric)))
real_Z <- matrix(NA, nrow = length(unique(dSC_Face$Dyad)), ncol = length(unique(dSC_Face$Trial_numeric)))

pb <- txtProgressBar(min=0, max = length(unique(dSC_Face$Dyad)), style=3)

for (i in 1: length(unique(dSC_Face$Dyad))){
  for (j in 1:length(unique(dSC_Face$Trial_numeric))) {
    t1 <- na.omit(dSC_Face$SC_PPN1_SkinConductance[which(dSC_Face$Dyad == unique(dSC_Face$Dyad)[i] & dSC_Face$Trial_numeric == j)])
    t2 <- na.omit(dSC_Face$SC_PPN2_SkinConductance[which(dSC_Face$Dyad == unique(dSC_Face$Dyad)[i] & dSC_Face$Trial_numeric == j)])
    if (length(t1) > 0 & length(t2) > 0){
      real_cor[i,j] <- cor(t1, t2)
      real_Z[i,j] <- fisherz(real_cor[i,j])
    } else{
      real_cor[i,j] <- NA
      real_Z[i,j] <- NA
    }
  }
  setTxtProgressBar(pb,i)
}
close(pb)


fake_Z_v <- as.vector(t(fake_Z))
fake_Z_v1 <- fake_Z_v
fake_Z_v1[is.infinite(fake_Z_v1)] <- NA 

real_Z_v <- as.vector(t(real_Z))
real_Z_v1 <- real_Z_v
real_Z_v1[is.infinite(real_Z_v1)] <- NA

t.test(na.omit(real_Z_v1), na.omit(fake_Z_v1)) # t(3622.7) = 8.06, p < .001 ; real Z mean estimate = .104 ; fake Z mean estimate = .026

plot(density(na.omit(real_Z_v1)))
lines(density(na.omit(fake_Z_v1)))

# Create a density plot comparing real vs. pseudo synchrony
plot(
  density(na.omit(real_Z_v1)), 
  col = "blue", 
  lwd = 2, 
  main = "Face-to-Face - Distribution of SCL Synchrony (Fisher-Z Transformed)", 
  xlab = "Fisher Z-Correlation", 
  ylim = c(0, max(density(na.omit(real_Z_v1))$y, density(na.omit(fake_Z_v1))$y))
)

lines(
  density(na.omit(fake_Z_v1)), 
  col = "red", 
  lwd = 2, 
  lty = 2
)

# Add legend
legend(
  "topright", 
  legend = c("Real Dyads", "Pseudo Dyads"), 
  col = c("blue", "red"), 
  lwd = 2, 
  lty = c(1, 2),
  bty = "n"
)
  # small jitter to the right

# Optionally, add vertical lines for means
abline(v = mean(na.omit(real_Z_v1)), col = "blue", lwd = 2, lty = 3)
abline(v = mean(na.omit(fake_Z_v1)), col = "red", lwd = 2, lty = 3)


abline(v = real_mean, col = "blue", lwd = 2, lty = 3)
abline(v = fake_mean, col = "red", lwd = 2, lty = 3)


### SC Face-Blocked
dSC_NoFace <- d_raw[d_raw$FaceCondition == 0,]
fake_cor_2 <- matrix(NA, nrow = length(unique(dSC_NoFace$Dyad)), ncol = length(unique(dSC_NoFace$Trial_numeric)))
fake_Z_2 <- matrix(NA, nrow = length(unique(dSC_NoFace$Dyad)), ncol = length(unique(dSC_NoFace$Trial_numeric)))
minlengthcount <- matrix(NA, nrow = length(unique(dSC_NoFace$Dyad)), ncol = length(unique(dSC_NoFace$Trial_numeric)))

pb <- txtProgressBar(min=0, max = length(unique(dSC_NoFace$Dyad)), style=3)

for (i in 1: length(unique(dSC_NoFace$Dyad))){
  for (j in 1:length(unique(dSC_NoFace$Trial_numeric))) {
    t1 <- na.omit(dSC_NoFace$SC_PPN1_SkinConductance[which(dSC_NoFace$Dyad == unique(dSC_NoFace$Dyad)[i] & dSC_NoFace$Trial_numeric == j)])
    t2 <- na.omit(dSC_NoFace$SC_PPN2_SkinConductance[which(dSC_NoFace$Dyad == unique(dSC_NoFace$Dyad)[i + 1] & dSC_NoFace$Trial_numeric == j)])
    minlength <- min(length(t1), length(t2))
    
    t1cut <- t1[1:minlength]
    t2cut <- t2[1:minlength]
    
    fake_cor_2[i,j] <- cor(t1cut, t2cut)
    fake_Z_2[i,j] <- fisherz(fake_cor_2[i,j])
    minlengthcount[i,j] <- minlength
    
  }
  setTxtProgressBar(pb,i)
}
close(pb)


real_cor_2 <- matrix(NA, nrow = length(unique(dSC_NoFace$Dyad)), ncol = length(unique(dSC_NoFace$Trial_numeric)))
real_Z_2 <- matrix(NA, nrow = length(unique(dSC_NoFace$Dyad)), ncol = length(unique(dSC_NoFace$Trial_numeric)))

pb <- txtProgressBar(min=0, max = length(unique(dSC_NoFace$Dyad)), style=3)

for (i in 1: length(unique(dSC_NoFace$Dyad))){
  for (j in 1:length(unique(dSC_NoFace$Trial_numeric))) {
    t1 <- na.omit(dSC_NoFace$SC_PPN1_SkinConductance[which(dSC_NoFace$Dyad == unique(dSC_NoFace$Dyad)[i] & dSC_NoFace$Trial_numeric == j)])
    t2 <- na.omit(dSC_NoFace$SC_PPN2_SkinConductance[which(dSC_NoFace$Dyad == unique(dSC_NoFace$Dyad)[i] & dSC_NoFace$Trial_numeric == j)])
    if (length(t1) > 0 & length(t2) > 0){
      real_cor_2[i,j] <- cor(t1, t2)
      real_Z_2[i,j] <- fisherz(real_cor_2[i,j])
    } else{
      real_cor_2[i,j] <- NA
      real_Z_2[i,j] <- NA
    }
  }
  setTxtProgressBar(pb,i)
}
close(pb)


fake_Z_2_v <- as.vector(t(fake_Z_2))
fake_Z_2_v1 <- fake_Z_2_v
fake_Z_2_v1[is.infinite(fake_Z_2_v1)] <- NA 

real_Z_2_v <- as.vector(t(real_Z_2))
real_Z_2_v1 <- real_Z_2_v
real_Z_2_v1[is.infinite(real_Z_2_v1)] <- NA

t.test(na.omit(real_Z_2_v1), na.omit(fake_Z_2_v1)) # t(3622.7) = 8.06, p < .001 ; real Z mean estimate = .104 ; fake Z mean estimate = .026

plot(density(na.omit(real_Z_2_v1)))
lines(density(na.omit(fake_Z_2_v1)))

# Create a density plot comparing real vs. pseudo synchrony
plot(
  density(na.omit(real_Z_2_v1)), 
  col = "blue", 
  lwd = 2, 
  main = "Face-Blocked - Distribution of SCL Synchrony (Fisher-Z Transformed)", 
  xlab = "Fisher Z-Correlation", 
  ylim = c(0, max(density(na.omit(real_Z_2_v1))$y, density(na.omit(fake_Z_2_v1))$y))
)

lines(
  density(na.omit(fake_Z_2_v1)), 
  col = "red", 
  lwd = 2, 
  lty = 2
)

# Add legend
legend(
  "topright", 
  legend = c("Real Dyads", "Pseudo Dyads"), 
  col = c("blue", "red"), 
  lwd = 2, 
  lty = c(1, 2),
  bty = "n"
)

# Optionally, add vertical lines for means
abline(v = mean(na.omit(real_Z_2_v1)), col = "blue", lwd = 2, lty = 3)
abline(v = mean(na.omit(fake_Z_2_v1)), col = "red", lwd = 2, lty = 3)

ks.test(na.omit(real_Z_v1), na.omit(fake_Z_v1))      # Face-to-Face
ks.test(na.omit(real_Z_2_v1), na.omit(fake_Z_2_v1))  # Face-Blocked
