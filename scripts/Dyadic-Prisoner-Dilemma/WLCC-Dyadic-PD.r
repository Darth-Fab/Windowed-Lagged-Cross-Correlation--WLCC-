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
library(tidyr)
library(viridis) #NOTE: added for colour palette
options(scipen = 999) # forces R to convert scientific notations to decimals

# set working directory
setwd('C:/Users/fabio/Documents/GitHub/WCC/Analysis_FBehrens_et_al/data_pd')

d_raw <- read.csv(file= 'ecg_pd_cond_trial_withexcl.csv')
#d_raw <- read.csv(file= 'eda_pd_cond_trial_withexcl.csv')

unique(d_raw$Condition)
#d_raw$FaceCondition <- ifelse(d_raw$Condition == "Cooperative", 1, 0)
d_raw$BlockCondition <- ifelse(d_raw$Condition == "Cooperative", 1, 
                              ifelse(d_raw$Condition == "Competitive", 0, 
                                     ifelse(d_raw$Condition == "Neutral", 2, NA)))

############### FUNCTIONS ####################
# Variance Inflation Factor ------
#to check if there are variable that are covariating a lot 

vif.mer <- function (fit) { # info about the VIF: https://onlinecourses.science.psu.edu/stat501/node/347/ 
  ## adapted from rms::vif ## function: https://github.com/aufrank/R-hacks/blob/master/mer-utils.R 
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

# Mean & SE per group for plotting ------

# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions 

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
} 

# Effect size f squared --------- 

fsquared.f <- function(model.incl, model.excl) { # Aiken and West (1991) as cited in Lorah (2018) - Large-Scale Assessments in Education 

  R.model.incl <- rsquared(model.incl)$Conditional
  R.model.excl <- rsquared(model.excl)$Conditional
  fsquared <- (R.model.incl - R.model.excl)/(1 - R.model.incl)
  fsquared 
} 
# small at a value of 0.02, medium at a value of 0.15, and large at a value of 0.35 (Cohen 1992).


# WCC function -----

WCC<- function(ts1,ts2,windowSize=wSize,windowInc=wInc,tauMax=tMax,tauInc=tInc){
  #Prep data
  #Turn vecors into row
  if(is.null((dim(ts1)))){
    ts1<-t(ts1)
    ts2<-t(ts2)
  }
  #build results matrix
  resultsMatrix<-matrix(NA,nrow=(floor((NCOL(ts1)-windowSize-tauMax)/windowInc)),ncol=((tauMax/tauInc)*2)+2) ###why it's +2 and +1??
  for ( i in 1:floor((NCOL(ts1)-windowSize-tauMax)/windowInc)){
    #Lag offset value
    tau<--tauMax
    #window overlap
    window<--windowInc+i*windowInc
    #fill in results matrix
    for ( j in 1:((2*(tauMax/tauInc))+1)){
      if (tau<=0) {
        Wx<-ts1[(1+tauMax+window):(1+tauMax+(windowSize-1)+window)]
        Wy<-ts2[(1+tauMax+tau+window):(1+tauMax+(windowSize-1)+tau+window)]
      } else {
        Wx<-ts1[(1+tauMax-tau+window):(1+tauMax+(windowSize-1)-tau+window)]
        Wy<-ts2[(1+tauMax+window):(1+tauMax+(windowSize-1)+window)]
      }
      resultsMatrix[i,j]<-(1/(windowSize-1))*sum(((Wx-mean(t(Wx)))*(Wy-mean(t(Wy))))/(sd(Wx)*sd(Wy)))
      resultsMatrix[]
      tau<-tau+tauInc
    }
  }
  return(resultsMatrix)
}


# Peak Picking function ------

peakpick<- function(tAllCor, Lsize=8, graphs=0, pspan=.25, type="Max",tFileName="peak") { 
  #check for validness of parameters
  colLen <- length(tAllCor[1,]) # col length --- number of columns ## number of lags, if tmax = 10, colLen = 21; 2*(tMax/tIncr)+1
  rowLen <- length(tAllCor[,1]) # row length --- number of rows   ## (number of observations - Wsize - tMax)/ wIncr  
  tLsize <- floor((1/2)*colLen) # maximum local search region  ## tMax
  
  if(Lsize<1 || Lsize>tLsize) { # Lsize too small or too large ## too large: you would search beyond the max lag; Lsize = 0 = correlations following the max cor would not have to be smaller (Lsize = how many neighbors need to be smaller than max?)
    errorStr<- paste("Lsize should be >0 and <= ", tLsize, sep="")
    stop(errorStr) # print error message and stop the program
  }
  if(graphs<0||graphs>rowLen) { # num of graphics to print is too small or large
    errorStr <- paste("graphs should be >=0 and <= ", rowLen, sep="")
    stop(errorStr) # print error message and stop the program
  }
  if(pspan<0 || pspan>1) { # invalid pspan value ## if 0 -> no smoothing; if 1 = estimated regression line based on all data points; no windows used to estimate small segments of time serie (more precise)
    stop("pspan should be >0 and <1\n") # print error message and stop program
  }
  if(type!="Max"&&type!="Min"&&type!="max"&&type!="min"){ # only two types
    stop("valid types are: max|Max or Min|min \n") # print error message and stop
  }
  #Initilization
  drawgraph <- 0 # graphics drawed
  colLen <- length(tAllCor[1,]) # col length  
  rowLen <- length(tAllCor[,1]) # row length 
  xSequence <- seq(-(colLen-1), (colLen-1), by=1) # X axis for each graph
  mx<- rep(NA, (2*colLen-1)) #vector for keeping temp peak value for a row
  #data points will be 2*colLen-1 after smooth
  tIndex <- rep(NA, rowLen) #vector of peak index---one peak index for each row
  tValue <- rep(NA, rowLen) #vector of peak value---one peak value for each row 
  
  #type == max or Max
  if(type=="Max"||type=="max") { #compute local maximum
    for(rowNo in c(1: rowLen)) { #access each row ## for each row, see whether and where a missing is in the columns; if so, go to the next row
      #eliminate missing value
      miss <- is.na(tAllCor[rowNo, ]) #transfer a row to be T--"NA"  and F not "NA"
      #initialize the position of NA in a row 
      missposition <- 0 
      for(mIndex in c(1:colLen)) { #evaluate an entire row
        missposition <- mIndex #the position of NA
        if(miss[missposition]) break #find one
        missposition <- missposition+1 #increase count
        if(missposition==colLen+1) break #No NA in this row
      }
      #cat("missposition=", missposition, "\n")
      #if has missing value 
      if(missposition <= colLen) next #skip a row with NA
      
      else { # no missing value
        drawgraph <- drawgraph+1 #number of graph to draw
        tCor <- tAllCor[rowNo, ] #number of columns ## (tMax*2/tIncr) + 1
        #smooth
        t1 <- predict(loess(tCor~c(1:colLen), degree=2, span=pspan )) ## smooth the correlations using a quadratic estimation function https://www.youtube.com/watch?v=Vf7oJ6z2LCc  
        ## use the 25% of data points that surround a particular data point, estimate a (quadratic) regression line for these point and use the predicted point for the point of interest; do that for every (original) data point
        ## points further away from point of interest have less weight; the new estimated data points are corrected in a 2nd step to decrease weight for points that have great deviations from their original data point
        #data points is set to n
        t2 <- spline(c(1:colLen), t1, n=(2*colLen-1))$y ## increase number of data points from colLen (=2*tMax/tIncr +1) to (2*colLen-1)
        ## - 1 to have an uneven number of points again
        # show calculate progress
        # cat("row=", rowNo, "\n")
        # process a row --find max value and max index 
        windowWidth <- 0  # searched region
        lookAhead <- 0  # look ahead data points
        for(j in 1:(colLen-1)) {  # search from 1 to colLen -1
          windowWidth <- windowWidth+1 # increase search ed region
          # select the search region, the center of search region 
          # is in the middle of t2, notice that t2 has 2*colLen-1
          # data points.
          tSelect <- (colLen - windowWidth):(colLen+windowWidth) #  ## because t2 has length 2*colLen-1, colLen is the middle of t2; from there you increase the search region by 1 on both sides (windowWidth)
          mx[j] <- max(t2[tSelect], na.rm=T) # store temp max value  ## find max correlation of the current search region
          if(j==1) mmx <- mx[j] # mmx is final local max value ## store the first max value in mmx, afterwards update the mmx if it's bigger than the current mmx
          #remember current max
          else { # if j != 1
            if(mx[j]>mmx) { # new temp max value , note only one max value in tSelect
              lookAhead <- 0 # set stop search criterion to 0,
              # the criterion is that if we find  ## there must be 8 data points following the current max that are less (lower correlation) than the current max
              # Lsize data points less than current ## due to the dubbling of data point (spline function above), the 8 data points translate to 4 steps in the original lag "steps" 
              # max, then stop and the current max  ## eg: tMax = 120; tIncr = 10 --> 25 lag steps, 12 in both directions --> if Lsize = 8, then the max can be at step 8 or earlier (12 - (8/2))
              # value is the local maximum we wanted ## (at step 10, you would have to have a Lsize of (12 - (4/2)) = 10)
              mmx <- mx[j] # update new max value 
            }
            else if(mx[j]<=mmx) { # if other values are less than
              # current maximum
              # increase the count---how many neighbor data point are less than current maximum 
              lookAhead <- lookAhead+1 
              if(lookAhead>=Lsize) break # meet criterion
            }
          }#else
        }#for j--max value and index for each row
        #use match function to find the index
        Index <- match(mmx, t2[tSelect])+tSelect[1]-1  ## find the value of mmx in the selected search region of t2 and give the position of that match in t2
        # tSelect[1] is the first index of the selected window
        
        # relative position to the middle point
        position <- Index -colLen ## the scale is still doubled; Index runs from 1:(2*colLen - 1); position runs from -2*(colLen-1) to 2*(colLen) in with only even numbers so that you have colLen steps again
        ## BUT: because of the Lsize criterion, there will only be positions between -(colLen-1 - Lsize) and (colLen-1 - Lsize)  
        #according to the local maximium definition 
        if(position >(colLen - Lsize - 1) || position < (-(colLen - Lsize -1))) { # fail ## if value falls beyond the boundaries of the position --> NA (no peak found)
          tIndex[rowNo] <- NA
          tValue[rowNo] <- NA
        }
        else { # found a local maximum
          tIndex[rowNo] <- position
          tValue[rowNo] <-mmx
        }
        
        #draw plots      
        if(drawgraph <= graphs) {
          # define graphic file name
          tepsfile <- paste(tFileName, "Max", rowNo, ".eps", sep="") 
          # title of the graph
          tmain <- paste("max Index", tFileName, "r", rowNo,"w", Lsize, sep="")
          # write to postscript format
          postscript(tepsfile, height=6.4, horizontal=F)
          
          # draw borders and their labels
          plot(c(-(colLen-1), (colLen-1)), c(-1,1), xlab="Lag", ylab="Cross Correlation", main=tmain, type="n")
          
          # draw the curve
          lines(xSequence, t2, type="l")
          
          # draw the local maximum
          lines(c(position,position), c(-1,1), type="l", lty=8)
          
          # draw the axies
          lines(c(0,0), c(-1,1), type="l", lty=4)
          lines(c(-(colLen-1), (colLen-1)), c(0,0), type="l", lty=4)
          dev.off() # term off other device in order run drawing procedure
        } # if drawgraph
      } #else no missing value
    } # for rowNo-- process each row 
    #end of process a row
    return(list(maxIndex=tIndex, maxValue=tValue))
  }#if type=max
  
  ###
  # type == Min or min
  ###
  else if(type=="Min"||type=="min") {
    for(rowNo in c(1: rowLen)) { 
      #eliminate missing value
      miss <- is.na(tAllCor[rowNo, ])
      missposition <- 0
      for(mIndex in c(1:colLen)) {
        missposition <- mIndex #the position of NA
        if(miss[missposition]) break #find one
        missposition <- missposition+1
        if(missposition==colLen+1) break
      }
      #cat("missposition=", missposition, "\n")
      #if has missing value 
      if(missposition <= colLen) next #skip a row with NA
      
      else { # no missing value
        drawgraph <- drawgraph+1
        tCor <- tAllCor[rowNo, ]
        t1 <- predict(loess(tCor~c(1:colLen), degree=2, span=pspan ))
        t2 <- spline(c(1:colLen), t1, n=(2*colLen-1))$y
        # show calculate progress
        #cat("row=", rowNo, "\n")
        
        #process a row --find max value and max index 
        
        windowWidth <- 0
        lookAhead <- 0
        for(j in 1:(colLen-1)) { 
          windowWidth <- windowWidth+1
          tSelect <- (colLen - windowWidth):(colLen+windowWidth)
          mx[j] <- min(t2[tSelect], na.rm=T)
          if(j==1) mmx <- mx[j]
          #remember current max
          else {
            if(mx[j]<mmx) { #only one value in tSelect
              lookAhead <- 0
              mmx <- mx[j]
            }
            else if(mx[j]>=mmx) {
              lookAhead <- lookAhead+1
              if(lookAhead>=Lsize) break
            }
          }#else
        }#for j--max value and index for each row
        #use match function to find the index
        Index <- match(mmx, t2[tSelect])+tSelect[1]-1
        position <- Index -colLen
        
        if(position >(colLen - Lsize -1) || position < (-(colLen -Lsize -1))) { # fail
          tIndex[rowNo] <- NA
          tValue[rowNo] <- NA
        }
        else { 
          tIndex[rowNo] <- position
          tValue[rowNo] <-mmx
        }
        
        #draw first 10 plots      
        if(drawgraph <= graphs) {
          tepsfile <- paste(tFileName, "Min", rowNo, ".eps", sep="") 
          tmain <- paste( "min Index", tFileName, "r", rowNo,"w", Lsize, sep="")
          postscript(tepsfile, height=6.4, horizontal=F)
          plot(c(-(colLen-1), (colLen-1)), c(-1,1), xlab="Lag", ylab="Cross Correlation", main=tmain, type="n")
          
          lines(xSequence, t2, type="l")
          lines(c(position,position), c(-1,1), type="l", lty=8)
          lines(c(0,0), c(-1,1), type="l", lty=4)
          lines(c(-(colLen-1), (colLen-1)), c(0,0), type="l", lty=4)
          dev.off()
        } # if drawgraph 
      } #else no missing value
    } # rowNo-- process each row 
    #end of process a row
    
    return(list(minIndex=tIndex, minValue=tValue))
  }#if type=mix
}


################ PERFORM WCC ANALYSIS ###############
d_raw_all<-d_raw
num = seq(1,5,by=1)
num2 = seq(11,15, by=1)
num3= seq(12,15, by=1)

#Condition Level
#d_raw <- d_raw_all[is.element(d_raw_all$Condition, c('Cooperative', 'Competitive')),]%>% droplevels

#Only selected dyad
#d_raw <- d_raw[is.element(d_raw$Dyad, "5"),]%>%droplevels
#d_raw <- d_raw[!is.element(d_raw$Dyad, "5"),] %>% droplevels

d_raw <- d_raw %>% drop_na()

#d_raw <- d_raw[is.element(d_raw$Dyad),]%>%droplevels
# to run the analysis, we first make lists containing matrices with the raw HR and SCL data per face condition for each dyad

numdyad = length(unique(d_raw$Dyad))

copHRlist <- vector("list", numdyad)
compHRlist <- vector("list", numdyad)
neuHRlist <-vector("list", numdyad)
copSClist <- vector("list", numdyad)
compSClist <- vector("list", numdyad)
neuSClist <-vector("list", numdyad)

copHRmatrix <- matrix(NA, nrow = 160000, ncol = ncol(d_raw))
compHRmatrix <- matrix(NA, nrow = 160000, ncol = ncol(d_raw))
neuHRmatrix <- matrix(NA, nrow = 160000, ncol = ncol(d_raw))
copSCmatrix <- matrix(NA, nrow = 160000, ncol = ncol(d_raw))
compSCmatrix <- matrix(NA, nrow = 160000, ncol = ncol(d_raw))
neuSCmatrix <- matrix(NA, nrow = 160000, ncol = ncol(d_raw))



# this is not working - take from d_raw instead of d (at least values in table then)

l = 1
for (i in unique(d_raw$Dyad)){
  copHRmatrix <- d_raw[which(d_raw$BlockCondition == 1 & d_raw$Dyad == i & d_raw$Excl_HR == 0),]
  compHRmatrix <- d_raw[which(d_raw$BlockCondition == 0 & d_raw$Dyad == i & d_raw$Excl_HR == 0),]
  neuHRmatrix <- d_raw[which(d_raw$BlockCondition == 2 & d_raw$Dyad == i & d_raw$Excl_HR == 0),]
  copSCmatrix <- d_raw[which(d_raw$BlockCondition == 1 & d_raw$Dyad == i & d_raw$Excl_SC == 0),]
  compSCmatrix <- d_raw[which(d_raw$BlockCondition == 0 & d_raw$Dyad == i & d_raw$Excl_SC == 0),]
  neuSCmatrix <- d_raw[which(d_raw$BlockCondition == 2 & d_raw$Dyad == i & d_raw$Excl_HR == 0),]
  copHRlist[[l]] <- copHRmatrix
  compHRlist[[l]] <- compHRmatrix
  neuHRlist[[l]] <- neuHRmatrix
  copSClist[[l]] <- copSCmatrix
  compSClist[[l]] <- compSCmatrix
  neuSClist[[l]] <- neuSCmatrix
  l = l + 1
}

numdyadWCC <- length(unique(d_raw$Dyad))

# Skin Conductance ------ 
copWCCList_SC<-vector("list", numdyadWCC)
copPeakList_SC<-vector("list", numdyadWCC)

compWCCList_SC <- vector("list", numdyadWCC)
compPeakList_SC <- vector("list", numdyadWCC)

neuWCCList_SC <- vector("list", numdyadWCC)
neuPeakList_SC <- vector("list", numdyadWCC)

# number of samples (we have a sample rate of 20 Hz, so 20 samples = 1 sec; 1 sample = 50 ms)
wSize <- 20*8 # 12 sec
tMax <- 20*4 # 12 sec
wInc <- 20*2 # 2 sec
tInc <- 2 # 100 ms
LoessSpan <- 0.25
L <- 25
SamplingRate <- 20 # Hz

# NOTE: error:  Error in col1(100) : could not find function "col1" - colour palette was missing

col1<-colorRampPalette(viridis(100)) # blue-yellow colors for WCC plots - NOTE: was missing!! viridis package also has to be installed + loaded

pb <- txtProgressBar(min = 0, max = numdyadWCC, style = 3)

for (i in 1:numdyadWCC){ 
  
  if (nrow(copSClist[[i]]) != 0) {
    
    
    print(Sys.time())
    # WCC analysis
    copWCCList_SC[[i]]<-WCC(copSClist[[i]][,"SC_PPN1_SkinConductance"],copSClist[[i]][,"SC_PPN2_SkinConductance"])
    copWCCList_SC[[i]] <- na.omit(t(copWCCList_SC[[i]]))
    copWCCList_SC[[i]] <- t(copWCCList_SC[[i]]) # columns: lags; rows: windows
    write.csv(copWCCList_SC[[i]], file = paste("SC_D", i, "_Cop_W", wSize, "_Incr", wInc, ".csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # peak picking
    copPeakList_SC[[i]] <- peakpick(copWCCList_SC[[i]], pspan = LoessSpan, graphs = 1)
    write.csv(copPeakList_SC[[i]], file = paste("SC_D", i, "_Cop_W", wSize, "_Incr", wInc, "_L", L, "_peak.csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # make the plot
    
    ## get correct markers and labels on the x-axis (transforming time in nrow to time in sec)
    labelsPlot_SC_cop <- seq(6, (((nrow(copWCCList_SC[[i]])-1)/(tInc/SamplingRate))+(tMax/SamplingRate)), by = ((98-1)/(tInc/SamplingRate))+(tMax/SamplingRate))
    atPlot_SC_cop <- seq(1, nrow(copWCCList_SC[[i]]), by = 98) ## 98 : t(sec) = ((t(nrow) -1)/tInc/SamplingRate) + tMax/SamplingRate; e.g. the equivalent of 200 sec are 98 rows in WCC matrix: 200 = ((t(nrow) - 1)/0.5) + 6
    if (length(atPlot_SC_cop) < length(labelsPlot_SC_cop)) { ## sometimes steps of the markers and labels are in such a way that one of them has one more step; in such cases the longer interval is trimmed to the length of the other
      labelsPlot_SC_cop <- labelsPlot_SC_cop[-length(labelsPlot_SC_cop)]
    }
    if (length(atPlot_SC_cop) > length(labelsPlot_SC_cop)) {
      atPlot_SC_cop <- atPlot_SC_cop[-length(atPlot_SC_cop)]
      
    } 
    
    colLen <- ncol(copWCCList_SC[[i]])
    
    pdf( paste("WCC_SC_D", i, "_Cop_W", wSize, "_Incr", wInc, "_L", L, ".pdf", sep = ""), width = 10, height = 5) 
    image(1:nrow(copWCCList_SC[[i]]), 1:ncol(copWCCList_SC[[i]]),copWCCList_SC[[i]], col = col1(100), xlab = "Elapsed Time (seconds)", ylab = "Time Lag (seconds)", axes = F, main = paste("Dyad ", i, " - Cooperative Condition - Skin Conductance"), cex.lab = 1.5, cex.main = 1.5)
    box()
    axis(2, at = seq(1, ncol(copWCCList_SC[[i]]), by = (ncol(copWCCList_SC[[i]])-1)/(tMax/SamplingRate)), labels = seq( (-tMax/SamplingRate), (tMax/SamplingRate), by = 2), las=1) ## by = 2 is only possible if tMax/samplingrate is an even number
    axis(1, at = atPlot_SC_cop, labels = labelsPlot_SC_cop)
    lines((copPeakList_SC[[i]]$maxIndex/(2*SamplingRate/tInc))+(colLen + 1)/2,lwd=3) ## (colLen + 1)/2 = if colLen = 25, then 13 is the middle point (12 steps to both sides)
    abline(h = (colLen + 1)/2)
    dev.off()
    
  }
  
  if (nrow(compSClist[[i]]) != 0) { 
    
    compWCCList_SC[[i]]<-WCC(compSClist[[i]][,"SC_PPN1_SkinConductance"],compSClist[[i]][,"SC_PPN2_SkinConductance"])
    compWCCList_SC[[i]] <- na.omit(t(compWCCList_SC[[i]]))
    compWCCList_SC[[i]] <- t(compWCCList_SC[[i]])
    write.csv(compWCCList_SC[[i]], file = paste("SC_D", i, "_Comp_W", wSize, "_Incr", wInc, ".csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # peak picking
    
    compPeakList_SC[[i]] <- peakpick(compWCCList_SC[[i]], pspan = LoessSpan)
    write.csv(compPeakList_SC[[i]], file = paste("SC_D", i, "_Comp_W", wSize, "_Incr", wInc, "_L", L, "_peak.csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # make the plot
    labelsPlot_SC_comp <- seq(6, (((nrow(compWCCList_SC[[i]])-1)/(tInc/SamplingRate))+(tMax/SamplingRate)), by = ((98-1)/(tInc/SamplingRate))+(tMax/SamplingRate))
    atPlot_SC_comp <- seq(1, nrow(compWCCList_SC[[i]]), by = 98) ## 98 : t(sec) = ((t(nrow) -1)/tInc/SamplingRate) + tMax/SamplingRate; e.g. the equivalent of 200 sec are 98 rows in WCC matrix: 200 = ((t(nrow) - 1)/0.5) + 6
    if (length(atPlot_SC_comp) < length(labelsPlot_SC_comp)) { ## sometimes steps of the markers and labels are in such a way that one of them has one more step; in such cases the longer interval is trimmed to the length of the other
      labelsPlot_SC_comp <- labelsPlot_SC_comp[-length(labelsPlot_SC_comp)]
    }
    if (length(atPlot_SC_comp) > length(labelsPlot_SC_comp)) {
      atPlot_SC_comp <- atPlot_SC_comp[-length(atPlot_SC_comp)]
    } 
    
    colLen <- ncol(compWCCList_SC[[i]])
    
    pdf( paste("WCC_SC_D", i, "_Comp_W", wSize, "_Incr", wInc, "_L", L, ".pdf", sep = ""), width = 10, height = 5) 
    image(1:nrow(compWCCList_SC[[i]]), 1:ncol(compWCCList_SC[[i]]),compWCCList_SC[[i]], col = col1(100), xlab = "Elapsed Time (seconds)", ylab = "Time Lag (seconds)", axes = F, main = paste("Dyad ", i, " - Competitive Condition - Skin Conductance"), cex.lab = 1.5, cex.main = 1.5)
    box()
    axis(2, at = seq(1, ncol(compWCCList_SC[[i]]), by = (ncol(compWCCList_SC[[i]])-1)/(tMax/SamplingRate)), labels = seq( (-tMax/SamplingRate), (tMax/SamplingRate), by = 2), las=1) ## by = 2 is only possible if tMax/samplingrate is an even number
    axis(1, at = atPlot_SC_comp, labels = labelsPlot_SC_comp)
    lines((compPeakList_SC[[i]]$maxIndex/(2*SamplingRate/tInc))+(colLen + 1)/2,lwd=3) ## (colLen + 1)/2 = if colLen = 25, then 13 is the middle point (12 steps to both sides)
    abline(h = (colLen + 1)/2)
    dev.off()
    
  }
  if (nrow(neuSClist[[i]]) != 0) {
    
    
    print(Sys.time())
    # WCC analysis
    neuWCCList_SC[[i]]<-WCC(neuSClist[[i]][,"SC_PPN1_SkinConductance"],neuSClist[[i]][,"SC_PPN2_SkinConductance"])
    neuWCCList_SC[[i]] <- na.omit(t(neuWCCList_SC[[i]]))
    neuWCCList_SC[[i]] <- t(neuWCCList_SC[[i]]) # columns: lags; rows: windows
    write.csv(neuWCCList_SC[[i]], file = paste("SC_D", i, "_Neu_W", wSize, "_Incr", wInc, ".csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # peak picking
    neuPeakList_SC[[i]] <- peakpick(neuWCCList_SC[[i]], pspan = LoessSpan, graphs = 1)
    write.csv(neuPeakList_SC[[i]], file = paste("SC_D", i, "_Neu_W", wSize, "_Incr", wInc, "_L", L, "_peak.csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # make the plot
    
    ## get correct markers and labels on the x-axis (transforming time in nrow to time in sec)
    labelsPlot_SC_neu <- seq(6, (((nrow(neuWCCList_SC[[i]])-1)/(tInc/SamplingRate))+(tMax/SamplingRate)), by = ((98-1)/(tInc/SamplingRate))+(tMax/SamplingRate))
    atPlot_SC_neu <- seq(1, nrow(neuWCCList_SC[[i]]), by = 98) ## 98 : t(sec) = ((t(nrow) -1)/tInc/SamplingRate) + tMax/SamplingRate; e.g. the equivalent of 200 sec are 98 rows in WCC matrix: 200 = ((t(nrow) - 1)/0.5) + 6
    if (length(atPlot_SC_neu) < length(labelsPlot_SC_neu)) { ## sometimes steps of the markers and labels are in such a way that one of them has one more step; in such cases the longer interval is trimmed to the length of the other
      labelsPlot_SC_neu <- labelsPlot_SC_neu[-length(labelsPlot_SC_neu)]
    }
    if (length(atPlot_SC_neu) > length(labelsPlot_SC_neu)) {
      atPlot_SC_neu <- atPlot_SC_neu[-length(atPlot_SC_neu)]
      
    } 
    
    colLen <- ncol(neuWCCList_SC[[i]])
    
    pdf( paste("WCC_SC_D", i, "_Neu_W", wSize, "_Incr", wInc, "_L", L, ".pdf", sep = ""), width = 10, height = 5) 
    image(1:nrow(neuWCCList_SC[[i]]), 1:ncol(neuWCCList_SC[[i]]),neuWCCList_SC[[i]], col = col1(100), xlab = "Elapsed Time (seconds)", ylab = "Time Lag (seconds)", axes = F, main = paste("Dyad ", i, " - Neutral Condition - Skin Conductance"), cex.lab = 1.5, cex.main = 1.5)
    box()
    axis(2, at = seq(1, ncol(neuWCCList_SC[[i]]), by = (ncol(neuWCCList_SC[[i]])-1)/(tMax/SamplingRate)), labels = seq( (-tMax/SamplingRate), (tMax/SamplingRate), by = 2), las=1) ## by = 2 is only possible if tMax/samplingrate is an even number
    axis(1, at = atPlot_SC_neu, labels = labelsPlot_SC_neu)
    lines((neuPeakList_SC[[i]]$maxIndex/(2*SamplingRate/tInc))+(colLen + 1)/2,lwd=3) ## (colLen + 1)/2 = if colLen = 25, then 13 is the middle point (12 steps to both sides)
    abline(h = (colLen + 1)/2)
    dev.off()
    
  }
  
  
  print(Sys.time())
  setTxtProgressBar(pb, i)
  
}    

close(pb)

save("copPeakList_SC", file = paste("SC_PeakList_Cop_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))
save("compPeakList_SC", file = paste("SC_PeakList_Comp_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))
save("neuPeakList_SC", file = paste("SC_PeakList_Neu_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))

# Heart Rate -------

copWCCList_HR<-vector("list", numdyadWCC)
copPeakList_HR<-vector("list", numdyadWCC)

compWCCList_HR <- vector("list", numdyadWCC)
compPeakList_HR <- vector("list", numdyadWCC)

neuWCCList_HR <- vector("list", numdyadWCC)
neuPeakList_HR <- vector("list", numdyadWCC)



# number of samples (we have a sample rate of 20 Hz, so 20 samples = 1 sec; 1 sample = 50 ms)
wSize <- 20*8 # 8 sec
tMax <- 20*4 # 4 sec
wInc <- 20*2 # 2 sec
tInc <- 2 # 100 ms
LoessSpan <- 0.25
L <- 25
SamplingRate <- 20 # Hz

numdyadWCC <- length(unique(d_raw$Dyad)) # NOTE: gave 0 wit d$dyad, therefore d adjusted to d$Dyad

col1<-colorRampPalette(viridis(100)) 
pb <- txtProgressBar(min = 0, max = numdyadWCC, style = 3)


for (i in 1:numdyadWCC){ 
  
  #cop
  if (nrow(copHRlist[[i]]) != 0) {
    
    #print('cop')
    #print(copHRlist[[i]][,"ECG_PPN1_HeartRate"])
    #print(copHRlist[[i]][,"ECG_PPN2_HeartRate"])
    
    # WCC analysis
    copWCCList_HR[[i]]<-WCC(copHRlist[[i]][,"ECG_PPN1_HeartRate"],copHRlist[[i]][,"ECG_PPN2_HeartRate"])
    copWCCList_HR[[i]] <- na.omit(t(copWCCList_HR[[i]]))
    copWCCList_HR[[i]] <- t(copWCCList_HR[[i]]) # columns: lags; rows: windows
    write.csv(copWCCList_HR[[i]], file = paste("HR_D", i, "_Cop_W", wSize, "_Incr", wInc, ".csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    
    # peak picking
    copPeakList_HR[[i]] <- peakpick(copWCCList_HR[[i]], pspan = LoessSpan, graphs = 1)
    write.csv(copPeakList_HR[[i]], file = paste("HR_D", i, "_Cop_W", wSize, "_Incr", wInc, "_L", L, "_peak.csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # make the plot
    
    ## get correct markers and labels on the x-axis (transforming time in nrow to time in sec)
    # cop
    labelsPlot_HR_cop <- seq(6, (((nrow(copWCCList_HR[[i]])-1)/(tInc/SamplingRate))+(tMax/SamplingRate)), by = ((98-1)/(tInc/SamplingRate))+(tMax/SamplingRate))
    atPlot_HR_cop <- seq(1, nrow(copWCCList_HR[[i]]), by = 98) ## 98 : t(sec) = ((t(nrow) -1)/tInc/SamplingRate) + tMax/SamplingRate; e.g. the equivalent of 200 sec are 98 rows in WCC matrix: 200 = ((t(nrow) - 1)/0.5) + 6
    if (length(atPlot_HR_cop) < length(labelsPlot_HR_cop)) { ## sometimes steps of the markers and labels are in such a way that one of them has one more step; in such cases the longer interval is trimmed to the length of the other
      labelsPlot_HR_cop <- labelsPlot_HR_cop[-length(labelsPlot_HR_cop)]
    }
    if (length(atPlot_HR_cop) > length(labelsPlot_HR_cop)) {
      atPlot_HR_cop <- atPlot_HR_cop[-length(atPlot_HR_cop)]
    } 
    
    colLen <- ncol(copWCCList_HR[[i]])
    
    pdf( paste("WCC_HR_D", i, "_Cop_W", wSize, "_Incr", wInc, "_L", L, ".pdf", sep = ""), width = 10, height = 5) 
    image(1:nrow(copWCCList_HR[[i]]), 1:ncol(copWCCList_HR[[i]]),copWCCList_HR[[i]], col = col1(100), xlab = "Elapsed Time (seconds)", ylab = "Time Lag (seconds)", axes = F, main = paste("Dyad ", i, " - Cooperative Condition - Heart Rate"), cex.lab = 1.5, cex.main = 1.5)
    box()
    axis(2, at = seq(1, ncol(copWCCList_HR[[i]]), by = (ncol(copWCCList_HR[[i]])-1)/(tMax/SamplingRate)), labels = seq( (-tMax/SamplingRate), (tMax/SamplingRate), by = 2), las=1) ## by = 2 is only possible if tMax/samplingrate is an even number
    #axis(1, at = atPlot_HR_cop, labels = labelsPlot_HR_cop)
    lines((copPeakList_HR[[i]]$maxIndex/(2*SamplingRate/tInc))+(colLen + 1)/2,lwd=3) ## (colLen + 1)/2 = if colLen = 25, then 13 is the middle point (12 steps to both sides)
    abline(h = (colLen + 1)/2)
    dev.off()
    
  }
  if (nrow(compHRlist[[i]]) != 0) {
    
    #WCC analysis
    compWCCList_HR[[i]]<-WCC(compHRlist[[i]][,"ECG_PPN1_HeartRate"],compHRlist[[i]][,"ECG_PPN2_HeartRate"])
    compWCCList_HR[[i]] <- na.omit(t(compWCCList_HR[[i]]))
    compWCCList_HR[[i]] <- t(compWCCList_HR[[i]])
    write.csv(compWCCList_HR[[i]], file = paste("HR_D", i, "_Comp_W", wSize, "_Incr", wInc, ".csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # peak picking
    
    compPeakList_HR[[i]] <- peakpick(compWCCList_HR[[i]], pspan = LoessSpan)
    write.csv(compPeakList_HR[[i]], file = paste("HR_D", i, "_Comp_W", wSize, "_Incr", wInc, "_L", L, "_peak.csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # make the plot
    
    ## get correct markers and labels on the x-axis (transforming time in nrow to time in sec)
    # comp
    labelsPlot_HR_comp <- seq(6, (((nrow(compWCCList_HR[[i]])-1)/(tInc/SamplingRate))+(tMax/SamplingRate)), by = ((98-1)/(tInc/SamplingRate))+(tMax/SamplingRate))
    atPlot_HR_comp <- seq(1, nrow(compWCCList_HR[[i]]), by = 98) ## 98 : t(sec) = ((t(nrow) -1)/tInc/SamplingRate) + tMax/SamplingRate; e.g. the equivalent of 200 sec are 98 rows in WCC matrix: 200 = ((t(nrow) - 1)/0.5) + 6
    if (length(atPlot_HR_comp) < length(labelsPlot_HR_comp)) { ## sometimes steps of the markers and labels are in such a way that one of them has one more step; in such cases the longer interval is trimmed to the length of the other
      labelsPlot_HR_comp <- labelsPlot_HR_comp[-length(labelsPlot_HR_comp)]
    }
    if (length(atPlot_HR_comp) > length(labelsPlot_HR_comp)) {
      atPlot_HR_comp <- atPlot_HR_comp[-length(atPlot_HR_comp)]
    } 
    
    colLen <- ncol(compWCCList_HR[[i]])
    
    pdf( paste("WCC_HR_D", i, "_Comp_W", wSize, "_Incr", wInc, "_L", L, ".pdf", sep = ""), width = 10, height = 5) 
    image(1:nrow(compWCCList_HR[[i]]), 1:ncol(compWCCList_HR[[i]]),compWCCList_HR[[i]], col = col1(100), xlab = "Elapsed Time (seconds)", ylab = "Time Lag (seconds)", axes = F, main = paste("Dyad ", i, " - Competitive Condition - Heart rate"), cex.lab = 1.5, cex.main = 1.5)
    box()
    axis(2, at = seq(1, ncol(compWCCList_HR[[i]]), by = (ncol(compWCCList_HR[[i]])-1)/(tMax/SamplingRate)), labels = seq( (-tMax/SamplingRate), (tMax/SamplingRate), by = 2), las=1) ## by = 2 is only possible if tMax/samplingrate is an even number
    #axis(1, at = atPlot_HR_comp, labels = labelsPlot_HR_comp)
    lines((compPeakList_HR[[i]]$maxIndex/(2*SamplingRate/tInc))+(colLen + 1)/2,lwd=3) ## (colLen + 1)/2 = if colLen = 25, then 13 is the middle point (12 steps to both sides)
    abline(h = (colLen + 1)/2)
    dev.off()
    
  }
  
  if (nrow(neuHRlist[[i]]) != 0) {
    
    #print('neu')
    #print(neuHRlist[[i]][,"ECG_PPN1_HeartRate"])
    #print(neuHRlist[[i]][,"ECG_PPN2_HeartRate"])
    
    # WCC analysis
    neuWCCList_HR[[i]]<-WCC(neuHRlist[[i]][,"ECG_PPN1_HeartRate"],neuHRlist[[i]][,"ECG_PPN2_HeartRate"])
    neuWCCList_HR[[i]] <- na.omit(t(neuWCCList_HR[[i]]))
    neuWCCList_HR[[i]] <- t(neuWCCList_HR[[i]]) # columns: lags; rows: windows
    write.csv(neuWCCList_HR[[i]], file = paste("HR_D", i, "_neu_W", wSize, "_Incr", wInc, ".csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    
    # peak picking
    neuPeakList_HR[[i]] <- peakpick(neuWCCList_HR[[i]], pspan = LoessSpan, graphs = 1)
    write.csv(neuPeakList_HR[[i]], file = paste("HR_D", i, "_neu_W", wSize, "_Incr", wInc, "_L", L, "_peak.csv", sep = ""), row.names = FALSE)
    print(Sys.time())
    
    # make the plot
    
    ## get correct markers and labels on the x-axis (transforming time in nrow to time in sec)
    # neu
    labelsPlot_HR_neu <- seq(6, (((nrow(neuWCCList_HR[[i]])-1)/(tInc/SamplingRate))+(tMax/SamplingRate)), by = ((98-1)/(tInc/SamplingRate))+(tMax/SamplingRate))
    atPlot_HR_neu <- seq(1, nrow(neuWCCList_HR[[i]]), by = 98) ## 98 : t(sec) = ((t(nrow) -1)/tInc/SamplingRate) + tMax/SamplingRate; e.g. the equivalent of 200 sec are 98 rows in WCC matrix: 200 = ((t(nrow) - 1)/0.5) + 6
    if (length(atPlot_HR_neu) < length(labelsPlot_HR_neu)) { ## sometimes steps of the markers and labels are in such a way that one of them has one more step; in such cases the longer interval is trimmed to the length of the other
      labelsPlot_HR_neu <- labelsPlot_HR_neu[-length(labelsPlot_HR_neu)]
    }
    if (length(atPlot_HR_neu) > length(labelsPlot_HR_neu)) {
      atPlot_HR_neu <- atPlot_HR_neu[-length(atPlot_HR_neu)]
    } 
    
    colLen <- ncol(neuWCCList_HR[[i]])
    
    pdf( paste("WCC_HR_D", i, "_neu_W", wSize, "_Incr", wInc, "_L", L, ".pdf", sep = ""), width = 10, height = 5) 
    image(1:nrow(neuWCCList_HR[[i]]), 1:ncol(neuWCCList_HR[[i]]),neuWCCList_HR[[i]], col = col1(100), xlab = "Elapsed Time (seconds)", ylab = "Time Lag (seconds)", axes = F, main = paste("Dyad ", i, " - Neutral Condition - Heart Rate"), cex.lab = 1.5, cex.main = 1.5)
    box()
    axis(2, at = seq(1, ncol(neuWCCList_HR[[i]]), by = (ncol(neuWCCList_HR[[i]])-1)/(tMax/SamplingRate)), labels = seq( (-tMax/SamplingRate), (tMax/SamplingRate), by = 2), las=1) ## by = 2 is only possible if tMax/samplingrate is an even number
    #axis(1, at = atPlot_HR_neu, labels = labelsPlot_HR_neu)
    lines((neuPeakList_HR[[i]]$maxIndex/(2*SamplingRate/tInc))+(colLen + 1)/2,lwd=3) ## (colLen + 1)/2 = if colLen = 25, then 13 is the middle point (12 steps to both sides)
    abline(h = (colLen + 1)/2)
    dev.off()
    
  }
  
  print(Sys.time())
  setTxtProgressBar(pb, i)
  
}    

close(pb)

save("copPeakList_HR", file = paste("HR_PeakList_Cop_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))
save("compPeakList_HR", file = paste("HR_PeakList_Comp_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))
save("neuPeakList_HR", file = paste("HR_PeakList_Neu_W", wSize, "_Incr", wInc, "_L", L, ".RData", sep = ""))



# Compute summary statistics and save results in new data file (d = DataFile_Analyses_FBehrens_et_al.RData) -----

#load the peak lists
# HR
load("~/GitHub/WCC/Analysis_FBehrens_et_al/data_pd/HR_PeakList_Comp_W160_Incr40_L25.RData")
load("~/GitHub/WCC/Analysis_FBehrens_et_al/data_pd/HR_PeakList_Cop_W160_Incr40_L25.RData")
load("~/GitHub/WCC/Analysis_FBehrens_et_al/data_pd/HR_PeakList_Neu_W160_Incr40_L25.RData")
tInc_HR = 20*0.5
#SC
load("~/GitHub/WCC/Analysis_FBehrens_et_al/data_pd/SC_PeakList_Comp_W160_Incr40_L25.RData")
load("~/GitHub/WCC/Analysis_FBehrens_et_al/data_pd/SC_PeakList_Cop_W160_Incr40_L25.RData")
load("~/GitHub/WCC/Analysis_FBehrens_et_al/data_pd/SC_PeakList_Neu_W160_Incr40_L25.RData")
tInc_SC = 20*2





#numdyadWCC = 64

SynchronyMeasures_matrix <- matrix(NA, nrow = numdyadWCC, ncol = 41)

colnames(SynchronyMeasures_matrix) <- c("dyad_i", "WCCmean_HR_Cop", "WCCsd_HR_Cop", "TLAmean_HR_Cop", "TLAsd_HR_Cop",
                                        "WCCmean_SC_Face", "WCCsd_SC_Face", "TLAmean_SC_Face", "TLAsd_SC_Face",
                                        "WCCmean_Corr_Face", "WCCsd_Corr_Face", "TLAmean_Corr_Face", "TLAsd_Corr_Face",
                                        "WCCmean_Zyg_Face", "WCCsd_Zyg_Face", "TLAmean_Zyg_Face", "TLAsd_Zyg_Face",
                                        "WCCmean_SKT_Face", "WCCsd_SKT_Face", "TLAmean_SKT_Face", "TLAsd_SKT_Face",
                                        "WCCmean_HR_NoFace", "WCCsd_HR_NoFace", "TLAmean_HR_NoFace", "TLAsd_HR_NoFace", 
                                        "WCCmean_SC_NoFace", "WCCsd_SC_NoFace", "TLAmean_SC_NoFace", "TLAsd_SC_NoFace", 
                                        "WCCmean_Corr_NoFace", "WCCsd_Corr_NoFace", "TLAmean_Corr_NoFace", "TLAsd_Corr_NoFace", 
                                        "WCCmean_Zyg_NoFace", "WCCsd_Zyg_NoFace", "TLAmean_Zyg_NoFace", "TLAsd_Zyg_NoFace", 
                                        "WCCmean_SKT_NoFace", "WCCsd_SKT_NoFace", "TLAmean_SKT_NoFace", "TLAsd_SKT_NoFace")

# in the current analysis, only the WCCmean measure is used quantifying the strength of synchrony in this manuscript and is included in the data file
# WCCsd variables represent the standard deviation of the peak correlations within a face condition per dyad 
# TLA refers to the corresponding time lag at which the peak correlation is identified 


## !scale of time lags is in seconds! 

SamplingRate = 20
#tInc = 2 # question: 2 means 0.5? adjust accordingly in loop for each measure

for (i in 1:numdyadWCC){
  SynchronyMeasures_matrix[i,1] <- i

#### here we have all the values of WCC for each condition, which is good if you want to do ANOVA or t-test
#### but if you want to do Linear Mixed Model as we want, then this structure is not okay 
  # Face
  SynchronyMeasures_matrix[i,2] <- mean(FacePeakList_HR[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,3] <- sd(FacePeakList_HR[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,4] <- mean(abs(FacePeakList_HR[[i]]$maxIndex/(2*SamplingRate/tInc_HR)), na.rm = T)
  SynchronyMeasures_matrix[i,5] <- sd(abs(FacePeakList_HR[[i]]$maxIndex/(2*SamplingRate/tInc_HR)), na.rm = T)
  
  SynchronyMeasures_matrix[i,6] <- mean(FacePeakList_SC[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,7] <- sd(FacePeakList_SC[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,8] <- mean(abs(FacePeakList_SC[[i]]$maxIndex/(2*SamplingRate/tInc_SC)), na.rm = T)
  SynchronyMeasures_matrix[i,9] <- sd(abs(FacePeakList_SC[[i]]$maxIndex/(2*SamplingRate/tInc_SC)), na.rm = T)
  
  SynchronyMeasures_matrix[i,10] <- mean(FacePeakList_Corr[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,11] <- sd(FacePeakList_Corr[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,12] <- mean(abs(FacePeakList_Corr[[i]]$maxIndex/(2*SamplingRate/tInc_EMG)), na.rm = T)
  SynchronyMeasures_matrix[i,13] <- sd(abs(FacePeakList_Corr[[i]]$maxIndex/(2*SamplingRate/tInc_EMG)), na.rm = T)
  
  SynchronyMeasures_matrix[i,14] <- mean(FacePeakList_Zyg[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,15] <- sd(FacePeakList_Zyg[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,16] <- mean(abs(FacePeakList_Zyg[[i]]$maxIndex/(2*SamplingRate/tInc_EMG)), na.rm = T)
  SynchronyMeasures_matrix[i,17] <- sd(abs(FacePeakList_Zyg[[i]]$maxIndex/(2*SamplingRate/tInc_EMG)), na.rm = T)
  
  SynchronyMeasures_matrix[i,18] <- mean(FacePeakList_SKT[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,19] <- sd(FacePeakList_SKT[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,20] <- mean(abs(FacePeakList_SKT[[i]]$maxIndex/(2*SamplingRate/tInc_SKT)), na.rm = T)
  SynchronyMeasures_matrix[i,21] <- sd(abs(FacePeakList_SKT[[i]]$maxIndex/(2*SamplingRate/tInc_SKT)), na.rm = T)
  
  
  # NoFace
  
  SynchronyMeasures_matrix[i,22] <- mean(NoFacePeakList_HR[[i]]$maxValue, na.rm = T) # 2 NA's = geen max WCC gevonden 
  SynchronyMeasures_matrix[i,23] <- sd(NoFacePeakList_HR[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,24] <- mean(abs(NoFacePeakList_HR[[i]]$maxIndex/(2*SamplingRate/tInc_HR)), na.rm = T)
  SynchronyMeasures_matrix[i,25] <- sd(abs(NoFacePeakList_HR[[i]]$maxIndex/(2*SamplingRate/tInc_HR)), na.rm = T)
  
  SynchronyMeasures_matrix[i,26] <- mean(NoFacePeakList_SC[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,27] <- sd(NoFacePeakList_SC[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,28] <- mean(abs(NoFacePeakList_SC[[i]]$maxIndex/(2*SamplingRate/tInc_SC)), na.rm = T)
  SynchronyMeasures_matrix[i,29] <- sd(abs(NoFacePeakList_SC[[i]]$maxIndex/(2*SamplingRate/tInc_SC)), na.rm = T)
  
  SynchronyMeasures_matrix[i,30] <- mean(NoFacePeakList_Corr[[i]]$maxValue, na.rm = T) # 2 NA's = geen max WCC gevonden 
  SynchronyMeasures_matrix[i,31] <- sd(NoFacePeakList_Corr[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,32] <- mean(abs(NoFacePeakList_Corr[[i]]$maxIndex/(2*SamplingRate/tInc_EMG)), na.rm = T)
  SynchronyMeasures_matrix[i,33] <- sd(abs(NoFacePeakList_Corr[[i]]$maxIndex/(2*SamplingRate/tInc_EMG)), na.rm = T)
  
  SynchronyMeasures_matrix[i,34] <- mean(NoFacePeakList_Zyg[[i]]$maxValue, na.rm = T) # 2 NA's = geen max WCC gevonden 
  SynchronyMeasures_matrix[i,35] <- sd(NoFacePeakList_Zyg[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,36] <- mean(abs(NoFacePeakList_Zyg[[i]]$maxIndex/(2*SamplingRate/tInc_EMG)), na.rm = T)
  SynchronyMeasures_matrix[i,37] <- sd(abs(NoFacePeakList_Zyg[[i]]$maxIndex/(2*SamplingRate/tInc_EMG)), na.rm = T)
  
  SynchronyMeasures_matrix[i,38] <- mean(NoFacePeakList_SKT[[i]]$maxValue, na.rm = T) # 2 NA's = geen max WCC gevonden 
  SynchronyMeasures_matrix[i,39] <- sd(NoFacePeakList_SKT[[i]]$maxValue, na.rm = T)
  SynchronyMeasures_matrix[i,40] <- mean(abs(NoFacePeakList_SKT[[i]]$maxIndex/(2*SamplingRate/tInc_SKT)), na.rm = T)
  SynchronyMeasures_matrix[i,41] <- sd(abs(NoFacePeakList_SKT[[i]]$maxIndex/(2*SamplingRate/tInc_SKT)), na.rm = T)
} ## warnings: there are no data for a dyad --> dyad is excluded from the analyses

d <- as.data.frame(SynchronyMeasures_matrix, colnames = names(SynchronyMeasures_matrix)) # beware that the measures not included in the manuscript are excluded from the data file (see above)

save("d", file = "Synchrony_summary_allmeasures_bigwindow.RData")

# NOTE: d newly defined here + doesn't contain info which is needed below (WCCmean_SC) 

load("DataFile_Analyses_FBehrens_et_al.RData") #datafile has to be loaded but also doesn't contain the info - other script: load("DataFile_FBehrens_et_al.RData")

