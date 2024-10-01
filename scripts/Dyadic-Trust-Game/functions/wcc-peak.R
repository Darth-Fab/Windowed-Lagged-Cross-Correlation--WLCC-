## Windowed Lagged Cross-Correlation function and Peak Picking Algorithm
## Friederike Behrens & Steven Boker, 2018

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
