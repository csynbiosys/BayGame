files<-c("Calibration_1", "Calibration_2","Calibration_3","Calibration_4","Calibration_5","Calibration_6", "BangBang_1", "BangBang_2",
         "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_4", "DynStim_5", "DynStim_6", "DynStim_7", "DynStim_8", "DynStim_9", "DynStim_10", 
         "DynStim_11", "DynStim_12", "DynStim_13", "DynStim_14", "DynStimTooFast", "PI_1", "PI_2", "PI_3")


data_extraction_multiexperiment <- function (fileNamesVector){
  
  ######################## OBSERVABLES ########################
  
  mdp <- c() # Maximum data point for each file
  for(i in fileNamesVector){ # Loop that opens a file for each iteration to extract the file with maximum data points
    fileObs <- paste(i, "_Observables.csv", sep="")
    obs <- read.csv(file=fileObs, header=TRUE, sep=",")
    mdp = c(mdp, length(obs[,1]))
  }
  mrow <- max(mdp) # Maximum number of rows
  mcol <- length(fileNamesVector) # Maximum number of columns
  
  # Definition of the matrices with 0s on them
  samplingT <- matrix(data=0, nrow=mrow, ncol=mcol)
  GFPmean <- matrix(data=0, nrow=mrow, ncol=mcol)
  GFPstd <- matrix(data=0, nrow=mrow, ncol=mcol)
  RFPmean <- matrix(data=0, nrow=mrow, ncol=mcol)
  RFPstd <- matrix(data=0, nrow=mrow, ncol=mcol)
  
  for(i in 1:length(fileNamesVector)){ # Loop that opens a file for each iteration to extract the data
    fileObs <- paste(fileNamesVector[i], "_Observables.csv", sep="")
    obs <- read.csv(file=fileObs, header=TRUE, sep=",")
      for(j in 1:length(obs[,1])){
        samplingT[j,i] <- round(obs[j,1])
        GFPmean[j,i] <- obs[j,2]
        GFPstd[j,i] <- obs[j,3]
        RFPmean[j,i] <- obs[j,5]
        RFPstd[j,i] <- obs[j,6]
      }
  }
  
  ######################## INPUTS ########################
  
  mdpI <- c() # Maximum data point for each file
  mtimes <- c() # Maximum time point for each file
  for(i in fileNamesVector){ # Loop that opens a file for each iteration to extract the file with maximum data points
    fileInp <- paste(i, "_Events_Inputs.csv", sep="")
    inp <- read.csv(file=fileInp, header=TRUE, sep=",")
    mdpI = c(mdpI, length(inp[,1]))
    mtimes = c(mtimes, round(inp[1,2]))
  }
  mrowI <- max(mdpI) # Maximum number of rows
  mcolI <- length(fileNamesVector) # Maximum number of columns
  mt <- max(mtimes) # Maximum number of time points

  evnT <- matrix(data=0, nrow=mrowI+1, ncol=mcolI)
  u_IPTG <- matrix(data=0, nrow=mrowI, ncol=mcolI)
  u_aTc <- matrix(data=0, nrow=mrowI, ncol=mcolI)
  preI <- matrix(data=0, nrow=1, ncol=mcolI)#<<-inputs[1,3]
  preA <- matrix(data=0, nrow=1, ncol=mcolI)#<<-inputs[1,4]
  time <- matrix(data=0, nrow=mt+1, ncol=mcolI)#<<- seq(1e-9, round(inputs[1,2]), length=round(inputs[1,2]))
  inps <- matrix(data=0, nrow=mrowI*2, ncol=mcolI)#<<-c()
  
  for(i in 1:length(fileNamesVector)){ # Loop that opens a file for each iteration to extract the data
    fileInp <- paste(fileNamesVector[i], "_Events_Inputs.csv", sep="")
    inp <- read.csv(file=fileInp, header=TRUE, sep=",")
    tempT <- seq(1e-9, round(inp[1,2]), length=round(inp[1,2])+1)
    ind = 1
    for(j in 1:length(inp[,1])){
      evnT[j,i] <- round(inp[j,1])
      u_IPTG[j,i] <- inp[j,5]
      u_aTc[j,i] <- inp[j,6]
      inps[ind,i] <- inp[j,5]
      inps[ind+1,i] <- inp[j,6]
      ind = ind+2
    }
    
    for(l in 1:length(tempT)){
      
      time[l,i] = tempT[l]
    }
    evnT[(length(inp[,1])+1),i] = round(inp[1,2])
    preI[,i] <- inp[1,3]
    preA[,i] <- inp[1,4]
     
  }

    toni <- seq(1e-9, 24*60) # Over night incuvation time
  
    data_multi <<- list (
          
          elm = mrowI, # Maximum length of the rows of the matrices except for time and evnT and pres
          tml = trunc(mt+1), # Maximum length of the rows for the time matrix
          ts = time+(1e-9),
          tsl = round(mtimes), # length of time series per event
          tsmax = round(mtimes)+1, # maximum time per event
          preIPTG = preI,
          preaTc = preA,
          IPTG = u_IPTG,
          aTc = u_aTc,
          Nsp = (mdpI+1), # length(evnT),
          inputs = inps,
          evnT = round(evnT),
          m = mcol, # Number of time series
          stsl = mdp, # Number of elements at each time series
          stslm = mrow,
          sts = trunc(samplingT),
          GFPmean = GFPmean,
          RFPmean = RFPmean,
          GFPstd = GFPstd,
          RFPstd = RFPstd,
          toni = toni,
          tonil = length(toni)

  )
}






