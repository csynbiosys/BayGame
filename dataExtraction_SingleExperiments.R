
# Function to extract all the data from the CSV files and and stored in a list that will be used as a data input for a RStan script

data_extraction <- function (fileName){
  
  # Set the name of the two (inputs and observables) files just by introducing the tag of the files
  fileInputs <- paste(fileName, "_Events_Inputs.csv", sep="")
  fileObservables <- paste(fileName, "_Observables.csv", sep="")
  fileInputsT <- paste(fileName, "_Inputs.csv", sep="")
  
  # Extract input data and stored into global variables
  inputs <<- read.csv(file=fileInputs, header=TRUE, sep=",")
  evnT <<- c(round(inputs[,1]),round(inputs[1,2]))
  u_IPTG<<- inputs[,5]
  u_aTc<<-inputs[,6]
  preI<<-inputs[1,3]
  preA<<-inputs[1,4]
  inp<<-c()
  # Inputs
  for (x in 1:length(u_IPTG)){
    inp <<- c(inp, u_IPTG[x], u_aTc[x])
  }
  time<- seq(1e-9, round(inputs[1,2]), length=round(inputs[1,2])+1)
  
  inputsT <- read.csv(file=fileInputsT, header=TRUE, sep=",")
  
  time2 <- inputsT[,1]
  time2[1] = 1e-9
  
  # Extract observables data and stored into global variables
  observables <<- read.csv(file=fileObservables, header=TRUE, sep=",")
  samplingT <<- round(observables[,1])
  GFPmean <<- observables[,2]
  GFPstd <<- observables[,3]
  RFPmean <<- observables[,5]
  RFPstd <<- observables[,6]
  
  # Time series for ON incubation 
  toni <- seq(1e-9, 24*60)
  
  # Store all data and modified data into a list for RSTan
  data_real <<- list (
 
                     ts = time2, # Total time vector
                     tsl = length(time2), # Length of time vector
                     tsmax = time2[length(time2)], # Last element of the time vector
                     time0 = time2[1], # First element of the time vector
                     preIPTG = preI, # Pres values
                     preaTc = preA,
                     IPTG = u_IPTG, # Input values at each event
                     aTc = u_aTc,
                     Nsp = length(evnT), # Number of event switching points (including initial and final)
                     inputs = inp+1e-5, # Inputs as IPTG, aTc, IPTG, aTc, ...
                     evnT = evnT, # Event switching points (including initial and final)
                     stsl = length(samplingT), # Number of sampling times
                     sts = samplingT, # Sampling times
                     GFPmean = GFPmean, # Mean and standard deviations for proteins
                     RFPmean = RFPmean,
                     GFPstd = GFPstd,
                     RFPstd = RFPstd,
                     toni = toni,
                     tonil = length(toni)
  )
  
}