
# Function to extract all the data from the CSV files and and stored in a list that will be used as a data input for a RStan script

data_extraction <- function (fileName){
  
  # Set the name of the two (inputs and observables) files just by introducing the tag of the files
  fileInputs <- paste(fileName, "_Events_Inputs.csv", sep="")
  fileObservables <- paste(fileName, "_Observables.csv", sep="")
  
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
  time<<- seq(1e-9, round(inputs[1,2]), length=round(inputs[1,2])+1)
  
  # Extract observables data and stored into global variables
  observables <<- read.csv(file=fileObservables, header=TRUE, sep=",")
  samplingT <<- observables[,1]
  GFPmean <<- observables[,2]
  GFPstd <<- observables[,3]
  RFPmean <<- observables[,5]
  RFPstd <<- observables[,6]
  
  # Store all data and modified data into a list for RSTan
  data_real <<- list (
                     ts = time,
                     ts2 = round(time),
                     tsl = length(time),
                     tsmax = time[length(time)],
                     time0 = time[1],
                     preIPTG = preI,
                     preaTc = preA,
                     IPTG = u_IPTG,
                     aTc = u_aTc,
                     Nsp = length(evnT),
                     inputs = inp,
                     evnT = evnT,
                     stsl = length(samplingT),
                     sts = trunc(samplingT),
                     GFPmean = GFPmean,
                     RFPmean = RFPmean,
                     GFPstd = GFPstd,
                     RFPstd = RFPstd
                     
  )
  
}