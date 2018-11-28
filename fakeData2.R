fakeData2 <- function(fileName, p2){
  
  expose_stan_functions("Model2_Function_FakeData.stan")
  
  # Set the name of the two (inputs and observables) files just by introducing the tag of the files
  fileInputs <- paste(fileName, "_Events_Inputs.csv", sep="")
  
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
  inp <- inp+1e-7
  
  et <- round(inputs[1,2])
  
  # Time series for ON incubation 
  toni <- seq(1e-9, 24*60)
  
  # Parameters
  p <- c(2.75e-2, 1.11e-1, 1.62e-1, 2e-2, (3.2e-2*9.726e-1), (8.3*9.726e-1), 30, 11.65, 2, 2, (1.19e-1*1.170), (2.06*1.170), 31.94, 9.06e-2, 2, 2)
  # p2 <- c(0.1024555, 0.007232621 , 0.03033617 , 0.3158936 , 0.004402344 , 82.65198, 7.402744, 9.963885 , 4.574192, 4.595242 , (1.19e-1*1.170), (2.06*1.170), 31.94, 9.06e-2, 2, 2)
  
  # Time series
  ts <- seq(1e-9, et, length=(et+1))
  
  fileObservables <- paste(fileName, "_Observables.csv", sep="")
  # Extract observables data and stored into global variables
  observables <<- read.csv(file=fileObservables, header=TRUE, sep=",")
  samplingT <<- round(observables[,1])
  GFPmean <<- observables[,2]
  GFPstd <<- observables[,3]
  RFPmean <<- observables[,5]
  RFPstd <<- observables[,6]
  
  ivss <- c(RFPmean[1], GFPmean[2])
  pre <- c(preI, preA)
  
  # GFP <- matrix(data=0, nrow=length(ts), ncol=10)
  # RFP <- matrix(data=0, nrow=length(ts), ncol=10)
  
  # for(x in 1:10){
  fd <<- solve_coupled_ode(ts, p, 0, 0, evnT, inp, toni,ivss, pre)
  fd2 <<- solve_coupled_ode(ts, p2, 0, 0, evnT, inp, toni,ivss, pre)
  RFP = fd[,3]
  GFP = fd[,4]
  
  # }
  RFP2 <<- c()
  GFP2 <<- c()
  RFPstdf <<- c()
  GFPstdf <<- c()
  ts2 <- c()
  for(x in seq(1,length(ts),10)){
    RFP2 <<- c(RFP2, RFP[x])
    GFP2 = c(GFP2, GFP[x])
    v1 = (RFP[x])*0.1
    v2 = (GFP[x])*0.1
    sd1 = sqrt(v1)
    sd2 = sqrt(v2)
    RFPstdf <<- c(RFPstdf, sd1)
    GFPstdf <<- c(GFPstdf, sd2)
    ts2 <- c(ts2, ts[x])
  }
  
  # GFPstdf <- c()
  # for(x in 1:length(GFP)){
  #   v = (GFP[x])*0.1
  #   sd = sqrt(v)
  #   GFPstdf = c(GFPstdf, sd)
  # }
  
  
  FO <- matrix(data=0, nrow=trunc(length(ts)/10)+1, ncol=6)
  FO[,1] = ts2
  FO[,2] = GFP2
  FO[,3] = GFPstdf
  FO[,4] = ts2
  FO[,5] = RFP2
  FO[,6] = RFPstdf
  
  
  fileObservablesOut <- paste(fileName, "_Fake_Observables.csv", sep="")
  
  write.table(FO, file = fileObservablesOut, row.names = FALSE, col.names = c("timeGFP", "GFPmean", "GFPstd", "timeRFP", "RFPmean", "RFPstd"), sep=",")
  
  fileInputsOut <- paste(fileName, "_Fake_Events_Inputs.csv", sep="")
  write.table(inputs, file = fileInputsOut, row.names = FALSE, sep=",")
  
  fileInputsT <- paste(fileName, "_Inputs.csv", sep="")
  inputsT <- read.csv(file=fileInputsT, header=TRUE, sep=",")
  fileInputsOutT <- paste(fileName, "_Fake_Inputs.csv", sep="")
  
  write.table(inputsT, file = fileInputsOutT, row.names = FALSE, sep=",")
  
  
}