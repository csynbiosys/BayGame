

pseudoData <- function(fileName, fit, fakes){
  
  ############## Theta ###############
  parN <- c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]","theta[6]","theta[7]",
            "theta[8]","theta[9]","theta[10]","theta[11]","theta[12]","theta[13]","theta[14]")
  
  q <- readRDS(paste(fit, ".rds", sep=""))
  
  fit_summary <- summary(q, pars = parN, probs = c(0.1, 0.5, 0.9))$summary
  
  median_fit <- fit_summary[, c("50%")]
  
  p <- c()
  
  for (x in 1:14){p <- c(p, median_fit[[x]])}
  
  
  ################ Pseudo-Data ################
  
  expose_stan_functions("Model3_Function_FakeData.stan")
  
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
  
  # Time series
  ts <<- seq(1e-9, et, length=(et+1))
  
  fileObservables <- paste(fileName, "_Observables.csv", sep="")
  # Extract observables data and stored into global variables
  observables <<- read.csv(file=fileObservables, header=TRUE, sep=",")
  samplingT <<- round(observables[,1])
  GFPmean <<- observables[,2]
  GFPstd <<- observables[,3]
  RFPmean <<- observables[,5]
  RFPstd <<- observables[,6]
  
  ivss <- c(preI, preA, RFPmean[1], GFPmean[1])
  pre <- c(preI, preA)
  
  fd <<- solve_coupled_ode(ts, p, 0, 0, evnT, inp, toni,ivss, pre)

  RFP = fd[,3]
  GFP = fd[,4]
  
  # }
  RFP2 <<- c()
  GFP2 <<- c()
  RFPstdf <<- c()
  GFPstdf <<- c()
  ts2 <- c()
  for(x in seq(1,length(ts),5)){
    RFP2 <<- c(RFP2, RFP[x])
    GFP2 <<- c(GFP2, GFP[x])
    v1 = (RFP[x])+rnorm(1,0,(RFP[x]*0.05))
    v2 = (GFP[x])+rnorm(1,0,(GFP[x]*0.05))
    sd1 = abs(RFP[x]-v1)
    sd2 = abs(GFP[x]-v2)
    RFPstdf <<- c(RFPstdf, sd1)
    GFPstdf <<- c(GFPstdf, sd2)
    ts2 <- c(ts2, ts[x])
  }
  
  
  FO <- matrix(data=0, nrow=trunc(length(ts)/5)+1, ncol=6)
  FO[,1] = ts2
  FO[,2] = GFP2
  FO[,3] = GFPstdf
  FO[,4] = ts2
  FO[,5] = RFP2
  FO[,6] = RFPstdf
  
  
  fileObservablesOut <- paste(fileName,"_", fakes, "_Observables.csv", sep="")
  
  write.table(FO, file = fileObservablesOut, row.names = FALSE, col.names = c("timeGFP", "GFPmean", "GFPstd", "timeRFP", "RFPmean", "RFPstd"), sep=",")
  
  fileInputsOut <- paste(fileName, "_", fakes, "_Events_Inputs.csv", sep="")
  write.table(inputs, file = fileInputsOut, row.names = FALSE, sep=",")
  
  fileInputsT <- paste(fileName, "_Inputs.csv", sep="")
  inputsT <- read.csv(file=fileInputsT, header=TRUE, sep=",")
  fileInputsOutT <- paste(fileName, "_", fakes, "_Inputs.csv", sep="")
  
  write.table(inputsT, file = fileInputsOutT, row.names = FALSE, sep=",")
  
  
}



