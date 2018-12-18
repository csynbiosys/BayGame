
files<-c(
  # "Calibration_1", "Calibration_2","Calibration_3",
  "Calibration_4","Calibration_5","Calibration_6", "BangBang_1", "BangBang_2",
  "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_4", "DynStim_5", "DynStim_6", "DynStim_7", "DynStim_8", "DynStim_9", "DynStim_10", 
  "DynStim_11", "DynStim_12", "DynStim_13", "DynStim_14", "DynStimTooFast", "PI_1", "PI_2", "PI_3")


runExp <- function (experiments){
  
  for(x in experiments){
    data_extraction_multiexperiment(x)
    inter <- paste("fit_", x, ".rds", sep="")
    fit <- stan(file="Model3.stan", data = data_multi, iter = 3000, warmup = 1000, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
    fit.dso <- new("cxxdso")
    saveRDS(fit, file = inter)
    
    mes <- paste("Fit for ", x , " has finished!", sep = "")
    print(mes)
  }
  
}