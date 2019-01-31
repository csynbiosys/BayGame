files2<-c("Calibration_4","Calibration_5","Calibration_6",
  "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_6", "DynStim_8",
  "DynStim_11", "DynStim_14")

runAllExp <- function (files, modelName, OptimVal){
  
    data_extraction_multiexperiment(files)
    
    # Initial value optimisation
    m <- stan_model(modelName)
    s1 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
    s2 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
    s3 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
    s4 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 50)
    
    while (s1$value < OptimVal || s2$value < OptimVal || s3$value < OptimVal || s4$value < OptimVal){
      if (s1$value < OptimVal){
        print("s1")
        s1 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
      }
      else if (s2$value < OptimVal){
        print("s2")
        s2 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
      }
      else if (s3$value < OptimVal){
        print("s3")
        s3 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
      }
      else if (s4$value < OptimVal){
        print("s4")
        s4 <- optimizing(m, data = data_multi, algorithm = 'Newton', iter = 200)
      }
      
    }
    
    # Formating initial values to be passed to the stan fit function
    p1 <- list()
    p2 <- list()
    p3 <- list()
    p4 <- list()
    
    q <- c()
    if (modelName == "Model1.stan"){
      q <- 1:13
    } else if (modelName == "Model2.stan"){
      q <- 1:16
    } else if (modelName == "Model3.stan"){
      q <- 1:14
    }
    
    for(x in q){
      val1 <- s1$par[[x]]
      nam1 <- names(s1$par[x])
      val2 <- s2$par[[x]]
      nam2 <- names(s2$par[x])
      val3 <- s3$par[[x]]
      nam3 <- names(s3$par[x])
      val4 <- s4$par[[x]]
      nam4 <- names(s4$par[x])
      
      p1[[nam1]] <- val1
      p2[[nam2]] <- val2
      p3[[nam3]] <- val3
      p4[[nam4]] <- val4
    }
    
    p <- list( p1=p1, p2=p2, p3=p3, p4=p4)
    
    # Inference
    
    inter <- paste("fit_ALL_", modelName, ".rds", sep="")
    fit <- stan(file=modelName, data = data_multi, iter = 3000, warmup = 1000, chains = 4, init = p, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
    fit.dso <- new("cxxdso")
    saveRDS(fit, file = inter)
    
  
}
