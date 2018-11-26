
log_lik_postProc <- function (fit, fileNamesVector, chains){
  
  # Oppen Loo package
  library(loo)
  
  # Extract log_likelihoods from the stan fit object
  logLikelihood_RFP_pre <<- extract_log_lik(fit, "logLikelihood_RFP")
  logLikelihood_GFP_pre <<- extract_log_lik(fit, "logLikelihood_GFP")
  
  
  mdp <- c() # Maximum data point for each file
  for(i in fileNamesVector){ # Loop that opens a file for each iteration to extract the file with maximum data points
    fileObs <- paste(i, "_Observables.csv", sep="")
    obs <- read.csv(file=fileObs, header=TRUE, sep=",")
    mdp = c(mdp, length(obs[,1]))
  }
  mrow <- max(mdp) # Maximum number of rows
  mcol <- length(fileNamesVector) # Maximum number of columns

  # Empty matrices that will store the log_likelihood valus without the introduced 0
  logLikelihood_RFP <<- matrix(data=0, nrow=length(logLikelihood_RFP_pre[,1]), ncol= sum(mdp))
  logLikelihood_GFP <<- matrix(data=0, nrow=length(logLikelihood_GFP_pre[,1]), ncol= sum(mdp))
  
  # Ad an initial point 1 to the lengths of each experiment sampling
  mdp <- c(1,mdp)
  
  # Exstraction of the values of interest and introduction to the empti matrices
  for(y in 1:length(logLikelihood_RFP_pre[,1])){
    i <- 1
    j <- mrow
    w <- 1
    e <- mdp[2]
    
    for(x in 1:(length(mdp)-1)){

      q1 <- logLikelihood_RFP_pre[y,][i:j]
      q2 <- logLikelihood_GFP_pre[y,][i:j]
      
      logLikelihood_RFP[y,][w:e] <<- q1[1:mdp[x+1]]
      logLikelihood_GFP[y,][w:e] <<- q2[1:mdp[x+1]]

      i = i+mrow
      j = j+mrow
      w = w+mdp[x+1]
      e = mdp[x+2]+e
      
    }
  }
  
  ##################################### DIVIDE BY EXPERIMENTS #####################################
  
  w <- 1
  e <- mdp[2]

  for(x in 1:(length(mdp)-1)){
    
    # Rows of interest (all)
    q1 <- logLikelihood_RFP[1:length(logLikelihood_RFP_pre[,1]),]
    q2 <- logLikelihood_GFP[1:length(logLikelihood_RFP_pre[,1]),]
    
    # Each iteration is gona be a matrix for with the log-likelihood for each of the experiments considered
    a1 <- q1[,w:e]
    a2 <- q2[,w:e]

    w = w+mdp[x+1]
    e = mdp[x+2]+e
    
  }

  ##################################### DIVIDE BY CHAINS #####################################
  
  # Iterations per chain
  ipc <- length(logLikelihood_RFP_pre[,1])/chains

  m <- 1
  n <- ipc

  for(x in 1:chains){
    
    # Rows of interest (each chain separately)
    q1 <- logLikelihood_RFP[m:n,]
    q2 <- logLikelihood_GFP[m:n,]
    
    m = m+ipc
    n = n+ipc

  }
  
  ##################################### DIVIDE BY BOTH #####################################
  
  m <- 1
  n <- ipc
  
  
  for(y in 1:chains){
    w <- 1
    e <- mdp[2]
    
    for(x in 1:(length(mdp)-1)){
      
      # Rows of interest (all)
      q1 <- logLikelihood_RFP[m:n,]
      q2 <- logLikelihood_GFP[m:n,]
      
      # Each iteration is gona be a matrix for with the log-likelihood for each of the experiments considered
      a1 <- q1[,w:e]
      a2 <- q2[,w:e]
      print(a1)
      
      w = w+mdp[x+1]
      e = mdp[x+2]+e
    }
    m = m+ipc
    n = n+ipc
  }
  
  
}
