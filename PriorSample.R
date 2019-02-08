
################# Define Priors as Gaussians (and hyperparameters) to sample from them (independent events)

MU <- matrix(c(-3.2188758248682,-2.30258509299405,-3.50655789731998,2.30258509299405,3.40119738166216,2.30258509299405,2.5,2.5,
               -2.30258509299405,2.22044604925031e-16,3.40119738166216,-2.30258509299405,2.5,2.5),
             nrow = 14,ncol = 1)

sd <- c(1.15129254649702,1.15129254649702,1.15129254649702,1.15129254649702,1.15129254649702,
        1.15129254649702,1.25,1.25,1.15129254649702,1.15129254649702,1.15129254649702,1.15129254649702,1.25,1.25)

vr <- (sd)^2

E <- diag(14)*vr

############### Sampling from prior distributions

library(MASS)
library(truncnorm)
samplesPre <- mvrnorm(10000, MU, E)

samplesPost <- samplesPre
samplesPost[,1:6] <- exp(samplesPost[,1:6])
samplesPost[,9:12] <- exp(samplesPost[,9:12])
n1 <- rtruncnorm(10000, a = 0, mean = 2.5, sd = 1.25)
n2 <- rtruncnorm(10000, a = 0, mean = 2.5, sd = 1.25)
n3 <- rtruncnorm(10000, a = 0, mean = 2.5, sd = 1.25)
n4 <- rtruncnorm(10000, a = 0, mean = 2.5, sd = 1.25)
samplesPost[,7] <- n1
samplesPost[,8] <- n2
samplesPost[,13] <- n3
samplesPost[,14] <- n4

############## Details of experiment to simulate

expose_stan_functions("Model3_Function_FakeData.stan")

fileName <- "DynStim_8"

fileInputs <- paste(fileName, "_Events_Inputs.csv", sep="")

# Extract input data and stored into global variables
inputs <- read.csv(file=fileInputs, header=TRUE, sep=",")
evnT <- c(round(inputs[,1]),round(inputs[1,2]))
u_IPTG<- inputs[,5]
u_aTc<-inputs[,6]
preI<-inputs[1,3]
preA<-inputs[1,4]
inp<-c()
# Inputs
for (x in 1:length(u_IPTG)){
  inp <- c(inp, u_IPTG[x], u_aTc[x])
}
inp <- inp+1e-7

et <- round(inputs[1,2])

# Time series for ON incubation 
toni <- seq(1e-9, 24*60)

# Time series
ts <- seq(1e-9, et, length=(et+1))

fileObservables <- paste(fileName, "_Observables.csv", sep="")
# Extract observables data and stored into global variables
observables <- read.csv(file=fileObservables, header=TRUE, sep=",")
samplingT <- round(observables[,1])
GFPmean <- observables[,2]
GFPstd <- observables[,3]
RFPmean <- observables[,5]
RFPstd <- observables[,6]

ivss <- c(preI, preA, RFPmean[1], GFPmean[2])
pre <- c(preI, preA)


############## Matrices where results for the simulations will be stored

chr <- matrix(data = 0, ncol = 10000, nrow = et+1)

chg <- matrix(data = 0, ncol = 10000, nrow = et+1)

######## Simulations
for(x in 1:1000) {
  
  q <- samplesPost[x,]
  
  fd2 <- solve_coupled_ode(ts, q, 0, 0, evnT, inp, toni,ivss, pre)
  chr[,x] <- fd2[,3]
  chg[,x] <- fd2[,4]
  
}  


############## Simulations by sampling times

chrST <- matrix(data = 0, ncol = 10000, nrow = (et/5)+1)

chgST <- matrix(data = 0, ncol = 10000, nrow = (et/5)+1)


for(q in 1:10000){
  i = 1
  for(y in seq(1, et+1, 5)){
    chgST[i,q] = chg[y,q]
    chrST[i,q] = chr[y,q]
    i = i+1
  }

}



################





############ LogLikelihood function

LL <- function(x, MUE, EE){
  
  logLg <- -log(sqrt(EE[1]))-0.5*log(2*pi)-((x[1]-MUE[1])^2)/2*EE[1]
  
  logLr <- -log(sqrt(EE[2]))-0.5*log(2*pi)-((x[2]-MUE[2])^2)/2*EE[2]
  
  logL <- logLg+logLr
  
  return(logL)
  
}

lik <- function(x, MUE, EE){
  
  a <- ((x[1]-MUE[1])^2)/EE[1]
  b <- ((x[2]-MUE[2])^2)/EE[2]
  
  Lg <- (1/(sqrt(EE[1])*sqrt(EE[2])*2*pi))*exp(-0.5*(a+b))
  
  return(Lg)
}

############### Mutual information function

MutInf <- function(N1, N2, OSolG, OSolR, sampl){
  
  llm <- matrix(data = 0, nrow = N1, ncol = 1)
  llm2 <- matrix(data = 1, nrow = N1, ncol = N2)
  co <- array(data = 0, dim = c(191, N1, N2))
  
  for(sim in 1:N1){
    
    m <- OSolG[,sim]
    s1 <- c()
    m2 <- OSolR[,sim]
    s2 <- c()
    
    
    
    for(ind in 1:length(m)){
      
      ev <- m[ind]+rnorm(1,0,(m[ind]*0.01))
      s1 <- c(s1, ev)
      ev2 <- m2[ind]+rnorm(1,0,(m2[ind]*0.01))
      s2 <- c(s2, ev2)
      
    }
    
    for(tp in 1:length(s1)){
      
      j <- lik(matrix(c(s1[tp], s2[tp]), nrow = 2, ncol = 1), 
              matrix(c(m[tp], m2[tp]), nrow = 2, ncol = 1), 
              matrix(c((s1[tp]-m[tp])^2, (s2[tp]-m2[tp])^2), nrow = 2, ncol = 1))
      
      # print(s1[tp])
      # print(s2[tp])
      # print(m[tp])
      # print(m2[tp])
      llm[sim] <- llm[sim] +log(j)
      

    }
    
    
    
    for(simm in (N1+1):(N1+N2)){
      
      m3 <- OSolG[,simm]
      m4 <- OSolR[,simm]
      g <- c()
      
      for(tp in 1:length(s1)){
       
        
        
        j2 <- lik(matrix(c(s1[tp], s2[tp]), nrow = 2, ncol = 1), 
                 matrix(c(m3[tp], m4[tp]), nrow = 2, ncol = 1), 
                 matrix(c((s1[tp]-m3[tp])^2, (s2[tp]-m4[tp])^2), nrow = 2, ncol = 1))
        g <- c(g, j2)
        
        # llm2[sim, simm-N1] <- llm2[sim, simm-N1]*j2
        
      }
    
    co[,sim, simm-N1] <- g
    
    # llm2[sim, simm-N1] <- prod(g/mg^length(g))
    # gM[sim, simm-N1] <- mg
    
    }
    
  }
  
  mcr <- c()
  comin <- co
  for(t in 1:N1){
    r <- min(co[,t,])
    mcr <- c(mcr, r)
    comin[,t,] <- comin[,t,]/r
  }
  
  comult <- matrix(data = 1, nrow = N1, ncol = N2)
  for(u in 1:N1){
    for(h in 1:N2){
    
    val <- prod(comin[,u,h])
    comult[u,h] <- val
    
    }
  }
  
  
  MI <- c()
  test <- c()
  counter <- 0
  for(w in 1:N1){
    
    inter <- log(sum(comult[w,]*is.finite(comult[w,]), na.rm = TRUE))+191*log(mcr[w])-log(sum(is.finite(comult[w,])))
    res <- llm[w]-inter
    if(is.finite(res)){
      MI <- c(MI, res)
      counter <- counter+1
    }
    
    test <- c(test, sum(is.finite(comult[w,])))
    
    
  }
  
  
  
  return(fMI)
}




















