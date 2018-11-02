
# To be able to use Stan Functions in R
expose_stan_functions("Model2_Function.stan")

# Inputs for the function
ts <- seq(1e-9, 60*24, length=(60*24)+1)
y0ss <- c(0,0,215.1551459, 782.4044739)
p <- c(2.75e-2, 1.11e-1, 1.62e-1, 2e-2, (3.2e-2*9.726e-1), (8.3*9.726e-1), 30, 11.65, 2, 2, (1.19e-1*1.170), (2.06*1.170), 31.94, 9.06e-2, 2, 2)
iss <- c(0,0)
x_i <- 0

# ODE input values
IPTG <- c(0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10)
aTc <- c(0,20,40,60,80,100)

# Calculation of ODE solutions and extraction of the last value at steady state
results1RFP <- matrix(data=NA,nrow=length(aTc),ncol=length(IPTG))
results1GFP <- matrix(data=NA,nrow=length(aTc),ncol=length(IPTG))

for(x in 1:length(IPTG)){
  for (y in 1:length(aTc)){
    inp <- c(IPTG[x],aTc[y])
    
    tps <- solve_coupled_ode(ts, y0ss, p, iss, inp, 0,0)
    
    results1RFP[y,x] = tps[length(tps[,3]),3]
    results1GFP[y,x] = tps[length(tps[,4]),4]
  
  }
}


aTc2 <- c(0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10)
IPTG2 <- c(0,0.20,0.40,0.60,0.80,1.00)

results2RFP <- matrix(data=NA,nrow=length(IPTG2),ncol=length(aTc2))
results2GFP <- matrix(data=NA,nrow=length(IPTG2),ncol=length(aTc2))



for(x in 1:length(aTc2)){
  
  for(y in 1:length(IPTG2)){
    inp <- c(IPTG2[y],aTc2[x])
    tps <- solve_coupled_ode(ts, y0ss, p, iss, inp, 0,0)
    results2RFP[y,x] = tps[length(tps[,3]),3]
    results2GFP[y,x] = tps[length(tps[,4]),4]
  }
  
}


# Write CSV files
write.csv(results1RFP, file = "RFP_IPTG_0.csv")
write.csv(results1GFP, file = "GFP_IPTG_0.csv")
write.csv(results2RFP, file = "RFP_aTc_0.csv")
write.csv(results2GFP, file = "GFP_aTc_0.csv")



# Calculation of ratios
ratio1RFP <- matrix(data=NA,nrow=length(aTc),ncol=length(IPTG))
ratio1GFP <- matrix(data=NA,nrow=length(aTc),ncol=length(IPTG))
ratio2RFP <- matrix(data=NA,nrow=length(IPTG2),ncol=length(aTc2))
ratio2GFP <- matrix(data=NA,nrow=length(IPTG2),ncol=length(aTc2))


for(x in 1:length(aTc)){
  r1 = results1RFP[x,]/results1RFP[x,1]
  r2 = results1GFP[x,]/results1GFP[x,1]
  r3 = results2RFP[x,]/results2RFP[x,1]
  r4 = results2GFP[x,]/results2GFP[x,1]
  
  ratio1RFP[x,] = r1
  ratio1GFP[x,] = r2
  ratio2RFP[x,] = r3
  ratio2GFP[x,] = r4
  
  
}


# Write CSV files
write.csv(ratio1RFP, file = "RFP_IPTG_0R.csv")
write.csv(ratio1GFP, file = "GFP_IPTG_0R.csv")
write.csv(ratio2RFP, file = "RFP_aTc_0R.csv")
write.csv(ratio2GFP, file = "GFP_aTc_0R.csv")


# Calculation of error in respect
ratio1RFPerror <- 1-ratio1RFP
ratio1GFPerror <- 1-ratio1GFP
ratio2RFPerror <- 1-ratio2RFP
ratio2GFPerror <- 1-ratio2GFP



















