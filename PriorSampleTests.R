N1 <- 1
N2 <- 1
OSolG <- chgST
OSolR <- chrST


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
  
  
  
  foreach(simm = (N1+1):(N1+N2)) %dopar%{
    
    m3 <- OSolG[,simm]
    m4 <- OSolR[,simm]
    g <- c()
    
    for(tp in 1:length(s1)){
      
      
      
      j2 <- lik(matrix(c(s1[tp], s2[tp]), nrow = 2, ncol = 1), 
                matrix(c(m[tp], m2[tp]), nrow = 2, ncol = 1), 
                matrix(c((s1[tp]-m[tp])^2, (s2[tp]-m2[tp])^2), nrow = 2, ncol = 1))
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
    # print(val)
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
  
  # print(log(sum(comult[w,]*is.finite(comult[w,]), na.rm = TRUE)))
}


fMI <- sum(MI)/counter










