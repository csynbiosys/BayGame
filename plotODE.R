fakeData("BangBang_1", c(1,1,1,1,1,1,1,1,1,1))

library(RColorBrewer)
cols<-brewer.pal(n=4,name="Set1")

IPTG <- fd[,1]
aTc <- fd[,2]
RFP <- fd[,3]
GFP <- fd[,4]
ts <- seq(1e-9, length(RFP)-1)
plot(ts, IPTG,col=cols[2])
legend("topright",legend="True Value",col=cols[2],pch=rep(c(16,18),each=4),bty="n",ncol=4,cex=0.7,pt.cex=0.7)
legend("topright",legend="Inference Results",col=cols[1],pch=rep(c(16,18),each=4),bty="n",ncol=1,cex=0.7,pt.cex=0.7)

RFPstdf <<- c()
GFPstdf <<- c()
for(x in seq(1,length(ts))){
  
  v1 = (RFP[x])*0.1
  v2 = (GFP[x])*0.1
  sd1 = sqrt(v1)
  sd2 = sqrt(v2)
  RFPstdf <<- c(RFPstdf, sd1)
  GFPstdf <<- c(GFPstdf, sd2)
  
}

# lines(ts[1:500], fd[,4][1:500]-GFPstdf[1:500])
# lines(ts[1:500], fd[,4][1:500]+GFPstdf[1:500])

for(x in 400:450){
  q <- c(past[x,1,11],past[x,1,12],past[x,1,13],past[x,1,14],past[x,1,15],past[x,1,16],past[x,1,17],past[x,1,18],past[x,1,19],past[x,1,20], (1.19e-1*1.170), (2.06*1.170), 31.94, 9.06e-2, 2, 2)
  fakeData2("BangBang_1", q)
  lines(ts, fd2[,1],col=cols[1])
  
  q <- c(past[x,2,11],past[x,2,12],past[x,2,13],past[x,2,14],past[x,2,15],past[x,2,16],past[x,2,17],past[x,2,18],past[x,2,19],past[x,2,20], (1.19e-1*1.170), (2.06*1.170), 31.94, 9.06e-2, 2, 2)
  fakeData2("BangBang_1", q)
  lines(ts, fd2[,1],col=cols[1])
  
  q <- c(past[x,3,11],past[x,3,12],past[x,3,13],past[x,3,14],past[x,3,15],past[x,3,16],past[x,3,17],past[x,3,18],past[x,3,19],past[x,3,20], (1.19e-1*1.170), (2.06*1.170), 31.94, 9.06e-2, 2, 2)
  fakeData2("BangBang_1", q)
  lines(ts, fd2[,1],col=cols[1])
  
  q <- c(past[x,4,11],past[x,4,12],past[x,4,13],past[x,4,14],past[x,4,15],past[x,4,16],past[x,4,17],past[x,4,18],past[x,4,19],past[x,4,20], (1.19e-1*1.170), (2.06*1.170), 31.94, 9.06e-2, 2, 2)
  fakeData2("BangBang_1", q)
  lines(ts, fd2[,1],col=cols[1])
  
  

}

lines(ts, IPTG,col=cols[2])




# pars = c("k_in_IPTG_raw", "k_in_aTc_raw", "k_L_pm0_raw", "k_L_pm_raw", "theta_T_raw", "theta_aTc_raw", "n_aTc_raw", "n_T_raw")
