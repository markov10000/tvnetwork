info=version
cat("start...")
timestart<-Sys.time()
print(timestart)
getwd()

# Get the run number
args <- commandArgs(trailingOnly=TRUE)
run_number <- as.integer(args[1])

source("code/SCAD.derivative.R",local=TRUE)
source("code/lagData.R",local=TRUE)
source("code/rewLASSO.R",local=TRUE)
source("code/get.residuals.R",local=TRUE)
source("code/get.predict.error.R",local=TRUE)
source("code/kernels.R",local=TRUE)
source("code/localpca.R",local=TRUE)
source("code/tvvarnets.R",local=TRUE)
source("code/tvCOV.R",local=TRUE)
source("code/calcLL.R",local=TRUE)

#Sys.setenv(PATH = paste("C:/Rtools40/bin", Sys.getenv("PATH"), sep=";"))
#Sys.setenv(BINPREF = "C:/Rtools40/mingw64/bin/")
#system("R CMD SHLIB ~/Desktop/network/R/rawfit_wglasso.c")
#dyn.load("rawfit_wglasso.dll")
#dyn.unload("rawfit_wglasso.dll")

library(mgm)
library(flare)
library(glmnet)
library(data.table)
library(parallel)
#library(plyr)
library(doParallel)
library(foreach)

#setup parallel backend to use many processors
cores=detectCores()
print(cores[1])
#cl <- makeCluster(1) #cores[1]-1 not to overload your computer
#registerDoParallel(cl)

dofun<- function(seed){
discard <- 0
clime_trim <- 0
lags <- c(1)
maxlag <- max(lags)


  dyn.load("code/rawfit_wglasso.so", local=TRUE)
  crinet <- array(rep(0,4*8), dim=c(4,8))
  dimnames(crinet)[[2]]<-c('FP','FN','TPR','TNR','PPV','NPV','F1 score','MCC')
  RMSE <- array(rep(0, 6), dim=c(6))
  ##estimation
  set.seed(seed)
  source("code/var.roots.R",local=TRUE)
  source("code/simUPFac.R",local=TRUE)
                
  bandwidth <- 0.75*(log(d)/n)^0.2 #LKZ15,AOS
  #estpoints <- c(0,0.5,1)
  estpoints <- seq(0, 1, length = n-lags) 
  #estpoints <- seq(0, 1, length = n-maxlag+1)[2:(n-maxlag+1)]#
  
  #Estimate Cnet
  # hat.X is estimate of X using localpca
  mvarmodel<- mvar(X,type=rep("g",d),lags=lags)
  netmodel<- tvvarnets(X,estpoints=estpoints,discard=discard, bandwidth = bandwidth, lags=c(1), lambdaSel='BIC')
  
  #rm(comb)
  
  netCnet<- Matrix(rep(FALSE,d*d),nrow=d,ncol=d, sparse=TRUE)
  constantCnet<- Matrix(mvarmodel$wadj[,,1]!=0,nrow=d,ncol=d, sparse=TRUE)
  
  for (i in 1:length(estpoints)){
    netCnet <- netCnet | (netmodel$glasso_wadj[,,1,i]!=0)
    #netPnet <- netPnet | (abs(netCLIME$iCOVoutput[[i]])>clime_trim)
  }
  #Summary
  ##netCnet
  TP <- sum(TrueCnet*netCnet)
  FP <- sum(netCnet)-TP
  FN <- sum(TrueCnet)-TP
  TN <- d^2 - TP - FP - FN
  TPR<- TP/(TP+FN);TNR<- TN/(TN+FP);PPV<-TP/(TP+FP);NPV<-TN/(TN+FN)
  crinet[1,1] <- FP
  crinet[1,2] <- FN
  crinet[1,3] <- TPR
  crinet[1,4] <- TNR
  crinet[1,5] <- PPV
  crinet[1,6] <- NPV
  crinet[1,7]  <- (2*TP)/(2*TP+FP+FN) # F1 score    
  #Matthews correlation coefficient (MCC)
  crinet[1,8] <- sqrt(TPR)*sqrt(TNR)*sqrt(PPV)*sqrt(NPV)-sqrt(1-TPR)*sqrt(1-TNR)*sqrt(1-PPV)*sqrt(1-NPV)
  
  ##constantCnet
  TP <- sum(TrueCnet*constantCnet)
  FP <- sum(constantCnet)-TP
  FN <- sum(TrueCnet)-TP
  TN <- d^2 - TP - FP - FN
  TPR<- TP/(TP+FN);TNR<- TN/(TN+FP);PPV<-TP/(TP+FP);NPV<-TN/(TN+FN)
  crinet[2,1] <- FP
  crinet[2,2] <- FN
  crinet[2,3] <- TPR
  crinet[2,4] <- TNR
  crinet[2,5] <- PPV
  crinet[2,6] <- NPV
  crinet[2,7] <- (2*TP)/(2*TP+FP+FN) # F1 score  
  #Matthews correlation coefficient (MCC)
  crinet[2,8] <- sqrt(TPR)*sqrt(TNR)*sqrt(PPV)*sqrt(NPV)-sqrt(1-TPR)*sqrt(1-TNR)*sqrt(1-PPV)*sqrt(1-NPV)

  #criteria R2, RMSE
  estpoints_res <- c(1:(n-maxlag)) / (n-maxlag)
  residuals <- get.residuals(data=X, estpoints=estpoints_res, lags=c(1:maxlag), coef=netmodel$glasso_wadj)
  residuals <-residuals[[1]]
  
  listC <- lapply(seq(dim(netmodel$glasso_wadj)[4]), function(x) netmodel$glasso_wadj[ , ,1, x])
  error1step <- get.predict.error(X, lags, discard=discard, coef=netmodel$glasso_wadj)[[1]]
  SST <- apply(X[-1,],2,var)* (n-1)  
  RMSE[1] <- 1-mean(apply((residuals)^2,2,sum)/SST)
  RMSE[2] <- mean(mapply( function(x,y)norm(x-y,"F"),listC,Ai))/sqrt(d) #RMSE
  RMSE[3] <- sqrt(mean((residuals-epsilon[-1,])^2))
  RMSE[4] <- sqrt(mean((error1step)^2))
  return(list(crinet=crinet,RMSE=RMSE, haterr=residuals))    
  # lambda=netCLIME$lambda_seq
}

results <- dofun(run_number)
saveRDS(results, sprintf("results/run_%d.rds", run_number))

#stop cluster
#stopCluster(cl)
timeend<-Sys.time()
print(timeend-timestart)

