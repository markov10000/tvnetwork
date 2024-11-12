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

#if (!require(fda)) install.packages('fda')
library(mgm)
library(flare)
library(glmnet)


dofun<- function(seed){

discard <- 0
clime_trim <- 0

lags <- c(1)
maxlag <- max(lags)


  crinet <- array(rep(0,4*8), dim=c(4,8))
  dimnames(crinet)[[2]]<-c('FP','FN','TPR','TNR','PPV','NPV','F1 score','MCC')
  RMSE <- array(rep(0, 2*6), dim=c(2,6))
  ##estimation
  set.seed(seed)
  source("code/var.roots.R",local=TRUE)
  #source("code/kernels.R",local=TRUE)
  source("code/simTOP.R",local=TRUE)
  n <- nrow(X)
  p <- ncol(X)              

  estpoints <- seq(0, 1, length = n-maxlag+1)[2:(n-maxlag+1)]#
  
  #Estimate Cnet
  eps <- 1e-5
  adj_GIC <- 0.3
  lambda2.fac <- 1 / (1 + seq(0, 20, 3)) #/ (3.7+1)#20 3
  mvarmodel<- mvar(X,type=rep("g",d),lags=lags)
  mvarmodel2<- rewLASSO(mvarmodel,X)[[1]]
  #rm(comb)

  constantCnet1<- Matrix(mvarmodel$wadj[,,1]!=0,nrow=d,ncol=d, sparse=TRUE)
  constantCnet2<- Matrix(mvarmodel2!=0,nrow=d,ncol=d, sparse=TRUE)
  
 
  #Summary
  ##constantCnet1
  TP <- sum(TrueCnet*constantCnet1)
  FP <- sum(constantCnet1)-TP
  FN <- sum(TrueCnet)-TP
  TN <- d^2 - TP - FP - FN
  TPR<- TP/(TP+FN);TNR<- TN/(TN+FP);PPV<-TP/(TP+FP);NPV<-TN/(TN+FN)
  crinet[1,1] <- FP
  crinet[1,2] <- FN
  crinet[1,3] <- TPR
  crinet[1,4] <- TNR
  crinet[1,5] <- PPV
  crinet[1,6] <- NPV
  crinet[1,7] <- (2*TP)/(2*TP+FP+FN) # F1 score  
  #Matthews correlation coefficient (MCC)
  crinet[1,8] <- sqrt(TPR)*sqrt(TNR)*sqrt(PPV)*sqrt(NPV)-sqrt(1-TPR)*sqrt(1-TNR)*sqrt(1-PPV)*sqrt(1-NPV)

  ##constantCnet
  TP <- sum(TrueCnet*constantCnet2)
  FP <- sum(constantCnet2)-TP
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
  estpoints_res <- c(1:(n-1)) / (n-1)
  coef1 <- array(rep(mvarmodel$wadj[,,1],n-1),dim=c(p,p,1,n-1))
  coef2 <- array(rep(mvarmodel2,n-1),dim=c(p,p,1,n-1))
  coef <- list(coef1,coef2)
  residuals <- get.residuals(data=X, estpoints=estpoints_res, lags=c(1:1), coef=coef)
  residuals1 <-residuals[[1]]
  residuals2 <-residuals[[2]]
  perr <- get.predict.error(X, lags, discard=discard, coef=coef)
  error1step1 <- perr[[1]]
  error1step2 <- perr[[2]]
  
  listC1 <- lapply(seq(dim(coef1)[4]), function(x) coef1[ , ,1, x])
  SST1 <- apply(X[-1,],2,var)* (n-1)  
  RMSE[1,1] <- 1-mean(apply((residuals1)^2,2,sum)/SST1)
  RMSE[1,2] <- mean(mapply( function(x,y)norm(x-y,"F"),listC1,Ai))/sqrt(d) #RMSE
  RMSE[1,3] <- sqrt(mean((residuals1-epsilon[-1,])^2))
  RMSE[1,4] <- sqrt(mean((error1step1)^2))
  
  listC2 <- lapply(seq(dim(coef2)[4]), function(x) coef2[ , ,1, x])
  SST2 <- apply(X[-1,],2,var)* (n-1)  
  RMSE[2,1] <- 1-mean(apply((residuals2)^2,2,sum)/SST2)
  RMSE[2,2] <- mean(mapply( function(x,y)norm(x-y,"F"),listC2,Ai))/sqrt(d) #RMSE
  RMSE[2,3] <- sqrt(mean((residuals2-epsilon[-1,])^2))
  RMSE[2,4] <- sqrt(mean((error1step2)^2))
  return(list(crinet=crinet,RMSE=RMSE))    
  # lambda=netCLIME$lambda_seq
}

results <- dofun(run_number)
saveRDS(results, sprintf("results/run_%d.rds", run_number))

timeend<-Sys.time()
print(timeend-timestart)





