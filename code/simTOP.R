library(Matrix)
library(MASS)
n = 400 #400 200
d = 100  #100 50

constant3 <- 0.4
constant4 <- 0.8
constant1 <- 0.1
constant2 <- 1
constant5 <- 1
epsilonscale<- sqrt(1)


Ai0 <- toeplitz((constant3)^(1+constant2*(0:(d-1))))
Ai20 <- toeplitz((constant4)^(0:(d-1)))

#if(!exists('rpti')) rpti<-1
#set.seed(rpti)
set.seed(seed)
trueCOV<-trueiCOV<-list() #t=2 to T, index from 1
Ai <- list()              #t=2 to T, index from 1
  epsilon <- matrix(0,nrow=n,ncol=d)
trueiCOV0 <- Ai20
trueCOV0 <- solve(trueiCOV0)
epsilon0 <- mvrnorm(n=1,mu=rep(0,d),Sigma=trueCOV0)*epsilonscale
epsilon[1,] <- mvrnorm(n=1,mu=rep(0,d),Sigma=trueCOV0)*epsilonscale
X <- matrix(0, nrow=n, ncol=d)




maxeigen <- 0
mineigen <- 1
X[1,] <- constant5 * as.numeric(Ai0 %*% epsilon0)  + epsilon[1,] 
for (t in (2:n))
{
  Ai[[t-1]] <- toeplitz((constant3-constant1*(t/n))^(1+constant2*(0:(d-1))))
  Ai2 <- toeplitz(((constant4)-constant1*(t/n))^(0:(d-1)))
  trueiCOV[[t-1]] <- Ai2
  trueCOV[[t-1]] <- chol2inv(chol(Ai2)) #sym inverse compare to solve
  epsilon[t,] <- mvrnorm(n=1,mu=rep(0,d),Sigma=trueCOV[[t-1]])*epsilonscale 
  X[t,] <- constant5 * as.numeric(Ai[[t-1]]%*%X[t-1,]) + epsilon[t,]
  maxeigen <- max(maxeigen,abs(var.roots(list(as.matrix(Ai[[t-1]]*constant5)))))
  maxeigen <- max(maxeigen,norm(Ai[[t-1]]*constant5,'2'))
  mineigen <- min(mineigen, eigen(trueCOV[[t-1]])$values)
}

TrueCnet<-(Ai0!=0)*1
TrueAiL2norm <- as.matrix(sqrt(Reduce("+", lapply(Ai,function(x)x^2))))
TruePnet<-(Ai20!=0)*1
diag(TruePnet) <- rep(0,d)
#heatmap(A01,Rowv=NA,Colv=NA)



#A <- matrix(c(0.7,0.7,0,0,0.7,0.7,0, 0,0.7),3,3)
#eigen(A,symmetric = FALSE)$values
#svd(A)$d
