library(Matrix)
library(MASS)
n = 200 #400 200
d = 100  #100 50

constant1 <- 0.8
constant2 <- 0
constant3 <- 1
rrho1 <- 0.2
rrho2 <- 2 #rrho controls the smoothness of function. rrho = 0 means canstant.
epsilonscale<- sqrt(1)


set.seed(1)
diag11<- sample(c(1, 2),size = d, prob = rep(c(1/2, 1/2)), replace = TRUE)
diag21<- sample(c(-1, 0, 1),size = d, prob = rep(c(1/4, 1/2, 1/4)), replace = TRUE)
A01 <- bandSparse(d, k = c(0), diag = list(diag11), symm=FALSE)
#A02 <- bandSparse(d, k = -c(0:1), diag = list(rep(0,d), diag21), symm=TRUE)
A02 <- bandSparse(d, k = -c(0:2), diag = list(rep(0,d), rep(c(1,0),d/2),rep(c(0,0,0),d/3)), symm=TRUE)
Ai0 <- (A01==1) * (constant2+pnorm(1/n-1/2,0,rrho1) * constant1) +  (A01==2) * (constant2+ (1-pnorm(1/n-1/2,0,rrho1)) * constant1)
Ai20 <- constant3 * ((A02==1) * (pnorm(1/n-1/2,0,rrho1)*1.4-0.7) +  (A02==2) * (pnorm(1/n-1/2,0,rrho1)*1.4-0.7)) + diag(d)*1
#Ai20 <- crossprod(Ai20)


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
X[1,] <- constant1 * as.numeric(Ai0 %*% epsilon0)  + epsilon[1,] 
for (t in (2:n))
{
  Ai[[t-1]] <- (A01==1) * (constant2+pnorm(t/n-1/2,0,rrho1) * constant1) +  (A01==2) * (constant2+(1-pnorm(1/n-1/2,0,rrho1)) * constant1)
  Ai2 <- constant3 * ((A02==1) * (pnorm(t/n-1/2,0,rrho1)* 1.4-0.7) +  (A02==2) *  (pnorm(t/n-1/2,0,rrho1)* 1.4-0.7))  + diag(d)*1
  trueiCOV[[t-1]] <- Ai2
  trueCOV[[t-1]] <- solve(trueiCOV[[t-1]])
  epsilon[t,] <- mvrnorm(n=1,mu=rep(0,d),Sigma=trueCOV[[t-1]])*epsilonscale 
  X[t,] <- constant1 * as.numeric(Ai[[t-1]]%*%X[t-1,]) + epsilon[t,]
  maxeigen <- max(maxeigen,abs(var.roots(list(as.matrix(Ai[[t-1]]*constant1)))))
  mineigen <- min(mineigen, eigen(trueCOV[[t-1]])$values)
}

TrueCnet<-(A01!=0)*1
TrueAiL2norm <- as.matrix(sqrt(Reduce("+", lapply(Ai,function(x)x^2))))
TruePnet<-Reduce("|", lapply(trueiCOV,function(x)x!=0))*1
diag(TruePnet) <- rep(0,d)
#heatmap(A01,Rowv=NA,Colv=NA)



#A <- matrix(c(0.7,0.7,0,0,0.7,0.7,0, 0,0.7),3,3)
#eigen(A,symmetric = FALSE)$values
#svd(A)$d
