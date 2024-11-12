localpca <- function(data,         # n x d data matrix
                     hat.r=NULL,
                   tkernel=NULL, 
                   bandwidth=NULL,    # choice of bandwidth
                   eps=1e-5,
                   doparallel = FALSE,
                   max.r = 20,
                   ...           # arguments passed to tvCOV
)
{
  #  require(flare)
 
  n <- nrow(data)
  d <- ncol(data)
  
  # -------------------- Input Checks -------------------
  
  args <- list(...)
  # ----- Fill in Defaults -----
  if(is.null(tkernel))tkernel<- "Epa" #"Gaussian" #"Epa" #"Triweight"
  if(!'saveData' %in% args) args$saveData <- FALSE
  if(!'standardize' %in% args) args$standardize <- FALSE
  
# ----- Basic Input Checks -----
if(is.null(bandwidth)) #bandwidth <- 0.75*(log(d)/n)^0.2
  bandwidth <- (2.35/sqrt(12))*n^(-1/5)*d^(-1/10)#Su and Wang (2017)
  if(bandwidth <= 0) stop('The bandwidth parameter has to be strictly positive')
  
  # ----- Create Output Object -----
  
  localpca_object <- list('call' = NULL,
                       'COVoutput'= NULL)
  
  
  # ----- Copy the Call -----

  localpca_object$call <- list('data' = NULL,
                             'bandwidth' = bandwidth,
                              'args' = args
                             )
  
  if(args$saveData) localpca_object$call$data <- data

  


 # -------------------- Compute Weightings -------------------
  timevec <- seq(0,1,length=n+1)[2:(n+1)]
  ##PCA
  eigen_Cov<- list()
  for(i in 1:n){
    # Compute kernel weights for each estimation point
    l_Weights <- .kernel(timevec-timevec[i], bandwidth, tkernel = tkernel)
    COVoutput <- crossprod(data*l_Weights,data)/sum(l_Weights)
    eigen_Cov[[i]] <- eigen(COVoutput,symmetric=TRUE)
  } # end for :i (estpoints)
   
## number of factors
  
hat.err <- matrix(0,n,d)
V2 <- c()
for(num.fac in 1:max.r){
  for(i in 1:n){
    hat.lam <- eigen_Cov[[i]]$vectors[,1:num.fac,drop=FALSE] * sqrt(d) ##sign not determined
    hat.err[i,] <- (diag(d)-crossprod(t(hat.lam))/ d)%*% t(data[i,,drop=FALSE])
  }
  V2[num.fac] <- mean(hat.err^2)
}
localpca_object$IC2 <- log(V2) + (1:max.r) * (d+n * bandwidth)/(d*n * bandwidth) *log(max(d,n * bandwidth))

localpca_object$IC2 <- c(log(mean(data^2)),localpca_object$IC2)
if(is.null(hat.r))hat.r <- which.min(localpca_object$IC2)-1

localpca_object$hat.r <- hat.r
localpca_object$COVoutput <- eigen_Cov
localpca_object$V2 <- V2

localpca_object$hat.lam <- array(0,c(d,hat.r,n))
localpca_object$hat.fac <- array(0,c(n,hat.r))
localpca_object$hat.err <- data


### calculate fac and loading, and determine the sign
if(hat.r!=0){
  hat.lam.current <- eigen_Cov[[1]]$vectors[,1:hat.r,drop=FALSE] * sqrt(d)
  localpca_object$hat.lam[,,1] <- hat.lam.current
  localpca_object$hat.fac[1,] <-  data[1,,drop=FALSE] %*% hat.lam.current/d
  localpca_object$hat.err[1,] <- (diag(d)-crossprod(t(hat.lam.current))/d)%*% t(data[1,,drop=FALSE])
  for(i in 2:n){
    hat.lam.last <- hat.lam.current
    hat.lam.current <- eigen_Cov[[i]]$vectors[,1:hat.r,drop=FALSE] * sqrt(d)
    for (j in 1:hat.r){
      if (sum(abs(hat.lam.last[,j] - hat.lam.current[,j]))>sum(abs(hat.lam.last[,j] + hat.lam.current[,j]))){
        hat.lam.current[,j] <- (-1)* hat.lam.current[,j]
      }
    }  
    localpca_object$hat.lam[,,i] <- hat.lam.current
    localpca_object$hat.fac[i,] <-  data[i,,drop=FALSE] %*% hat.lam.current/d
    localpca_object$hat.err[i,] <- (diag(d)-crossprod(t(hat.lam.current))/ d)%*% t(data[i,,drop=FALSE])
  }
}

   # -------------------- Output -------------------
   return(localpca_object)
  

}
