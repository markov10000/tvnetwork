tvCOV <- function(data,         # n x p data matrix
                   tkernel=NULL, 
                   timepoints=NULL,
                   estpoints=NULL,    # vector of estimation points in [0,1]
                   bandwidth=NULL,    # choice of bandwidth
                   eps=1e-5,
                   ...           # arguments passed to tvCOV
)
{
  require(flare)
  n <- nrow(data)
  p <- ncol(data)
  
  # -------------------- Input Checks -------------------
  
  args <- list(...)
  # ----- Fill in Defaults -----
  if(is.null(tkernel))tkernel<- "Epa" #"Gaussian" #"Epa" #"Triweight"
  if(!'saveData' %in% args) args$saveData <- FALSE
  if(!'lambda' %in% args) args$lambda <- NULL
  if(!'lambdaGam' %in% args) args$lambdaGam <- 0.25
  if(!'nlambda' %in% args) args$nlambda<- 20
  if(!'lambda.min.ratio' %in% args)args$lambda.min.ratio<- 0.2
  if(!'rho' %in% args) args$rho<- NULL
  if(!'method' %in% args) args$method <- "clime"
  if(!'sym' %in% args) args$sym <- "and"
  if(!'shrink' %in% args) args$shrink <- NULL
  if(!'prec' %in% args) args$prec <- 1e-4
  if(!'max.ite' %in% args) args$max.ite <- 1e4
  if(!'standardize' %in% args) args$standardize <- FALSE
  if(!'perturb' %in% args) args$perturb <- FALSE
  if(!'verbose' %in% args) args$verbose <- FALSE
  if(!'warnings' %in% args) args$warnings <- FALSE
  if(!'verbatim' %in% args) args$verbatim <- FALSE
  if(!args$verbatim) args$warnings <- FALSE
  
  # Switch all warnings off
  if(!args$warnings) {
    oldw <- getOption("warn")
    options(warn = -1)
  }
  
# ----- Basic Input Checks -----
  if(is.null(timepoints)) timepoints <- seq(0, 1, length = n)
  if(is.null(estpoints)) estpoints <- seq(0, 1, length = n) 
  if(any(estpoints < 0 | estpoints > 1)) stop('Estimation points have to be specified on the unit interval [0,1].')
  if(is.null(bandwidth)) bandwidth <- 0.75*(log(d)/n)^0.2
  if(bandwidth <= 0) stop('The bandwidth parameter has to be strictly positive')
  
  # ----- Create Output Object -----
  
  tvCOV_object <- list('call' = NULL,
                       'COVoutput'= NULL,
                       'COVoutput'= NULL)
  
  
  # ----- Copy the Call -----

  tvCOV_object$call <- list('data' = NULL,
                             'timepoints' = timepoints,
                             'estpoints' = estpoints,
                             'estpointsNorm' = NULL,
                             'bandwidth' = bandwidth,
                              'args' = args
                             )
  
  if(args$saveData) tvCOV_object$call$data <- data

  


 # -------------------- Compute Weightings -------------------

  # Define time vector: if not provided, assume equally spaced time points
  tvCOV_object$call$timepoints <- timepoints
  # normalize to [0,1]
  timevec <- (timepoints - min(timepoints)) / max(timepoints)
  tvCOV_object$call$timepoints_cut <- timevec

  # Normalize time estimation points to interval [0,1]
  estpoints_norm <- estpoints
  tvCOV_object$call$estpointsNorm <- estpoints_norm
  no_estpoints <- length(estpoints_norm)

  lambda.max <- pi*sqrt(log(d)/(2*n*bandwidth)) * min(apply(data,2,sd))
  lambda.min <- args[['lambda.min.ratio']] * lambda.max
  lambda_seq <- exp(seq(log(lambda.max), log(lambda.min), 
                   length = args[['nlambda']]))
  
# Compute kernel weights for each estimation point
  COVoutput<-iCOVoutput<-list()
  
   lambda <- EBIC <- rep(0,no_estpoints)
  for(i in 1:no_estpoints) {
    K0 <- .kernel(timevec-estpoints_norm[i], bandwidth, tkernel = tkernel)
    K1 <- K0 * (timevec-estpoints_norm[i])/bandwidth
    K2 <- K1 * (timevec-estpoints_norm[i])/bandwidth
    s1 <- sum(K1)
    s2 <- sum(K2)
    l_Weights <- K0 * s2 - K1 * s1
    n_eff <- n * (min(estpoints_norm[i] + bandwidth, 1) - max(estpoints_norm[i] - bandwidth, 0))
    n_eff <- sum(l_Weights)/max(l_Weights)
    perturb <- sqrt(1/(2*n*bandwidth))
    COVoutput[[i]] <- crossprod(data*l_Weights,data)/sum(l_Weights)
    eigen_iCov <- eigen(COVoutput[[i]],symmetric=TRUE)
    #eigen_iCov$vectors %*% (sapply(eigen_iCov$values,function(x)sqrt(max(x,0))) * t(eigen_iCov$vectors))
    COVoutput[[i]] <- eigen_iCov$vectors %*% (sapply(eigen_iCov$values,function(x)sqrt(max(x,0))) * t(eigen_iCov$vectors))#

    
    EBICicov<-c()
    #Package flare
    out_e <- flare::sugm(COVoutput[[i]] + diag(d) * perturb, lambda = lambda_seq, rho = args[['rho']], method = args[['method']], sym = args[['sym']], shrink=args[['shrink']], prec = args[['prec']], max.ite = args[['max.ite']], standardize = args[['standardize']], perturb = args[['perturb']], verbose = args[['verbose']])
    for(j in 1:length(out_e$icov)){
      #EBICicov[j] <- n_eff * (-log(det(out_e$icov[[j]]))+sum(diag(out_e$icov[[j]]%*%COVoutput[[i]]))) + (log(n_eff) + 4*args$lambdaGam*log(d)) * (sum(abs(out_e$icov[[j]])>out_e$lambda[j])-sum(abs(diag(out_e$icov[[j]]))>out_e$lambda[j]))/2
      EBICicov[j] <- n_eff * (-log(det(out_e$icov[[j]]))+sum(diag(out_e$icov[[j]]%*%COVoutput[[i]]))) + (log(n_eff) + 4*args$lambdaGam*log(d)) * (sum(abs(out_e$icov[[j]])>0)-d)/2
    }
    EBIC.min.idx <- which.min(EBICicov)
    EBIC[i] <- EBICicov[EBIC.min.idx]
    lambda[i] <- lambda_seq[EBIC.min.idx]
    iCOVoutput[[i]] <- out_e$icov[[EBIC.min.idx]]
    #Package fastclime seems not robust to covariance matrix with negative eigenvalues  
     #out1 <- fastclime::fastclime(COVoutput[[i]], lambda.min = lambda.min, nlambda = args$nlambda)
     #for(j in 1:length(lambda_seq)){
    # out2 <- fastclime.selector(out1$lambdamtx, out1$icovlist,lambda_seq[j])
    #   EBICicov[j] <- n_eff * (-log(det(out2$icov))+sum(diag(out2$icov%*%COVoutput[[i]]))) + (log(n_eff) + 4*args$lambdaGam*log(d)) * (sum(out2$icov!=0)-d)/2
     #}
     #EBIC.min.idx <- which.min(EBICicov)
     #EBIC[i] <- EBICicov[EBIC.min.idx]
     #lambda[i] <- lambda_seq[EBIC.min.idx]
     #iCOVoutput[[i]] <- fastclime.selector(out1$lambdamtx, out1$icovlist,lambda_seq[EBIC.min.idx])$icov
  } # end for:i (estpoints)

   Pnet<- Matrix(rep(FALSE,d*d),nrow=d,ncol=d, sparse = TRUE)
   for (i in 1:no_estpoints){
     Pnet <- Pnet | (abs(iCOVoutput[[i]])>lambda[i])
   }
  diag(Pnet)<- rep(FALSE,d)
   # -------------------- Output -------------------
  tvCOV_object$COVoutput <- COVoutput
  tvCOV_object$iCOVoutput <- iCOVoutput
  tvCOV_object$Pnet <- Pnet
  tvCOV_object$perturb <- perturb
  tvCOV_object$EBIC <- EBIC
  tvCOV_object$lambda_seq <- lambda_seq
  tvCOV_object$lambda <- lambda
  return(tvCOV_object)
  

}
