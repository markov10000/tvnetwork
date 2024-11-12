var.roots<-function (coef.list, modulus = TRUE) 
{

  K <- dim(coef.list[[1]])
  p <- length(coef.list)
  A <- unlist(coef.list)
  companion <- matrix(0, nrow = K * p, ncol = K * p)
  companion[1:K, 1:(K * p)] <- A
  if(p>1) diag(companion[(K+1):(K * p), 1:(K * p)]) <- 1
  roots <- eigen(companion)$values
  if (modulus) 
    roots <- Mod(roots)
  return(roots)
}