rewLASSO <- function(mvar_object, data) {
  # -------------------- Local Estimation Loop over Estpoints -------------------
  n <- nrow(data)
  p <- ncol(data)
  lags <- 1
  lambda2.fac <- 1 / (1 + seq(0, 20, 1))
  
  data_lagged <- lagData(data = data,
                         lags = lags)
  
  data_response <- data_lagged$data_response
  l_data_lags <- data_lagged$l_data_lags
  # delete rows that cannot be predicted
  data_response <- data_response[data_lagged$included,]
  l_data_lags <-
    lapply(l_data_lags, function(x)
      x[data_lagged$included,])
  
  SSE <- rep(0, length(lambda2.fac))
  GIC <- rep(0, length(lambda2.fac))
  gnspa <- rep(0, length(lambda2.fac))
  adlasso_wadj_pre <- array(0, c(p, p, 1, length(lambda2.fac)))
  adlasso_wadj <- array(0, c(p, p, 1))
  
  
  for (j in 0:(p - 1)) {
    if (sum(mvar_object$wadj[j + 1, , 1]!=0)>0){
    for (ind.lambda2 in 1:length(lambda2.fac)) {
      #glmnet_fit1 <- glmnet(x=data_localXX, y=data_response[, j+1], weights=l_weights[[i]], standardize = args$standardize,gamma=args$gamma,lambda = lambda1[j+1,i])
      penalty.factor <-
        SCAD.derivative(mvar_object$wadj[j + 1, , 1], lambda =  max(abs(mvar_object$wadj[j + 1, , 1])) * lambda2.fac[ind.lambda2]) #n^(7/30)
      glmnet_fit2 <-
        glmnet(
          x = l_data_lags[[1]],
          y = data_response[, j + 1],
          lambda = 1,
          penalty.factor = penalty.factor
        )
      for (lagi in 1:1) {
           adlasso_wadj_pre[j + 1, , lagi, ind.lambda2] <- coef(glmnet_fit2)[((lagi - 1) * p + 1 + 1):(lagi * p + 1)]
      }
      err <-
        data_response[, j + 1] - l_data_lags[[lagi]] %*% adlasso_wadj_pre[j + 1, , lagi,ind.lambda2]
      SSE[ind.lambda2] <- mean((err - mean(err)) ^ 2)
      gnspa[ind.lambda2] <- sum(adlasso_wadj_pre[j + 1, , lagi,ind.lambda2]!=0)
    }
    adj_GIC <- 1
    GIC <-
      (n - 1) * log(SSE) + adj_GIC * log(log(n-1)) * log(p) * gnspa
    #GIC <- n_var * log(SSE) +  log(n_var) * 1.02857/bandwidth/n_var * gnspa
    GIC.min <- which.min(GIC)
    adlasso_wadj[j + 1, , 1] <- adlasso_wadj_pre[j + 1, , 1, GIC.min]
    }
    else{adlasso_wadj[j + 1, , 1] <- 0}
    }
  
  return(list(adlasso_wadj=adlasso_wadj))
}