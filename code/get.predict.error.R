get.predict.error <- function(data,         # n x p data matrix
                      lags=1,         # lag of VAR model
                      pre_lag =1,
                      discard=0,
                      coef){
  if(!is.list(coef))coef <- list(coef=coef)
  modeln <- length(coef)
  n_lags <- length(lags)
  n <- nrow(data)
  p <- ncol(data)
  n_var <- n - max(lags) # this is how many rows there are after transforming the data
  
  # ----- Give Names to Variables & create data.frame -----
  
  colnames(data)[1:p] <- paste("V", 1:p, '.', sep = "")
  data <- as.data.frame(data)
  
  
  # ----- Split up predictor Sets by lags -----
  
  # Divide Data in several parts: one response set, and one set of predictors for each lag
  data_lagged <- lagData(data = data, 
                         lags = lags)
  
  data_response <- data_lagged$data_response
  l_data_lags <- data_lagged$l_data_lags
  
  # delete rows that cannot be predicted
  data_response <- data_response[data_lagged$included, ]
  l_data_lags <- lapply(l_data_lags, function(x) x[data_lagged$included, ])
  
  
  # calculate residuals
  res <- list()
  for(modeli in 1:modeln) res[[modeli]]<-data_response[(1+discard+pre_lag):(n_var-discard),]
  for(i in 1:(n_var-2*discard-pre_lag)){
    for(lagi in 1:n_lags){
      for(modeli in 1:modeln){
        res[[modeli]][i,] <- res[[modeli]][i,] - coef[[modeli]][,,lagi,i+discard]%*%l_data_lags[[lagi]][i+discard+pre_lag,]
      }
    }
  }
  return(res)
}