get.residuals <- function(data,         # n x p data matrix
                      estpoints=NULL,    # vector of estimation points in [0,1]
                      lags=1,         # lag of VAR model
                      coef){
  if(!is.list(coef))coef <- list(coef=coef)
  modeln <- length(coef)
  n_lags <- length(lags)
  n <- nrow(data)
  p <- ncol(data)
  n_var <- n - max(lags) # this is how many rows there are after transforming the data
  
  if(is.null(estpoints)) estpoints <- seq(0, 1, length = n_var)
  n_est <- length(estpoints)
  estidx <- floor(n_var * (estpoints))
  estpoints <- estidx/n_var
  
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
  for(modeli in 1:modeln) res[[modeli]]<-data_response[estidx,]
  print(length(estidx))
  print(dim(data_response[estidx,]))
  for(i in 1:n_est){
    for(lagi in 1:n_lags){
      for(modeli in 1:modeln){
        res[[modeli]][i,] <- res[[modeli]][i,] - coef[[modeli]][,,lagi,i]%*%l_data_lags[[lagi]][i,]
      }
    }
  }
  return(res)
}