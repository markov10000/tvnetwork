

calcLL <- function(X,
                   y,
                   fit, #glmnet fit object
                   weights,
                   lambda,
                   LLtype = 'model')


  {

  # This function calculates three different LL:
  # 1) LLtype = 'model': The LL of a given model via fit
  # 2) LLtype = 'nullmodel': The LL of the Null (Intercept) Model
  # 3) LLtype = 'saturated': The LL of the saturated model

  n <- nrow(X)

  if(LLtype == 'model') {
    beta_vector <- matrix(coef(fit, s = lambda), ncol = 1)
    predicted_mean <- cbind(rep(1, n), X) %*% as.vector(beta_vector)
    LL_model <- dnorm(y, mean = predicted_mean, sd = 1, log = TRUE)
    mean_LL_model <- sum(LL_model * weights)
  }

  if(LLtype == 'nullmodel') {
    beta_vector <- matrix(coef(fit, s = 1)[1], ncol = 1) # only intercept here
    predicted_mean <- rep(1, n) * as.vector(beta_vector)
    LL_model <- dnorm(y, mean = predicted_mean, sd = 1, log = TRUE)
    mean_LL_model <- sum(LL_model * weights)  
  }


  if(LLtype == 'saturated') {
    predicted_mean <- y
    LL_model <- dnorm(y, mean = predicted_mean, sd = 1, log = TRUE)
    mean_LL_model <- sum(LL_model * weights)
  }

  return(mean_LL_model)

}
