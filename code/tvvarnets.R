tvvarnets <- function(data,         # n x p data matrix
                   lags=1,         # lag of VAR model
                   timepoints=NULL,
                   estpoints=NULL,    # vector of estimation points in [0,1]
                   bandwidth=NULL,    # choice of bandwidth
                   lambda=NULL,    # choice of lambda
                   lambda2.fac=NULL,
                   adj_GIC=1,
                   eps=1e-5,
                   max.iter=1000,
                   do.parallel = FALSE,
                   ...           # arguments passed to tvvarnets
)
  

  {
  n_lags <- length(lags)
  max_lag <- max(lags)
  n <- nrow(data)
  p <- ncol(data)
  n_var <-
    n - max_lag # this is how many rows there are after transforming the data
  
  # -------------------- Input Checks -------------------
  
  args <- list(...)
  # ----- Fill in Defaults -----
  if (is.null(args$lambdaSeq))
    args$lambdaSeq <- NULL
  if (is.null(args$lambdaSel))
    args$lambdaSel <- 'AIC' # AIC, BIC or EBIC
  if (is.null(args$lambdaGam))
    args$lambdaGam <- .25
  if (is.null(args$glmGam))
    args$glmGam <- 1
  if (is.null(args$tkernel))
    args$tkernel <- "Epa" #"Gaussian" #"Epa" #"Triweight"
  if (is.null(args$standardize))
    args$standardize <- FALSE
  if (is.null(args$intercept))
    args$intercept <- FALSE
  if (is.null(args$warnings))
    args$warnings <- FALSE
  if (is.null(args$saveModels))
    args$saveModels <- FALSE
  if (is.null(args$saveData))
    args$saveData <- FALSE
  if (is.null(args$pbar))
    args$pbar <- pbar <-  FALSE
  if (is.null(args$verbatim))
    args$verbatim <- FALSE
  if (args$verbatim)
    args$pbar <- FALSE
  if (args$verbatim)
    args$warnings <- FALSE
  
    
    # Switch all warnings off
    if (!args$warnings) {
      oldw <- getOption("warn")
      options(warn = -1)
    }
    
    # ----- Basic Input Checks -----
    if (is.null(timepoints))
      timepoints <- seq(0, 1, length = n)
    if (is.null(estpoints))
      estpoints <- seq(0, 1, length = n_var)
    n_est <- length(estpoints)
    estidx <- floor(n_var * (estpoints))
    estpoints <- estidx / n_var
    if (any(estpoints < 0 |
            estpoints > 1))
      stop('Estimation points have to be specified on the unit interval [0,1].')
    if (is.null(bandwidth))
      bandwidth <- 0.75 * (log(d) / n) ^ 0.2
    if (bandwidth <= 0)
      stop('The bandwidth parameter has to be strictly positive')
    if (is.null(lambda))
      lambda <- list()
    if (length(lambda$lambda1) == 1)
      lambda1 <- array(lambda$lambda1, c(p, n_est))
    if (length(lambda$lambda1) == p)
      lambda1 <- lambda$lambda1
    if (is.null(lambda2.fac))
      lambda2.fac <- 1 / (1 + seq(0, 20, 1))
    
    # ----- Create Output Object -----
    
    tvvarnets_object <- list(
      'call' = NULL,
      'wadj' = NULL,
      'signs' = NULL,
      'edgecolor' = NULL,
      'intercepts' = NULL,
      'tvmodels' = NULL
    )
    
    
    # ----- Copy the Call -----
    
    call <- list(
      'data' = NULL,
      'timepoints' = timepoints,
      'estpoints' = estpoints,
      'estpointsNorm' = NULL,
      'bandwidth' = bandwidth,
      'lags' = lags,
      'lambdaSeq' = args$lambdaSeq,
      'lambdaSel' = args$lambdaSel,
      'lambdaFolds' = args$lambdaFolds,
      'lambdaGam' = args$lambdaGam,
      'warnings' = args$warnings,
      'saveModels' = args$saveModels,
      'saveData' = args$saveData
    )
    
    
    if (args$saveData)
      call$data <- data
    
    
    
    
    # -------------------- Compute Weightings -------------------
    
    # Define time vector: if not provided, assume equally spaced time points
    # normalize to [0,1]
    timepoints <-
      timepoints[-(1:max(lags))] # delete first x rows that have to be exluded by definition of VAR model
    timepoints <- timepoints - min(timepoints)
    timevec <- timepoints / max(timepoints)
    
    # Normalize time estimation points to interval [0,1]
    estpoints_norm <- estpoints
    call$estpointsNorm <- estpoints_norm
    n_est <- length(estpoints_norm)
    
    # Compute kernel weights for each estimation point
    l_weights <- list()
    for (i in 1:n_est) {
      l_weights[[i]] <-
        .kernel(timevec - estpoints_norm[i], bandwidth, tkernel = args$tkernel) /
        bandwidth
      # Normalize to [x,1]
      l_weights[[i]] <- l_weights[[i]] / sum(l_weights[[i]])
    }
    n_wadj <- sapply(l_weights, sum)
    
    # -------------------- Create VAR data structure -------------------
    
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
    data_response <- data_response[data_lagged$included,]
    l_data_lags <-
      lapply(l_data_lags, function(x)
        x[data_lagged$included,])
    rm(data_lagged)
    
    # Create design matrix
    
    data_localXX <- list()
    #lasso_models <- list()
    #adlasso_models <- list()
    sdy <- apply(data_response, 2, sd)
    data_localX <- do.call(cbind, l_data_lags)
    lasso_wadj <- lasso_slope <- array(0, c(p, p, length(lags), n_est))
    adlasso_wadj <- adlasso_slope <- array(0, c(p, p, length(lags), n_est))
    lasso_wadj_s <- array(0, c(p, p, length(lags), n_est))
    for (i in 1:n_est) {
      data_localXX[[i]] <-
        cbind(data_localX,
              data_localX * (timevec - estpoints_norm[i]) / bandwidth)
    }
    p_design <- ncol(data_localX) * 2 ##local linear
    
    # Select lambda 1
    # y is standardized by sd(y)
    lambda1_Sel <- NULL
    if (is.null(lambda$lambda1)) {
      # Create Progress Bar
      if (args$pbar == TRUE)
        pb <-
          txtProgressBar(
            min = 0,
            max = p - 1,
            initial = 0,
            char = "-",
            style = 3
          )
      if (args$lambdaSel == 'EBIC') {
        parres <- foreach(
          j = 0:(p - 1),
          .combine = 'comb',
          .multicombine = TRUE,
          .init = list(list(), list(), list()),
          .packages = c('glmnet')
        ) %dopar% {
          lambda1_j <- array(0, c(n_est))
          lasso_wadj_j <-
            lasso_slope_j <- array(0, c(p, length(lags), n_est))
          for (i in 1:n_est) {
            glmnet_fit0 <-
              glmnet(
                x = data_localXX[[i]],
                y = data_response[, j + 1],
                weights = l_weights[[i]],
                standardize = args$standardize,
                lambda = args$lambdaSeq,
                gamma = args$glmGam
              )
            EBIC_lambda <-
              deviance(glmnet_fit0) + glmnet_fit0$df * (log(n_wadj[i]) + 2 * args$lambdaGam * log(p_design))
            lambda1_j[i] <- glmnet_fit0$lambda[which.min(EBIC_lambda)]
            lambad_min_model <- coef(glmnet_fit0, s = lambda1_j[i])
            for (lagi in 1:n_lags) {
              lasso_wadj_j[, lagi, i] <-
                lambad_min_model[((lagi - 1) * p + 1 + 1):(lagi * p + 1)]
              lasso_slope_j[, lagi, i] <-
                lambad_min_model[((lagi - 1) * p + 1 + 1 + ncol(data_localX)):(lagi * p +
                                                                                 1 + ncol(data_localX))]
            }
            # Update Progress Bar
            if (args$pbar == TRUE)
              setTxtProgressBar(pb, j)
            return(
              list(
                lasso_wadj = lasso_wadj_j,
                lasso_slope = lasso_slope_j,
                lambda1_Sel = lambda1_j
              )
            )
          }
        }#end parallel
        lasso_wadj <- parres[[1]]
        lasso_slope <- parres[[2]]
        lambda1_Sel <- parres[[3]]
        lasso_wadj <- aperm(lasso_wadj, c(4, 1, 2, 3))
        lasso_slope <- aperm(lasso_slope, c(4, 1, 2, 3))
        lambda1_Sel <- aperm(lambda1_Sel, c(2, 1))
      }
      if (args$lambdaSel == 'AIC') {
        parres <- foreach(
          j = 0:(p - 1),
          .combine = 'comb',
          .multicombine = TRUE,
          .init = list(list(), list(), list()),
          .packages = c('glmnet')
        ) %dopar% {
          lambda1_j <- array(0, c(n_est))
          lasso_wadj_j <-
            lasso_slope_j <- array(0, c(p, length(lags), n_est))
          fit.list <- list()
          for (i in 1:n_est) {
            glmnet_fit0 <-
              glmnet(
                x = data_localXX[[i]],
                y = data_response[, j + 1],
                weights = l_weights[[i]],
                standardize = args$standardize,
                intercept = args$intercept,
                lambda = args$lambdaSeq,
                gamma = args$glmGam
              )
            #n_eff <- n_var * (min(estpoints_norm[i] + bandwidth, 1) - max(estpoints_norm[i] - bandwidth, 0))
            n_eff <- sum(l_weights[[i]]) / max(l_weights[[i]])
            AIC_lambda <-
              deviance(glmnet_fit0) + 2 / n_eff * glmnet_fit0$df
            lambda1_j[i] <- glmnet_fit0$lambda[which.min(AIC_lambda)]
            lambad_min_model <- coef(glmnet_fit0, s = lambda1_j[i])
            for (lagi in 1:n_lags) {
              lasso_wadj_j[, lagi, i] <-
                lambad_min_model[((lagi - 1) * p + 1 + 1):(lagi * p + 1)]
              lasso_slope_j[, lagi, i] <-
                lambad_min_model[((lagi - 1) * p + 1 + 1 + ncol(data_localX)):(lagi * p +
                                                                                 1 + ncol(data_localX))]
            }
          }
          # Update Progress Bar
          if (args$pbar == TRUE)
            setTxtProgressBar(pb, j)
          return(
            list(
              lasso_wadj = lasso_wadj_j,
              lasso_slope = lasso_slope_j,
              lambda1_Sel = lambda1_j
            )
          )
        }
        lasso_wadj <- simplify2array(parres[[1]])
        lasso_slope <- simplify2array(parres[[2]])
        lambda1_Sel <- simplify2array(parres[[3]])
        lasso_wadj <- aperm(lasso_wadj, c(4, 1, 2, 3))
        lasso_slope <- aperm(lasso_slope, c(4, 1, 2, 3))
        lambda1_Sel <- aperm(lambda1_Sel, c(2, 1))
      }
      if (args$lambdaSel == 'BIC') {
        parres <- foreach(
          j = 0:(p - 1),
          .combine = 'comb',
          .multicombine = TRUE,
          .init = list(list(), list(), list()),
          .packages = c('glmnet')
        ) %dopar% {
          lambda1_j <- array(0, c(n_est))
          lasso_wadj_j <-
            lasso_slope_j <- array(0, c(p, length(lags), n_est))
          fit.list <- list()
          for (i in 1:n_est) {
            glmnet_fit0 <-
              glmnet(
                x = data_localXX[[i]],
                y = data_response[, j + 1],
                weights = l_weights[[i]],
                standardize = args$standardize,
                intercept = args$intercept,
                lambda = args$lambdaSeq,
                gamma = args$glmGam
              )
            #n_eff <- n_var * (min(estpoints_norm[i] + bandwidth, 1) - max(estpoints_norm[i] - bandwidth, 0))
            n_eff <- sum(l_weights[[i]]) / max(l_weights[[i]])
            BIC_lambda <-
              deviance(glmnet_fit0) + log(n_eff) / n_eff * glmnet_fit0$df
            lambda1_j[i] <- glmnet_fit0$lambda[which.min(BIC_lambda)]
            lambad_min_model <- coef(glmnet_fit0, s = lambda1_j[i])
            for (lagi in 1:n_lags) {
              lasso_wadj_j[, lagi, i] <-
                lambad_min_model[((lagi - 1) * p + 1 + 1):(lagi * p + 1)]
              lasso_slope_j[, lagi, i] <-
                lambad_min_model[((lagi - 1) * p + 1 + 1 + ncol(data_localX)):(lagi * p +
                                                                                 1 + ncol(data_localX))]
            }
          }
          # Update Progress Bar
          if (args$pbar == TRUE)
            setTxtProgressBar(pb, j)
          return(
            list(
              lasso_wadj = lasso_wadj_j,
              lasso_slope = lasso_slope_j,
              lambda1_Sel = lambda1_j
            )
          )
        }
        lasso_wadj <- simplify2array(parres[[1]])
        lasso_slope <- simplify2array(parres[[2]])
        lambda1_Sel <- simplify2array(parres[[3]])
        lasso_wadj <- aperm(lasso_wadj, c(4, 1, 2, 3))
        lasso_slope <- aperm(lasso_slope, c(4, 1, 2, 3))
        lambda1_Sel <- aperm(lambda1_Sel, c(2, 1))
      }
      
    }
    L2norm_alpha <-
      apply(lasso_wadj, c(1, 2, 3), function(x)
        sqrt(sum(x ^ 2)))
    sd_alpha <-
      apply(lasso_wadj, c(1, 2, 3), function(x)
        sd(x) * sqrt(n_est - 1))
    L2norm_alpha <- t(apply(L2norm_alpha, 1, rbind))
    sd_alpha <-  t(apply(sd_alpha, 1, rbind))
    
    
    
    # -------------------- Global Estimation Combining Estpoints -------------------
    glasso_wadj <- glasso_slope <- array(0, c(p, p, length(lags), n_est))
    glasso_models <- list()
    loss <- c()
    iter <- c()
    
    # Create Progress Bar
    if (args$pbar == TRUE)
      pb2 <-
      txtProgressBar(
        min = 0,
        max = p - 1,
        initial = 0,
        char = "-",
        style = 3
      )
    
    r <- rep(NA_real_, n_est * n_est)
    group <- rep(1:p_design, n_est)
    XX <-
      do.call(cbind, mapply('*', data_localXX, lapply(l_weights, sqrt), SIMPLIFY = FALSE))
    
    lambda2_Sel <- rep(0, p)
    parres <- foreach(
      j = 0:(p - 1),
      .combine = 'comb',
      .multicombine = TRUE,
      .init = list(list(), list(), list(), list(), list(), list())
    ) %dopar% {
      #source('code/SCAD.derivative.R')
      #dyn.load("code/rawfit_wglasso.so")
      lambda1_j <- array(0, c(n_est))
      glasso_wadj <- glasso_slope <- array(0, c(p, length(lags), n_est))
      if (max(L2norm_alpha[j + 1, ]) == 0) {
        for (i in 1:n_est) {
          for (lagi in 1:n_lags) {
            glasso_wadj[, lagi, i] <- rep(0, p)
            glasso_slope[, lagi, i] <- rep(0, p)
          }
        }
        glasso_models <- list()
        loss <- c(loss, NA)
        GIC.min <- NA
        GIC <- NA
        iter <- c(iter, 0)
      } else{
        YY <- c(sapply(l_weights, function(x)
          data_response[, j + 1] * sqrt(x)))
        #glasso_models[[j*n_est+i]] <- gglasso(x=do.call(blkdiag,XX),y=do.call(c, YY),group=rep(1:p_design,n_est))
        #ncvreg(do.call(blkdiag,XX), do.call(c, YY), penalty="lasso", lambda=lambda[1])
        ginit <- list()
        for (l in 1:n_est) {
          ginit[[l]] <- c(lasso_wadj[j + 1, , , l], lasso_slope[j + 1, , , l])
        }
        ginit <- do.call(c, ginit)
        maxeigen <- rep(0, p_design) #p_design=p*n_lags*2
        for (k in 1:p_design) {
          XXk <- XX[, seq(k, (n_est - 1) * p_design + k, by = p_design)]
          maxeigen[k] <-
            max(apply(XXk, 2, function(x)
              sum(x ^ 2))) / n_est / n_var
        }

        gmodel <- list()
        
        penalty.factor <-
          SCAD.derivative(c(L2norm_alpha[j + 1, ], sd_alpha[j + 1, ]), max(L2norm_alpha[j +
                                                                                          1, ]))
        gmodel[[1]] <-
          .Call(
            "rawfit_glasso",
            XX,
            YY,
            ginit,
            n_est,
            r,
            maxeigen,
            1,
            eps,
            as.integer(max.iter),
            penalty.factor
          )
        gnspa <- c(sum(matrix(gmodel[[1]][[1]],nrow=p_design)[1:p,]!=0))
        #gnspa <-
        #  c(sum(apply(matrix(gmodel[[1]][[1]], nrow = p_design)[1:(p * n_lags), ] !=
        #                       0, 1, function(x)any(x))))
        # calculate residuals
        gres <- data_response[estidx, j + 1]
        for (i in 1:(n_var)) {
          for (lagi in 1:n_lags) {
            gres[i] <-
              gres[i] - gmodel[[1]][[1]][((i - 1) * p_design + (lagi - 1) * p + 1):((i -
                                                                                       1) * p_design + lagi * p)] %*% l_data_lags[[lagi]][i, ]
          }
        }
        SSE <- c(mean((gres - mean(gres)) ^ 2))
        #SSE <- gmodel[[1]][[2]]/n_var

        for (ind.lambda2 in 2:length(lambda2.fac)) {
          penalty.factor <-
            SCAD.derivative(c(L2norm_alpha[j + 1, ], sd_alpha[j + 1, ]),
                            max(L2norm_alpha[j + 1, ]) * lambda2.fac[ind.lambda2])
          gmodel[[ind.lambda2]] <-
            .Call(
              "rawfit_glasso",
              XX,
              YY,
              gmodel[[ind.lambda2 - 1]][[1]],
              n_est,
              r,
              maxeigen,
              1,
              eps,
              as.integer(max.iter),
              penalty.factor
            )
          #gnspa <- c(gnspa, sum(apply(matrix(gmodel[[ind.lambda2]][[1]], nrow = p_design)[1:(p * n_lags), ] !=0,
          #  1, function(x)any(x))))
          gnspa <- c(gnspa,sum(matrix(gmodel[[ind.lambda2]][[1]],nrow=p_design)[1:(p*n_lags),]!=0))
          # calculate residuals
          gres <- data_response[, j + 1]
          for (i in 1:n_est) {
            for (lagi in 1:n_lags) {
              gres[i] <-
                gres[i] - gmodel[[ind.lambda2]][[1]][((i - 1) * p_design + (lagi - 1) *
                                                       p + 1):((i - 1) * p_design + lagi * p)] %*% l_data_lags[[lagi]][i, ]
            }
          }
          SSE <- c(SSE, mean((gres - mean(gres)) ^ 2))
          #SSE <- c(SSE, gmodel[[ind.lambda2]][[2]]/n_var)
        }
        gnspa <- gnspa/n_var
        GIC <-
          n_var * log(SSE) + adj_GIC * log(log(n_var)) * log(36 / 35 * p * n_lags / bandwidth) * 36 /
          35 / bandwidth * gnspa
        #GIC <- n_var * log(SSE) +  log(n_var) * 1.02857/bandwidth/n_var * gnspa
        GIC.min <- which.min(GIC)
        glasso_models <- gmodel[[GIC.min]]
        loss <- glasso_models[[2]]
        iter <- glasso_models[[3]]
        lambda2_Sel <- max(L2norm_alpha[j + 1, ]) * lambda2.fac[GIC.min]
        rm(gmodel)
        for (i in 1:n_est) {
          for (lagi in 1:n_lags) {
            glasso_wadj[, lagi, i] <-
              glasso_models[[1]][((i - 1) * p_design + (lagi - 1) * p + 1):((i - 1) *
                                                                              p_design + lagi * p)]
            glasso_slope[, lagi, i] <-
              glasso_models[[1]][((i - 1) * p_design + (lagi - 1) * p + 1 + ncol(data_localX)):((i -
                                                                                                   1) * p_design + lagi * p + ncol(data_localX))]
          }
        }
        
        # Update Progress Bar
        if (args$pbar == TRUE)
          setTxtProgressBar(pb2, j + 1)
      }
      return(
        list(
          glasso_wadj = glasso_wadj,
          glasso_slope = glasso_slope,
          lambda2_Sel = lambda2_Sel,
          loss = loss,
          GIC = GIC[GIC.min],
          iter = iter
        )
      )#glasso_models=glasso_models
    }#end parallel
    glasso_wadj <- simplify2array(parres[[1]])
    glasso_slope <- simplify2array(parres[[2]])
    glasso_wadj <- aperm(glasso_wadj, c(4, 1, 2, 3))
    glasso_slope <- aperm(glasso_slope, c(4, 1, 2, 3))
    lambda2_Sel <- simplify2array(parres[[3]])
    loss <- simplify2array(parres[[4]])
    GIC.min <- simplify2array(parres[[5]])
    iter <- simplify2array(parres[[6]])
    
    names(loss) <- paste("V", 1:p, '', sep = "")
    
    # -------------------- Local Estimation of Precision Matrix over Estpoints -------------------
    lasso_residuals <- adlasso_residuals <- glasso_residuals <- NA
    
    #check!
    #residuals <- get.residuals(data=data, estpoints=estpoints, lags=lags, coef=list(lasso_wadj,glasso_wadj))
    #lasso_residuals <- residuals[[1]]
    #glasso_residuals <- residuals[[2]]
    
    
    
    # -------------------- Process Output -------------------
    # Restructure output to make time-varying model easier accessible
    
    
    #if (warn & sum(iter)==max.iter) warning("Maximum number of iterations reached")
    
    
    
    # -------------------- Output -------------------
    tvmvar_object <- list(
      lasso_wadj = lasso_wadj,
      lasso_slope = lasso_slope,
      lasso_residuals = lasso_residuals,
      glasso_wadj = glasso_wadj,
      glasso_slope = glasso_slope,
      glasso_residuals = glasso_residuals,
      L2norm_alpha = L2norm_alpha,
      sd_alpha = sd_alpha,
      loss = loss,
      gic = GIC.min,
      iter = iter,
      sdy = sdy,
      lambda1 = lambda1_Sel,
      lambda2 = lambda2_Sel,
      call = call
    )
    if (args$saveModels) {
      tvmvar_object$glasso_models <- parres[[7]]
      names(v) <- paste("V", 1:p, '', sep = "")
    }
    class(tvmvar_object) <- c('mgm', 'tvmvar')
    
    return(tvmvar_object)
    
    
  }
  
  Lambdas <- function(...) {
    cv <- cv.glmnet(...)
    return(data.table(cvm = cv$cvm, lambda = cv$lambda))
  }
  
  OptimLambda <- function(k, ...) {
    # Returns optimal lambda for glmnet.
    #
    # Args:
    #   k: # times to loop through cv.glmnet.
    #   ...: Other args passed to cv.glmnet.
    #
    # Returns:
    #   Lambda associated with minimum average CV error over runs.
    #
    # Example:
    #   OptimLambda(k=100, y=y, x=x, alpha=alpha, nfolds=k)
    #
    require(parallel)
    require(data.table)
    require(plyr)
    MSEs <-
      data.table(rbind.fill(mclapply(seq(k), function(dummy)
        Lambdas(...))))
    return(MSEs[, list(mean.cvm = mean(cvm)), lambda][order(mean.cvm)][1]$lambda)
  }
  
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i)
             c(x[[i]], lapply(list(...), function(y)
               y[[i]])))
  }