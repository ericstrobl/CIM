
gcm.test <- function(X, Y, Z = NULL, alpha = 0.05, regr.method = "xgboost", regr.pars = list(),
                     plot.residuals = FALSE, nsim=499L, resid.XonZ=NULL, resid.YonZ=NULL) {

  X = matrix2(X); Y=matrix2(Y);
  r = nrow(X)
  X = rank(X,ties.method="random")/r - 0.5; X = matrix2(X)
  Y = rank(Y,ties.method="random")/r - 0.5; Y = matrix2(Y)
  if (!is.null(Z)){
    Z=matrix2(Z)
    for (i in 1:ncol(Z)){
      Z[,i] = rank(Z[,i],ties.method="random")/r - 0.5
    }
  }
  # X = normalize01(X);
  # Y = normalize01(Y);
  # Z = normalize01(Z);

  if (is.null(Z)) {
    resid.XonZ <- X
    resid.YonZ <- Y
  } else {
    if (is.null(resid.XonZ)) {
      resid_func <- function(V) comp.resids(V, Z, regr.pars, regr.method)
      if (is.matrix(X)) {
        # For KRR this approach should be much faster as the kernel matrix doesn't change
        # for each of the regressions
        resid.XonZ <- apply(X, 2, resid_func)
      } else {
        resid.XonZ <- resid_func(X)
      }
    }
    if (is.null(resid.YonZ)) {
      resid_func <- function(V) comp.resids(V, Z, regr.pars, regr.method)
      if (is.matrix(Y)) {
        resid.YonZ <- apply(Y, 2, resid_func)
      } else {
        resid.YonZ <- resid_func(Y)
      }
    }
  }
  nn <- NA
  if (NCOL(resid.XonZ) > 1 || NCOL(resid.YonZ) > 1) {
    d_X <- NCOL(resid.XonZ); d_Y <- NCOL(resid.YonZ)
    #n <- nrow(resid.XonZ)
    nn <- NROW(resid.XonZ)
    #R_mat <- rep(resid.XonZ, times=d_Y) * as.numeric(resid.YonZ[, rep(seq_len(d_X), each=d_Y)])
    R_mat <- rep(resid.XonZ, times=d_Y) * as.numeric(as.matrix(resid.YonZ)[, rep(seq_len(d_Y), each=d_X)])
    dim(R_mat) <- c(nn, d_X*d_Y)
    R_mat <- t(R_mat)
    #R_mat < R_mat / sqrt((rowMeans(R_mat^2) - rowMeans(R_mat)^2))
    R_mat <- R_mat / sqrt((rowMeans(R_mat^2) - rowMeans(R_mat)^2))

    # there are faster approaches if nsim > nn

    test.statistic <- max(abs(rowMeans(R_mat))) * sqrt(nn)
    test.statistic.sim <- apply(abs(R_mat %*% matrix(rnorm(nn*nsim), nn, nsim)), 2, max) / sqrt(nn)
    p.value <- (sum(test.statistic.sim >= test.statistic)+1) / (nsim+1)

    #plotting
    if(plot.residuals){
      par(mfrow = c(NCOL(resid.XonZ), NCOL(resid.YonZ)))
      for(ii in 1:NCOL(resid.XonZ)){
        for(jj in 1:NCOL(resid.YonZ)){
          plot(resid.XonZ[,ii], resid.YonZ[,jj], main = "scatter plot of residuals")
        }
      }
      par(mfrow = c(1,1))
    }
  } else {
    nn <- ifelse(is.null(dim(resid.XonZ)), length(resid.XonZ), dim(resid.XonZ)[1])
    R <- resid.XonZ * resid.YonZ
    R.sq <- R^2
    meanR <- mean(R)
    test.statistic <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
    p.value <- 2 * pnorm(abs(test.statistic), lower.tail = FALSE)
    if(plot.residuals){
      plot(resid.XonZ, resid.YonZ, main = "scatter plot of residuals")
    }
  }

  return(list(p.value = p.value, test.statistic = test.statistic, reject = (p.value < alpha)) )

}

comp.resids <- function(V, Z, regr.pars, regr.method) {
  V <- as.numeric(V)
  switch(regr.method,
         "gam" = {
           mod.VonZ <- train.gam(Z, V, pars = regr.pars)
         },
         "xgboost" = {
           mod.VonZ <- train.xgboost(Z, V, pars = regr.pars)
         },
         "kernel.ridge" = {
           mod.VonZ <- train.krr(Z, V, pars = regr.pars)
           #         },
           #         "nystrom" = {
           #           mod.VonZ <- train.nystrom(Z, V)
         }
  )
  return(mod.VonZ$residuals)
}



