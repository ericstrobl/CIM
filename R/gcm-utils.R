
require(mgcv)

train.krr <- function(X,y,pars = list()){
  if(!exists("lambda",pars)){
    pars$lambda <- NA
  }
  if(!exists("sigma",pars)){
    pars$sigma <- 1
  }
  if(!exists("sigmas.CV",pars)){
    pars$sigmas.CV <- NA
  }
  if(!exists("lambdas.CV",pars)){
    pars$lambdas.CV <- 10^(-6:4)
  }


  if(is.null(X) || dim(as.matrix(X))[2] == 0){
    result <- list()
    result$Yfit <- as.matrix(rep(mean(y), length(y)))
    result$residuals <- as.matrix(y - result$Yfit)
    result$model <- NA
    result$df <- NA
    result$edf <- NA
    result$edf1 <- NA
    result$p.values <- NA
  } else {
    data.set <- list(x = X, y = y)
    class(data.set) <- "CVST.data"
    if(is.na(pars$sigma)){
      sigmas <- pars$sigmas.CV
    } else {
      sigmas <- pars$sigma
    }
    if(is.na(pars$lambda)){
      lambdas <- pars$lambdas.CV
    } else {
      lambdas <- pars$lambda
    }
    krr = constructKRRLearner()
    params <- constructParams(kernel="rbfdot", sigma=sigmas, lambda=lambdas)
    opt <- CV(data.set, krr, params, fold=10, verbose=FALSE)
    # show the lambda that was chosen by CV
    # show(opt[[1]]$lambda)
    param <- list(kernel="rbfdot", sigma=opt[[1]]$sigma, lambda=opt[[1]]$lambda)
    krr.model <- krr$learn(data.set,param)
    pred <- krr$predict(krr.model,data.set)
    result <- list()
    result$Yfit <- pred
    result$residuals <- as.matrix(y - result$Yfit)
    result$model <- NA
    result$df <- NA
    result$edf <- NA
    result$edf1 <- NA
    result$p.values <- NA
  }
  return(result)
}



train.xgboost <- function(X,y,pars = list()){
  n <- length(y)
  if(!exists("nrounds",pars)){
    #pars$nrounds <- round(3*log(length(y))/log(10))
    #pars$nrounds <- max(5,round(sqrt(length(y))/3))
    pars$nrounds <- 50
  }
  if(!exists("max_depth",pars)){
    #pars$max_depth <- max(5,round(sqrt(length(y))/3))
    pars$max_depth <- c(1,3,4,5,6)
    # pars$max_depth = c(1,3,5)
  }
  if(!exists("CV.folds",pars)){
    pars$CV.folds <- 10
  }
  if(!exists("ncores",pars)){
    pars$ncores <- 1
  }
  if(!exists("early_stopping",pars)){
    pars$early_stopping <- 10
  }
  if(!exists("silent",pars)){
    pars$silent <- TRUE
  }
  if(is.null(X) || dim(as.matrix(X))[2] == 0){
    result <- list()
    result$Yfit <- as.matrix(rep(mean(y), length(y)))
    result$residuals <- as.matrix(y - result$Yfit)
    result$model <- NA
    result$df <- NA
    result$edf <- NA
    result$edf1 <- NA
    result$p.values <- NA
  } else {
    X <- as.matrix(X)
    # if(!is.na(pars$CV.folds)){
    #   num.folds <- pars$CV.folds
    #   rmse <- matrix(0,pars$nrounds, length(pars$max_depth))
    #   whichfold <- sample(rep(1:num.folds, length.out = n))
    #   for(j in 1:length(pars$max_depth)){
    #     max_depth <- pars$max_depth[j]
    #     for(i in 1:num.folds){
    #       dtrain <- xgb.DMatrix(data = data.matrix(X[whichfold != i,]), label=y[whichfold != i])
    #       dtest <- xgb.DMatrix(data = data.matrix(X[whichfold == i,]), label=y[whichfold == i])
    #       watchlist <- list(train = dtrain, test = dtest)
    #       if(pars$ncores > 1){
    #         bst <- xgb.train(data = dtrain, nthread = pars$ncores, watchlist = watchlist, nrounds = pars$nrounds, max_depth = max_depth, verbose = FALSE, early_stopping_rounds = pars$early_stopping, callbacks = list(cb.evaluation.log()))
    #       } else{
    #         bst <- xgb.train(data = dtrain, nthread = 1, watchlist = watchlist, nrounds = pars$nrounds, max_depth = max_depth, verbose = FALSE, early_stopping_rounds = pars$early_stopping, callbacks = list(cb.evaluation.log()))
    #       }
    #       newscore <- (bst$evaluation_log$test_rmse)^1
    #       # print(bst$evaluation_log$test_rmse)
    #       if(length(newscore) < pars$nrounds){
    #         newscore <- c(newscore, rep(Inf, pars$nrounds - length(newscore)))
    #       }
    #       rmse[,j] <- rmse[,j] + newscore
    #     }
    #   }
    #   mins <- arrayInd(which.min(rmse),.dim = dim(rmse))
    #   if(!pars$silent){
    #     show(rmse)
    #     show(mins)
    #     if( (mins[1] == 1) | (mins[1] == pars$nrounds) | (mins[2] == 1) | (mins[2] == length(pars$max_depth)) ){
    #       show("There have been parameters selected that were the most extreme of the CV values")
    #       show(mins)
    #     }
    #   }
    #   final.nrounds <- mins[1]
    #   final.max_depth <- pars$max_depth[mins[2]]
    # } else{
    #   if(length(pars$max_depth) > 1){stop("providing a vector of parameters must be used with CV")}
    #   final.max_depth <- pars$max_depth
    #   final.nrounds <- pars$nrounds
    # }
    dtrain <- xgb.DMatrix(data = data.matrix(X), label=y)
    final.nrounds = 50; final.max_depth = 6;
    bstY <- xgb.train(data = dtrain, nrounds = final.nrounds, max_depth = final.max_depth,
                      verbose = !pars$silent)
    result <- list()
    result$Yfit <- predict(bstY, data.matrix(X))
    result$residuals <- as.matrix(y - result$Yfit)
    result$model <- NA
    result$df <- NA
    result$edf <- NA
    result$edf1 <- NA
    result$p.values <- NA
  }
  return(result)
}



train.gam <- function(X,y,pars = list()){
  if(!exists("numBasisFcts",pars)){
    pars$numBasisFcts <- 100
    # pars$numBasisFcts <- seq(5,100,length.out = 10)
  }
  if(!exists("staysilent",pars)){
    pars$staysilent <- TRUE
  }
  if(!exists("CV.folds",pars)){
    pars$CV.folds <- NA
  }
  if(is.null(X) || dim(as.matrix(X))[2] == 0){
    result <- list()
    result$Yfit <- as.matrix(rep(mean(y), length(y)))
    result$residuals <- as.matrix(y - result$Yfit)
    result$model <- NA
    result$df <- NA
    result$edf <- NA
    result$edf1 <- NA
    result$p.values <- NA
  } else {
    p <- dim(as.matrix(X))
    if(!is.na(pars$CV.folds)){
      # This is not necessary -- gam has a built-in CV
      num.folds <- pars$CV.folds
      rmse <- Inf
      whichfold <- sample(rep(1:num.folds, length.out = p[1]))
      for(j in 1:length(pars$numBasisFcts)){
        mod <- train.gam(as.matrix(X)[whichfold == j,], y[whichfold == j], pars = list(numBasisFcts = pars$numBasisFcts, CV.folds = NA))
        datframe <- data.frame(as.matrix(X)[whichfold != j,])
        names(datframe) <- paste("var", 2:p[2], sep = "")
        rmse.tmp <- sum((predict(mod$model, datframe) - y[whichfold != j])^2)
        if(rmse.tmp < rmse){
          rmse <- rmse.tmp
          final.numBasisFcts <- pars$numBasisFcts[j]
        }
      }
    } else {
      final.numBasisFcts <- pars$numBasisFcts
    }

    if(p[1]/p[2] < 3*final.numBasisFcts)
    {
      final.numBasisFcts <- ceiling(p[1]/(3*p[2]))
      if(pars$staysilent == FALSE){
        cat("changed number of basis functions to    ", final.numBasisFcts, "    in order to have enough samples per basis function\n")
      }
    }
    dat <- data.frame(as.matrix(y),as.matrix(X))
    coln <- rep("null",p[2]+1)
    for(i in 1:(p[2]+1))
    {
      coln[i] <- paste("var",i,sep="")
    }
    colnames(dat) <- coln
    labs<-"var1 ~ "
    if(p[2] > 1)
    {
      for(i in 2:p[2])
      {
        labs<-paste(labs,"s(var",i,",k = ",final.numBasisFcts,") + ",sep="")
        #      labs<-paste(labs,"s(var",i,") + ",sep="")
        # labs<-paste(labs,"lo(var",i,") + ",sep="")
      }
    }
    labs<-paste(labs,"s(var",p[2]+1,",k = ",final.numBasisFcts,")",sep="")
    # labs<-paste(labs,"s(var",p[2]+1,", bs = "cc")",sep="") # factor 2 faster
    # labs<-paste(labs,"s(var",p[2]+1,", bs = "cr")",sep="") # factor 2 + eps faster
    mod_gam <- FALSE
    try(mod_gam <- gam(formula=formula(labs), data=dat),silent = TRUE)
    if(typeof(mod_gam) == "logical")
    {
      cat("There was some error with gam. The smoothing parameter is set to zero.\n")
      labs<-"var1 ~ "
      if(p[2] > 1)
      {
        for(i in 2:p[2])
        {
          labs<-paste(labs,"s(var",i,",k = ",final.numBasisFcts,",sp=0) + ",sep="")
        }
      }
      labs<-paste(labs,"s(var",p[2]+1,",k = ",final.numBasisFcts,",sp=0)",sep="")
      mod_gam <- gam(formula=formula(labs), data=dat)
    }
    result <- list()
    result$Yfit <- as.matrix(mod_gam$fitted.values)
    result$residuals <- as.matrix(mod_gam$residuals)
    result$model <- mod_gam
    result$df <- mod_gam$df.residual
    result$edf <- mod_gam$edf
    result$edf1 <- mod_gam$edf1
    result$p.values <- summary.gam(mod_gam)$s.pv
    # for degree of freedom see mod_gam$df.residual
    # for aic see mod_gam$aic
  }
  return(result)
}


