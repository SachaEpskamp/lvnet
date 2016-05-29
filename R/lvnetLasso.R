# Fun:
lassoSearchFun <- function(i, tuning, Init, args, verbose, lassoMatrix,nTuning,lassoTol){

  if (verbose){
    cat(paste("\rIteration",i, "of ",nTuning))
    # flush.console() 
  }
  
  t <- tuning[i]
  args$lasso <- t
  
  Res <- list()
  
  # Make output:
  Res <- list(tuning=t)
  
  # Fit model:
  Res$res <- suppressWarnings(
      do.call(lvnet,c(args,list(fitInd = Init$mxResults$independence,
                                fitSat = Init$mxResults$saturated,
                                startValues = Init))))
        
        #data,lassoMatrix=lassoMatrix,lassoTol=lassoTol,lasso=tuning,

  
  # Extract fit indices:
  Res$fit <- Res$res$fitMeasures
  
  # Count free parameters in the model matrices:
  Res$nPar <- sum(sapply(lassoMatrix, function(m){
    mat <- Res$res$matrices[[m]]    
    if (m %in% c("lambda","beta")){
      ix <- matrix(TRUE,nrow(mat),ncol(mat))
    } else {
      ix <- upper.tri(mat,diag=FALSE)
    }
    sum(abs(mat[ix]) > lassoTol)
  }))
  
  return(Res)
}

### Search function using mutiple lasso estimates:
lvnetLasso <- function(
  data, # Data to use
  lassoMatrix, # vector of matrices to apply LASSO to 
  lassoTol = 1e-4,
  nTuning = 20,
  tuning.min = 0.01,
  tuning.max = 0.5,
  criterion = c("bic","aic","ebic"),
  verbose = TRUE,
  refit = TRUE,
  nCores = 1, # Set to > 1 to use parallel computing
  ... # lvnet arguments
){
  criterion <- match.arg(criterion)
  criterion <- switch(criterion,
                      bic = "bic",
                      aic = "aic",
                      ebic = "ebic")
  # Full results list:
  Results <- list()
  
  # Tuning sequence:
  tuning = exp(seq(log(tuning.min), log(tuning.max), length = nTuning))
  
  # Fit inital model to obtain start values and ind/sat:
  if (verbose){
    cat("Fitting initial model to obtain start-values and independence/saturated models.\n")
  }
  
  if (missing(lassoMatrix)){
    stop("'lassoMatrix' must be assigned")
  }
  if (length(lassoMatrix) > 1){
    warning("Multiple matrices in LASSO is not recommended. Use sequential estimation (e.g., first, omega_theta, then omega_psi).")
  }

  Init <- suppressWarnings(lvnet(data,...))
  
  ###
  args <- c(list(data=data,lassoMatrix=lassoMatrix,lassoTol=lassoTol),list(...))
  
  ### MAIN LOOP ###
  if (nCores == 1){
    
    Results <- lapply(seq_len(nTuning),lassoSearchFun,tuning=tuning, Init=Init, args=args, verbose=verbose, lassoMatrix=lassoMatrix,nTuning=nTuning,lassoTol=lassoTol)
  
  } else {
    # Number of clusters:
    nClust <- nCores - 1
    
    # Make cluster:
    cl <- makePSOCKcluster(nClust)  
    
    if (verbose){
      cat("Estimating LASSO penalized models.\n")
    }
    
    # Run loop:
    Results <- parallel::parLapply(cl,seq_len(nTuning),lassoSearchFun,tuning=tuning, Init=Init, args=args, verbose=FALSE, lassoMatrix=lassoMatrix,nTuning=nTuning,lassoTol=lassoTol)
    
    # Stop cluster:
    stopCluster(cl)
  }


  # Create fit table:
  Fits <- as.data.frame(do.call(rbind,lapply(Results,function(x)unlist(x[['fit']]))))
    
  # Select best:
  if (!criterion %in% names(Fits)){
    stop("Criterion is not supported")
  }
  best <- which.min(Fits[[criterion]][Fits[["df"]]>=0])
    
  if (length(best) == 0){
    stop("No identified model found.")
  }
  
  # Refit best model without lasso:
  dots <- args <- list(...)
  
  if (refit){
    newMod <- lapply(lassoMatrix, function(m){
      mat <- Results[[best]]$res$matrices[[m]]
      ifelse(abs(mat) > lassoTol,NA,0)
    })
    names(newMod) <- lassoMatrix

    for (i in seq_along(lassoMatrix)){
      dots[[lassoMatrix[[i]]]] <- newMod[[lassoMatrix[[i]]]] 
    }
    
    bestModel <- do.call(lvnet,c(list(data=data,
          fitInd = Init$mxResults$independence,
          fitSat = Init$mxResults$saturated,
          startValues = Init),
          dots))

  } else {
    bestModel <- Results[[best]]$res
  }
  
  Output <- list(
   best = bestModel,
   modList = lapply(Results,'[[','res'),
   tuning = sapply(Results,'[[','tuning'),
   lassoMatrix = lassoMatrix,
   args = c(list(data=data,
                 fitInd = Init$mxResults$independence,
                 fitSat = Init$mxResults$saturated,
                 startValues = Init),
            args),
   criterion = criterion,
   bestID = best
  )
  
  class(Output) <- "lvnetLasso"
  
  return(Output)
}