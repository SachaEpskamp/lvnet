# Some helper functions:
curMat2modMat <- function(x, matrix){
  x <- ifelse(x,NA,0)
  if (grepl("omega",matrix)){
    diag(x) <- 0
  } else {
    if (matrix=="psi"){
      diag(x) <- 1
    } else {
      diag(x) <- NA
    }
  }
  return(x)
}

maxNull <- function(x){
  if (length(x[!is.na(x)])==0) return(0) else return(max(x,na.rm=TRUE))
}

# Search logic: TO DOCUMENT: 
# Start with initial maxChange
# Change min(all, maxChange) improving edges. 
# Set maxChange to max(# changed edges - 1,1)
# repeat until convergence

lvnetSearch <- function(
  data,
  matrix = c("omega_theta","omega_psi","theta","psi"), # Matrix to optimize
  criterion = c("bic", "ebic","chisq","aic"), # Chisquare will attempt to remove edge with no sig difference, and otherwise add edge with sig difference.
  start = c("default","empty","full"), # CAN ALSO BE MATRIX "glasso" & "lvglasso" currently disabled. glasso runs glasso on Psi or misfit, after running CFA
  alpha = 0.05,
  lambda,
  sampleSize,
  maxIter,
  nCores = 1, # Set to > 1 to use parallel computing
  maxChange = 1, # Set by default to 1 if start = "empty" and all possible edges if start = "full".
  ..., # Arguments sent to lvnet
  # lvglassoArgs = list(gamma = 0, nRho = 20), # Arguments sent to ebiclvglasso
  # glassoArgs = list(gamma = 0.5, nlambda = 100), # Arguments sent to ebicglasso
  verbose = TRUE,
  file, # If not missing, reads file to continue and stores results to file.
  startValues = list()
){
  if (!is.matrix(start)){
    start <- start[1]
    stopifnot(start %in% c("default","empty","full"))
  }
  matrix <- match.arg(matrix)
  criterion <- match.arg(criterion)
  
  if (criterion %in% c("aic","bic","ebic")){
    criterion <- toupper(criterion)
  }
  
  # start <- match.arg(start)
  
  if (ncol(data) == nrow(data) && isSymmetric(unname(data))){
    if (missing(sampleSize)){
      stop("sampleSize needs to be assigned if input is covariance matrix.")
    }
    
    covmat <- data * (sampleSize - 1)/sampleSize
    rownames(covmat) <- colnames(covmat)
    
  } else {
    sampleSize <- nrow(data)   
    
    data <- as.matrix(data)
    covmat <- cov(data, use = "pairwise.complete.obs")* (sampleSize - 1)/sampleSize
  }
  
  if (missing(lambda)){
    message("Fitting network without latent variables")
    lambda <- matrix(,ncol(covmat),0)
  }
  
  if (is.matrix(start)){
    curMat <- start 
  } else {
    if (start == "lvglasso" & matrix != "omega_theta"){
      stop("start = 'lvglasso' only supports matrix = 'omega_theta'")
    }
    
    Nvar <- nrow(lambda)
    Nlat <- ncol(lambda)
    
    # Select start:
    if(start=="default"){
      #     if (matrix=="omega_theta"){
      #       start <- "empty"
      # #       if (Nlat > 0){
      # #         start <- "lvglasso"
      # #       } else {
      # #         start <- "glasso"
      # #       }
      #     } else {
      if (matrix %in% c("psi","omega_psi")){
        start <- "full"
      } else start <- "empty"
      # }
    }
    
    if (matrix %in% c("omega_theta","theta")){
      
      curMat <- matrix(start == "full", Nvar, Nvar)
      
    } else {
      
      curMat <- matrix(start == "full", Nlat, Nlat)
    }
  }

  
  #   if (start == "lvglasso"){
  #     stop("'lvglasso' start not supported")
  #     if (verbose){
  #       message("Estimating optimal lvglasso result")
  #     }
  #     
  #     lvglassoRes <- do.call("ebiclvglasso", c(list(S=covmat, n = sampleSize, nLatents = Nlat), lvglassoArgs  ))
  #     startValues$omega_theta <- lvglassoRes$omega_theta
  #     curMat <- lvglassoRes$omega_theta!=0
  #     
  #   } else if (start == "glasso"){
  #     browser()
  #     
  #     if (matrix %in% c("theta","psi")){
  #       stop("'glasso' is not a valid start for optimizing theta or psi.")
  #     }
  #     
  #     if (verbose){
  #       message("Estimating starting matrix")
  #     }
  #     
  #     # Start args:
  #     lvnetArgs_start <- list(...)
  #     lvnetArgs_start$data <- covmat
  #     lvnetArgs_start$sampleSize <- sampleSize
  #     lvnetArgs_start$lambda <- lambda
  #     lvnetArgs_start$startValues <- startValues
  #     
  #     if (verbose){
  #       message("Estimating initial lvnet model")
  #     }
  #     
  # 
  #     if (matrix == "omega_theta"){
  #       lvnetArgs_start[[matrix]] <- curMat2modMat(matrix(0, Nvar, Nvar), matrix)
  #       startMod <- do.call("lvnet", lvnetArgs_start)
  #       
  #       # Misfit:
  #       misFit <- startMod$sampleStats$covMat - startMod$matrices$sigma + startMod$matrices$theta
  #       
  #       # make positive definite:
  #       if (any(eigen(misFit)$values < 0)){
  #         misFit <- misFit - diag(min(eigen(misFit)$values) - 0.0001, ncol(misFit))
  #       }
  #       
  #       # Run glasso:
  #       glassoRes <- do.call(qgraph::ebicglasso, c(list(S=misFit, n = sampleSize), glassoArgs ))
  #       
  #     }
  #     
  #     
  # 
  #     
  #     glassoRes <- do.call(qgraph::ebicglasso, c(list(S=covmat, n = sampleSize), glassoArgs ))
  #     startValues$omega_theta <- glassoRes
  #     curMat <- glassoRes!=0
  #     
  #   } else 
  

  
  if (missing(maxIter)) maxIter <- ncol(curMat) * (ncol(curMat)-1) / 2
  # Empty model list:
  # modList <- list()
  
  # Compute first model:
  lvnetArgs <- list(...)
  lvnetArgs$data <- covmat
  lvnetArgs$sampleSize <- sampleSize
  lvnetArgs$lambda <- lambda
  lvnetArgs[[matrix]] <- curMat2modMat(curMat, matrix)
  lvnetArgs$startValues <- startValues
  
  if (verbose){
    message("Estimating initial lvnet model")
  }
  
  curMod <- do.call("lvnet", lvnetArgs)
  it <- 0

  if (missing(maxChange)){
    if (all(curMat)){
      maxChange <- Inf
    } else {
      maxChange <- 1
    }
  }
  
  lvnetArgs$fitInd <- curMod$mxResults$independence
  lvnetArgs$fitSat <- curMod$mxResults$saturated
  
  upTriElements <- which(upper.tri(curMat, diag=FALSE), arr.ind=TRUE)
  
  repeat{
    curEst <- curMod$matrices[[matrix]]
    it <- it + 1  
    if (it > maxIter){
      warning("Maximum number of iterations reached")
      break
    }
    # modList <- c(modList,list(curMod))
    if (!missing(file)){
      save(it,curMod,lvnetArgs,curMat,file=file)
    }
    
    propModels <- vector("list", nrow(upTriElements))
    
    ### Now using lapply:
    if (nCores == 1){
      if (verbose){
        # message(paste("Iteration:",it))
        pb <- txtProgressBar(0, nrow(upTriElements), style = 3)
      }
      
      for (i in seq_len(nrow(upTriElements))){
        propMat <- curMat
        propMat[upTriElements[i,1],upTriElements[i,2]] <- propMat[upTriElements[i,2],upTriElements[i,1]] <- 
          !curMat[upTriElements[i,1],upTriElements[i,2]]
        
        lvnetArgs[[matrix]] <- curMat2modMat(propMat, matrix)
        lvnetArgs$startValues[[matrix]] <- curEst * propMat
        propModels[[i]] <- try(do.call("lvnet", lvnetArgs),silent = TRUE)
        if (is(propModels[[i]],"try-error")){
          propModels[[i]] <- curMod
          propModels[[i]]$fitMeasures[c("aic","bic","ebic")] <- Inf
          if (curMat[upTriElements[i,1],upTriElements[i,2]]){
            propModels[[i]]$fitMeasures[c("chisq")] <- 0 
          } else {
            propModels[[i]]$fitMeasures[c("chisq")] <- Inf
          }
        }
        
        
        if (verbose){
          setTxtProgressBar(pb, i)
        }      
      }
      if (verbose) close(pb)
    } else {
      if (verbose){
        message(paste("Iteration:",it))
      }
      # Number of clusters:
      nClust <- nCores - 1
      
      # Make cluster:
      cl <- makePSOCKcluster(nClust)  
      
      # Export stuff:
      clusterExport(cl, c("curMat","upTriElements","lvnetArgs","curMat2modMat",
                          "matrix","lvnet","curEst","lvnetArgs"), envir = environment())
      
      propModels <- parLapply(cl, seq_len(nrow(upTriElements)), function(i){
        propMat <- curMat
        propMat[upTriElements[i,1],upTriElements[i,2]] <- propMat[upTriElements[i,2],upTriElements[i,1]] <- 
          !curMat[upTriElements[i,1],upTriElements[i,2]]
        
        lvnetArgs[[matrix]] <- curMat2modMat(propMat, matrix)
        lvnetArgs$startValues[[matrix]] <- curEst * propMat
        propModels[[i]] <- try(do.call("lvnet", lvnetArgs),silent = TRUE)
        if (is(propModels[[i]],"try-error")){
          propModels[[i]] <- curMod
          propModels[[i]]$fitMeasures[c("aic","bic","ebic")] <- Inf
          if (curMat[upTriElements[i,1],upTriElements[i,2]]){
            propModels[[i]]$fitMeasures[c("chisq")] <- 0 
          } else {
            propModels[[i]]$fitMeasures[c("chisq")] <- Inf
          }
        }
      })
      
      # Stop cluster:
      stopCluster(cl)
    }
    

    
    # Create table:
    origFit <- anova(curMod)[-1,,drop=FALSE]
    fits <- do.call(lvnetCompare,propModels)[-1,,drop=FALSE]
    
    
    if (criterion %in% c("AIC","BIC","EBIC")){
      fits[[criterion]][is.na(fits[[criterion]])] <- Inf
      
      # Any is better?
      if (!any(fits[[criterion]] < origFit[[criterion]])){
        break
      } else {
        
#         # Which best?
#         best <- which.min(fits[[criterion]])
        # Select set of best edges:

        nImprove <- sum(fits[[criterion]] < origFit[[criterion]])
        best <- order(fits[[criterion]],decreasing=FALSE)[1:min(nImprove,maxChange)]
      }
    } else {
      # Test if parameter is currently an edge or not:
      curEdge <- curMat[upper.tri(curMat,diag=FALSE)]
      
      # Obtain the Chi-square of current model:
      curChisq <- curMod$fitMeasures$chisq
      curDF <- curMod$fitMeasures$df
      
      # Obtain the chi-squares of proposed models:
      propChisq <- fits$Chisq
      propDF <- fits$Df
      
      # Compute the p-values of chi-square difference tests:
      Pvals <- pchisq(abs(curChisq-propChisq), abs(curDF - propDF), lower.tail=FALSE)
      
      # If not currently edge, adding an edge should significantly improve fit (p < 0.05).
      # If currently an edge, removing that edge should *not* significantly worsen fit (p > 0.05)
      # Prioritize removing edges
      
      PvalsNAmax <- ifelse(is.na(Pvals),1,Pvals)
      PvalsNAmin <- ifelse(is.na(Pvals),0,Pvals)
      
      # Edges that can be removed:
      improveRemoved <- which(curEdge & PvalsNAmin > alpha)
      
      # Relative rank how well they can be removed:
      scoreRemoved <- rank(-PvalsNAmin[improveRemoved],ties.method = "random")
      
      # edges that can be added:
      improveAdded <- which(!curEdge & PvalsNAmax < alpha)
      
      
      # Relative rank how well they can be added:
      scoreAdded <- rank(PvalsNAmax[improveAdded],ties.method = "random") + maxNull(scoreRemoved)
      
      # Combine:
      improve <- c(improveRemoved,improveAdded)
      score <- c(scoreRemoved,scoreAdded)
      
      # Number that can be improved:
      nImprove <- length(improve)
      
      if (nImprove==0){
        break
      }
      best <- improve[order(score,decreasing=FALSE)][1:min(nImprove,maxChange)]
    }
    
    # Number of edges to change:
    nChange <- length(best)
    
    if (verbose){
      if (nChange > 1)  {
        message(paste("Changing",nChange,"edges"))
      } else {
        message(paste("Changing",nChange,"edge"))
      }
    }
    
    # Update the current matrix:
    for (b in best){
      curMat[upTriElements[b,1],upTriElements[b,2]] <- curMat[upTriElements[b,2],upTriElements[b,1]] <- 
        !curMat[upTriElements[b,1],upTriElements[b,2]]
      

    } 

    # Set new model:
    if (nChange > 1){
      lvnetArgs[[matrix]] <- curMat2modMat(curMat, matrix)
      lvnetArgs$startValues[[matrix]] <- curEst * curMat
      curMod <- do.call("lvnet", lvnetArgs)
    } else {
      # Compute new current model:
      curMod <- propModels[[best]]
    }

    
    # Set the maxChange counter:
    maxChange <- max(length(best)-1,1)
    
  }
  
  Results <- list(
    best = curMod,
    # modList = modList,
    niter = it,
    criterion = criterion)
  
  class(Results) <- c("lvnetSearch","list")
  return(Results)
}
