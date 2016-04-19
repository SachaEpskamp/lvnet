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

lvnetSearch <- function(
  data,
  matrix = c("omega_theta","omega_psi","theta","psi"), # Matrix to optimize
  criterion = c("chisq", "BIC", "AIC"), # Chisquare will attempt to remove edge with no sig difference, and otherwise add edge with sig difference.
  start = c("default","empty","full","lvglasso","glasso"),
  alpha = 0.05,
  lambda,
  sampleSize,
  maxIter,
  ..., # Arguments sent to lvnet
  lvglassoArgs = list(gamma = 0, nRho = 20), # Arguments sent to EBIClvglasso
  glassoArgs = list(gamma = 0, nlambda = 100), # Arguments sent to EBICglasso
  verbose = TRUE,
  file, # If not missing, reads file to continue and stores results to file.
  startValues = list()
){
  matrix <- match.arg(matrix)
  criterion <- toupper(match.arg(criterion))
  start <- match.arg(start)
  
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
  
  if (start == "lvglasso" & matrix != "omega_theta"){
    stop("start = 'lvglasso' only supports matrix = 'omega_theta'")
  }
  
  Nvar <- nrow(lambda)
  Nlat <- ncol(lambda)
  
  # Select start:
  if(start=="default"){
    if (matrix=="omega_theta"){
      if (Nlat > 0){
        start <- "lvglasso"
      } else {
        start <- "glasso"
      }
    } else {
      if (matrix %in% c("psi","omega_psi")){
        start <- "full"
      } else start <- "empty"
    }
  }
  
  if (start == "lvglasso"){
    
    if (verbose){
      message("Estimating optimal lvglasso result")
    }
    
    lvglassoRes <- do.call("EBIClvglasso", c(list(S=covmat, n = sampleSize, nLatents = Nlat), lvglassoArgs  ))
    startValues$omega_theta <- lvglassoRes$omega_theta
    curMat <- lvglassoRes$omega_theta!=0
    
  } else if (start == "glasso"){
    
    if (verbose){
      message("Estimating optimal glasso result")
    }
    
    glassoRes <- do.call(qgraph::EBICglasso, c(list(S=covmat, n = sampleSize), glassoArgs ))
    startValues$omega_theta <- glassoRes
    curMat <- glassoRes!=0
    
  } else if (matrix %in% c("omega_theta","theta")){
    
    curMat <- matrix(start == "full", Nvar, Nvar)
    
  } else {
    
    curMat <- matrix(start == "full", Nlat, Nlat)
  }
  
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
    
    if (verbose){
      message(paste("Iteration:",it))
      pb <- txtProgressBar(0, nrow(upTriElements), style = 3)
    }
    
    for (i in seq_len(nrow(upTriElements))){
      propMat <- curMat
      propMat[upTriElements[i,1],upTriElements[i,2]] <- propMat[upTriElements[i,2],upTriElements[i,1]] <- 
        !curMat[upTriElements[i,1],upTriElements[i,2]]
      
      lvnetArgs[[matrix]] <- curMat2modMat(propMat, matrix)
      lvnetArgs$startValues[[matrix]] <- curEst * propMat
      propModels[[i]] <- do.call("lvnet", lvnetArgs)
      
      
      if (verbose){
        setTxtProgressBar(pb, i)
      }      
    }
    if (verbose) close(pb)
    
    # Create table:
    origFit <- anova(curMod)[-1,,drop=FALSE]

    fits <- do.call(lvnetCompare,propModels)[-1,,drop=FALSE]
    

    if (criterion %in% c("AIC","BIC")){
      fits <- fits[rowSums(is.na(fits))==0,]
      if (nrow(fits)==0) break
      
      # Any is better?
      if (!any(fits[[criterion]] < origFit[[criterion]])){
        break
      } else {
        
        # Which best?
        best <- which.min(fits[[criterion]])

        
      }
    } else {

      # Significance testing!
      curEdge <- curMat[upper.tri(curMat,diag=FALSE)]

      curChisq <- curMod$fitMeasures$chisq
      curDF <- curMod$fitMeasures$df
      
      propChisq <- fits$Chisq
      propDF <- fits$Df
      # First try to add an edge that improves fit significantly and the best:
      Pvals <- pchisq(abs(curChisq-propChisq), abs(curDF - propDF), lower.tail=FALSE)
   
      PvalsNAmax <- ifelse(is.na(Pvals),1,Pvals)
      PvalsNAmin <- ifelse(is.na(Pvals),0,Pvals)
      if (any(PvalsNAmax < alpha & !curEdge)){
        best <- which(PvalsNAmax == min(PvalsNAmax[!curEdge]))[1]
      } else if (any(PvalsNAmin > alpha & curEdge)){
        best <- which(PvalsNAmin == max(PvalsNAmin[curEdge]))[1]
      } else {
        break
      }
    }

    curMat[upTriElements[best,1],upTriElements[best,2]] <- curMat[upTriElements[best,2],upTriElements[best,1]] <- 
      !curMat[upTriElements[best,1],upTriElements[best,2]]
    curMod <- propModels[[best]]
    
  }
  
    Results <- list(
      best = curMod,
      # modList = modList,
      niter = it)
    
    class(Results) <- c("lvnetSearch","list")
    return(Results)
}
