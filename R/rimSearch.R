## This function searches fo residual interactions/correlations given a starting structure

rimSearch <- function(
  data, # Raw data or a covariance matrix
  lambda, # Lambda design matrix. NA indicates free parameters. If missing and psi is missing, defaults to identity matrix with warning
  omega, # Network. If missing, defaults to matrix of zeroes
  psi, # Latent variance-covariance matrix. If missing, defaults to free
  beta, # Structural matrix. If missing, defaults to zero.
  delta, # Scaling matrix, can be missing
  theta, # Used if model = "sem". Defaults to diagonal
  sampleSize,
  model = c("rim","sem"),
  method = c(
    "chisq", # will test for significance and stop if no significant improve can be found
    "bic", # Will minimize bic
    "aic"), # Will minimize aic
  alpha = 0.05,
  verbose = TRUE
){
  if (missing(omega)){
    omega <- matrix(0, ncol(data), ncol(data))
  }
  if (missing(theta)){
    theta <- diag(NA, ncol(data))
  }
  
  if (model[[1]]=="rim"){
    optMat <- omega
  } else {
    optMat <- theta
  }
  
  modList <- list()
  
  # Compute first model:
  curMod <- rim(data=data,lambda=lambda,omega=omega,psi=psi,beta=beta,delta=delta,theta=theta,
                sampleSize=sampleSize,model=model)
  it <- 0
  
  repeat{
    it <- it + 1  
    modList <- c(modList,list(curMod))
    
    # Compute proposals:
    proposals <- which(!is.na(optMat) & lower.tri(optMat,diag=FALSE), arr.ind=TRUE)
    if (nrow(proposals)==0){
      break
    }
    
    # Proposed models:
    propModels <- list()
    
    if (verbose){
      message(paste("Iteration:",it))
      pb <- txtProgressBar(0, nrow(proposals), style = 3)
    }
    
    for (i in seq_len(nrow(proposals))){
      mat <- optMat
      mat[proposals[i,1],proposals[i,2]] <- mat[proposals[i,2],proposals[i,1]] <- NA
      if (model[[1]]=="rim"){
        propModels[[i]] <- rim(data=data,lambda=lambda,omega=mat,psi=psi,beta=beta,delta=delta,theta=theta,
                               sampleSize=sampleSize,model=model)
      } else {
        propModels[[i]] <- rim(data=data,lambda=lambda,omega=omega,psi=psi,beta=beta,delta=delta,theta=mat,
                               sampleSize=sampleSize,model=model)
      }
      
      if (verbose){
        setTxtProgressBar(pb, i)
      }
    }
    if (verbose) close(pb)
    
    curTab <- rimCompare(curMod)
    propTab <- do.call(rimCompare,propModels)[-1,]
    
    # Optimize:
    if (method[[1]] == "chisq"){
      chisqDiff <- curTab$Chisq[[2]] - propTab$Chisq
      dfDiff <- curTab$Df[[2]] - propTab$Df
      pvals <- pchisq(chisqDiff,dfDiff,lower.tail=FALSE)
      
      if (!any(pvals < alpha)){
        break
      } else {
        best <- which.min(pvals)
      }
      
    } else  if (method[[1]] == "bic"){
      bics <- propTab$BIC
      
      if (!any(bics < curTab$BIC[[2]])){
        break
      } else {
        best <- which.min(bics)
      } 
    } else  if (method[[1]] == "aic"){
      aics <- propTab$AIC
      
      if (!any(aics < curTab$AIC[[2]])){
        break
      } else {
        best <- which.min(aics)
      }
    } else stop(paste("Method",method,"not supported."))
    
    optMat[proposals[best,1],proposals[best,2]] <- optMat[proposals[best,2],proposals[best,1]] <- NA
    if (model[[1]]=="rim"){
      omega <- optMat
    } else {
      theta <- optMat
    }
    curMod <- propModels[[best]]
  }
  
  Results <- list(
    best = curMod,
    modList = modList,
    niter = it)
  
  class(Results) <- c("rimSearch","list")
  return(Results)
}