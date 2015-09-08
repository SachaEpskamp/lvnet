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

rimSearch <- function(
  matrix = c("omega_theta","omega_psi","theta","psi"), # Matrix to optimize
  criterion = c("chisq", "BIC", "AIC"), # Chisquare will attempt to remove edge with no sig difference, and otherwise add edge with sig difference.
  start = c("Za","empty","full","lvglasso","glasso"),
  alpha = 0.05,
  lambda,
  covmat,
  sampleSize,
  maxIter,
  ..., # Arguments sent to rim
  lvglassoArgs = list(gamma = 0, nLambda = 20), # Arguments sent to EBIClvglasso
  glassoArgs = list(gamma = 0, nlambda = 100), # Arguments sent to EBICglasso
  verbose = TRUE,
  file, # If not missing, reads file to continue and stores results to file.
  startValues = list()
){
  matrix <- match.arg(matrix)
  criterion <- toupper(match.arg(criterion))
  start <- match.arg(start)

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
  modList <- list()
  
  # Compute first model:
  rimArgs <- list(...)
  rimArgs$data <- covmat
  rimArgs$sampleSize <- sampleSize
  rimArgs$lambda <- lambda
  rimArgs[[matrix]] <- curMat2modMat(curMat, matrix)
  rimArgs$startValues <- startValues
  
  if (verbose){
    message("Estimating initial RIM model")
  }

  curMod <- do.call("rim", rimArgs)
  it <- 0
  
  rimArgs$fitInd <- curMod$mxResults$independence
  rimArgs$fitSat <- curMod$mxResults$saturated
  
  upTriElements <- which(upper.tri(curMat, diag=FALSE), arr.ind=TRUE)
  
  repeat{
    curEst <- curMod$matrices[[matrix]]
    it <- it + 1  
    if (it > maxIter){
      warning("Maximum number of iterations reached")
      break
    }
    modList <- c(modList,list(curMod))
    if (!missing(file)){
      save(modList,it,curMod,rimArgs,curMat,file=file)
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
      
      rimArgs[[matrix]] <- curMat2modMat(propMat, matrix)
      rimArgs$startValues[[matrix]] <- curEst * propMat
      propModels[[i]] <- do.call("rim", rimArgs)
      
      
      if (verbose){
        setTxtProgressBar(pb, i)
      }      
    }
    if (verbose) close(pb)
    
    # Create table:
    origFit <- anova(curMod)[-1,,drop=FALSE]
    fits <- do.call(rimCompare,propModels)[-1,,drop=FALSE]
    

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
      modList = modList,
      niter = it)
    
    class(Results) <- c("rimSearch","list")
    return(Results)
}



# ## This function searches fo residual interactions/correlations given a starting structure
# 
# rimSearch <- function(
#   data, # Raw data or a covariance matrix
#   lambda, # Lambda design matrix. NA indicates free parameters. If missing and psi is missing, defaults to identity matrix with warning
#   beta, # Structural matrix. If missing, defaults to zero.
#   omega_theta, # Observed residual network. If missing, defaults to matrix of zeroes
#   delta_theta, # Scaling matrix, can be missing
#   omega_psi, # Latent residual network. If missing, defaults to matrix of zeroes
#   delta_psi, # Scaling matrix, can be missing
#   psi, # Latent variance-covariance matrix. If missing, defaults to free
#   theta, # Used if model = "sem". Defaults to diagonal
#   sampleSize,
#   model = c("rim","sem"),
#   method = c(
#     "chisq", # will test for significance and stop if no significant improve can be found
#     "bic", # Will minimize bic
#     "aic"), # Will minimize aic
#   alpha = 0.05,
#   verbose = TRUE
# ){
#   if (missing(omega_theta)){
#     omega_theta <- matrix(0, ncol(data), ncol(data))
#   }
#   if (missing(theta)){
#     theta <- diag(NA, ncol(data))
#   }
#   
#   if (model[[1]]=="rim"){
#     optMat <- omega_theta
#   } else {
#     optMat <- theta
#   }
#   
#   modList <- list()
#   
#   # Compute first model:
#   curMod <- rim(data=data,lambda=lambda,omega_theta=omega_theta,omega_psi = omega_psi,psi=psi,beta=beta,delta_theta=delta_theta,delta_psi=delta_psi,theta=theta,
#                 sampleSize=sampleSize,model=model)
#   it <- 0
#   
#   repeat{
#     it <- it + 1  
#     modList <- c(modList,list(curMod))
#     
#     # Compute proposals:
#     proposals <- which(!is.na(optMat) & lower.tri(optMat,diag=FALSE), arr.ind=TRUE)
#     if (nrow(proposals)==0){
#       break
#     }
#     
#     # Proposed models:
#     propModels <- list()
#     
#     if (verbose){
#       message(paste("Iteration:",it))
#       pb <- txtProgressBar(0, nrow(proposals), style = 3)
#     }
#     
#     for (i in seq_len(nrow(proposals))){
#       mat <- optMat
#       mat[proposals[i,1],proposals[i,2]] <- mat[proposals[i,2],proposals[i,1]] <- NA
#       if (model[[1]]=="rim"){
#         propModels[[i]] <- rim(data=data,lambda=lambda,omega_theta=mat,delta_theta=delta_theta,delta_psi=delta_psi,psi=psi,beta=beta,theta=theta,
#                                sampleSize=sampleSize,model=model)
#       } else {
#         propModels[[i]] <- rim(data=data,lambda=lambda,omega_theta=omega_theta,omega_psi=omega_psi,delta_theta=delta_theta,delta_psi=delta_psi,psi=psi,beta=beta,theta=mat,
#                                sampleSize=sampleSize,model=model)
#       }
#       
#       if (verbose){
#         setTxtProgressBar(pb, i)
#       }
#     }
#     if (verbose) close(pb)
#     
#     curTab <- rimCompare(curMod)
#     propTab <- do.call(rimCompare,propModels)[-1,]
#     
#     # Optimize:
#     if (method[[1]] == "chisq"){
#       chisqDiff <- curTab$Chisq[[2]] - propTab$Chisq
#       dfDiff <- curTab$Df[[2]] - propTab$Df
#       pvals <- pchisq(chisqDiff,dfDiff,lower.tail=FALSE)
#       
#       if (!any(pvals < alpha)){
#         break
#       } else {
#         best <- which.min(pvals)
#       }
#       
#     } else  if (method[[1]] == "bic"){
#       bics <- propTab$BIC
#       
#       if (!any(bics < curTab$BIC[[2]])){
#         break
#       } else {
#         best <- which.min(bics)
#       } 
#     } else  if (method[[1]] == "aic"){
#       aics <- propTab$AIC
#       
#       if (!any(aics < curTab$AIC[[2]])){
#         break
#       } else {
#         best <- which.min(aics)
#       }
#     } else stop(paste("Method",method,"not supported."))
#     
#     optMat[proposals[best,1],proposals[best,2]] <- optMat[proposals[best,2],proposals[best,1]] <- NA
#     if (model[[1]]=="rim"){
#       omega_theta <- optMat
#     } else {
#       theta <- optMat
#     }
#     curMod <- propModels[[best]]
#   }
#   
#   Results <- list(
#     best = curMod,
#     modList = modList,
#     niter = it)
#   
#   class(Results) <- c("rimSearch","list")
#   return(Results)
# }