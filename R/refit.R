lvnetRefit <- function(
  lvnetObject, # lvnet object
  data, # New covariance matrix
  sampleSize # New sample size
){
  # If lvnetLasso or lvnetSearch, extract:
  if (is(lvnetObject,"lvnetLasso") | is(lvnetObject,"lvnetSearch")){
    lvnetObject <- lvnetObject$best
  }
  
  # Check if proper:
  if (!is(lvnetObject,"lvnet")){
    stop("Input must be a 'lvnet' object.")
  }
  
  # If dataset, compute cov:
  if (nrow(data) > ncol(data)){
    if (missing(sampleSize)){
      sampleSize <- nrow(data)
    }
    data <- cov(data, use = "pairwise.complete.obs")
  }
  
  # Check if samplesize is missing:
  if (missing(sampleSize)){
    stop("sampleSize may not be missing")
  }
  
  # Extract inverse var-cov:
  sigmaInverse <- corpcor::pseudoinverse(lvnetObject$matrices$sigma_positive)
  
  # Obtain fit measures:
  Res <- ggmFit(covMat = data, # sample variance-covariance matrix
                sampleSize = sampleSize, # Sample sample-size
                invSigma = sigmaInverse, # inverse variance covariance matrix instead of pcor
                refit = FALSE,ebicTuning = lvnetObject$fitMeasures$ebicTuning,
                nPar = lvnetObject$fitMeasures$npar)
  
  return(Res)
}