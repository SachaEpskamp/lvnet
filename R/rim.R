# Main function for confirmatory RIM

rim <- function(
  data, # Raw data or a covariance matrix
  lambda, # Lambda design matrix. NA indicates free parameters. If missing and psi is missing, defaults to identity matrix with warning
  beta, # Structural matrix. If missing, defaults to zero.
  omega_theta, # Observed residual network. If missing, defaults to matrix of zeroes
  delta_theta, # Scaling matrix, can be missing
  omega_psi, # Latent residual network. If missing, defaults to matrix of zeroes
  delta_psi, # Scaling matrix, can be missing
  psi, # Latent variance-covariance matrix. If missing, defaults to free
  theta, # Used if model = "sem". Defaults to diagonal
  sampleSize,
  fitInd,
  fitSat,
  startValues=list() # Named list of starting values
  ){
  
  Nvar <- ncol(data)
  
  ### Generate model:
  mod <- generateRIMmodel(
    data = data,
    lambda = lambda,
    omega_psi = omega_psi,
    omega_theta = omega_theta,
    delta_psi = delta_psi,
    delta_theta = delta_theta,
    psi = psi,
    beta = beta,
    theta = theta,
    sampleSize = sampleSize,
    name = "model",
    startValues=startValues)

  fitMod <- mxRun(mod, silent = TRUE,
                  suppressWarnings = TRUE)
  
  if (missing(fitSat)){
    # Saturated model:
    satMod <- generateRIMmodel(
      data = data, 
      lambda = diag(Nvar), 
      psi = matrix(NA,Nvar,Nvar), 
      theta = matrix(0, Nvar,Nvar), 
      name = "saturated",
      sampleSize = sampleSize
    )
    
    
    fitSat <- mxRun(satMod, silent = TRUE,
                    suppressWarnings = TRUE)  
  }
  
  if (missing(fitInd)){
    # Independence model:
    indMod <- generateRIMmodel(
      data = data, 
      lambda = diag(Nvar), 
      psi = diag(NA, Nvar, Nvar), 
      theta = matrix(0, Nvar, Nvar), 
      name = "independence",
      sampleSize = sampleSize
    )
    
    fitInd <- mxRun(indMod, silent = TRUE,
                    suppressWarnings = TRUE)
  }

  if (missing(sampleSize)){
    sampleSize <- nrow(data)
  }
  
  # Estract estimated matrices:
  Matrices <- c(lapply(fitMod$matrices,'slot','values'),
  lapply(fitMod$algebras,'slot','result'))
  
  ### COMPUTE RESULTS ###
  Results <- list(
    matrices = Matrices,
    sampleStats = list(
      covMat = fitMod$data@observed,
      sampleSize = sampleSize
    ),
    mxResults = list(
      model = fitMod,
      independence = fitInd,
      saturated = fitSat),
    fitMeasures = list()
  )
  
  sigma <- Results$matrices$sigma
  S <- Results$sampleStats$covMat
  
  # Compute chi-square:
  Results$fitMeasures$npar <- summary(fitMod)$estimatedParameters
  Results$fitMeasures$fmin <- (sum(diag(S %*% solve(sigma)))- log(det(S %*% solve(sigma))) - Nvar)/2
  Results$fitMeasures$chisq <- 2 * (sampleSize - 1) * Results$fitMeasures$fmin
  Results$fitMeasures$df <- summary(fitMod)$degreesOfFreedom
  Results$fitMeasures$pvalue <- pchisq(Results$fitMeasures$chisq, Results$fitMeasures$df, lower.tail = FALSE)
  
  # Baseline model:
  sigmaBase <- fitInd$algebras$sigma$result
  Results$fitMeasures$baseline.chisq <- (sampleSize - 1) * (sum(diag(S %*% solve(sigmaBase)))- log(det(S %*% solve(sigmaBase))) - Nvar)
  Results$fitMeasures$baseline.df <- summary(fitInd)$degreesOfFreedom
  Results$fitMeasures$baseline.pvalue <- pchisq(Results$fitMeasures$baseline.chisq, Results$fitMeasures$baseline.df, lower.tail = FALSE)
  
  # Incremental Fit Indices
  Tb <- Results$fitMeasures$baseline.chisq
  Tm <- Results$fitMeasures$chisq
  
  dfb <- Results$fitMeasures$baseline.df
  dfm <- Results$fitMeasures$df
  
  Results$fitMeasures$nfi <- (Tb - Tm) / Tb
  Results$fitMeasures$tli <-  (Tb/dfb - Tm/dfm) / (Tb/dfb - 1) 
  Results$fitMeasures$rfi <-  (Tb/dfb - Tm/dfm) / (Tb/dfb ) 
  Results$fitMeasures$ifi <-  (Tb - Tm) / (Tb - dfm)
  Results$fitMeasures$rni <-  ((Tb- dfb) - (Tm - dfm)) / (Tb - dfb)
  Results$fitMeasures$cfi <- ifelse(dfm > Tm, 1, 1 - (Tm - dfm)/(Tb - dfb))
  
  # RMSEA
  Results$fitMeasures$rmsea <- sqrt( (Tm - dfm) / ((sampleSize - 1) * dfm))
  
  # Codes for rmsea confidence interval taken from lavaan:
  lower.lambda <- function(lambda) {
    (pchisq(Tm, df=dfm, ncp=lambda) - 0.95)
  }
  if(is.na(Tm) || is.na(dfm)) {
    Results$fitMeasures$rmsea.ci.lower <- NA
  } else if(dfm < 1 || lower.lambda(0) < 0.0) {
    Results$fitMeasures$rmsea.ci.lower <- 0
  } else {
    lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=Tm)$root,
                    silent=TRUE)
    if(inherits(lambda.l, "try-error")) { lambda.l <- NA }
    Results$fitMeasures$rmsea.ci.lower <- sqrt( lambda.l/(sampleSize*dfm) )
  }
  
  N.RMSEA <- max(sampleSize, Tm*4) 
  upper.lambda <- function(lambda) {
    (pchisq(Tm, df=dfm, ncp=lambda) - 0.05)
  }
  if(is.na(Tm) || is.na(dfm)) {
    Results$fitMeasures$rmsea.ci.upper <- NA
  } else if(dfm < 1 || upper.lambda(N.RMSEA) > 0 || upper.lambda(0) < 0) {
    Results$fitMeasures$rmsea.ci.upper <- 0
  } else {
    sink(tempfile())
    lambda.u <- try(uniroot(f=upper.lambda, lower=0,upper=N.RMSEA)$root,
                    silent=TRUE)
    sink()
    if(inherits(lambda.u, "try-error")) {lambda.u <- NA }
    
    Results$fitMeasures$rmsea.ci.upper <- sqrt( lambda.u/(sampleSize*dfm) )
  }
  
  Results$fitMeasures$rmsea.pvalue <- 
    1 - pchisq(Tm, df=dfm, ncp=(sampleSize*dfm*0.05^2))
  
  # RMR:
  sqrt.d <- 1/sqrt(diag(S))
  D <- diag(sqrt.d, ncol=length(sqrt.d))
  R <- D %*% (S - sigma) %*% D
  RR <- (S - sigma)
  e <- Nvar*(Nvar+1)/2 + Nvar
  
  Results$fitMeasures$rmr <- sqrt( sum(RR[lower.tri(RR, diag=TRUE)]^2) / e )
  Results$fitMeasures$srmr <-  sqrt( sum(R[lower.tri(R, diag=TRUE)]^2) / e )
  
  
  # information criteria:
  
  # Saturated log-likelihood:
  c <- sampleSize*Nvar/2 * log(2 * pi)
  satLL <- ( -c -(sampleSize/2) * log(det(S)) - (sampleSize/2)*Nvar )
  
  # log likelihood:
  LL <-  -sampleSize * (Results$fitMeasures$fmin -  satLL/sampleSize)
  
  Results$fitMeasures$logl <- LL
  Results$fitMeasures$unrestricted.logl <- satLL
  
  Results$fitMeasures$aic <-  -2*LL + 2* Results$fitMeasures$npar
  
  BIC <- -2*LL + Results$fitMeasures$npar * log(sampleSize)
  Results$fitMeasures$bic <- BIC
  
  # add sample-size adjusted bic
  N.star <- (sampleSize + 2) / 24
  BIC2 <- -2*LL + Results$fitMeasures$npar * log(N.star)
  Results$fitMeasures$bic2 <- BIC2
  
  
  
  class(Results) <- "rim"
  
  return(Results)
}