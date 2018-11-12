# Main function for confirmatory lvnet
setSym <- function(x) {
  if (is.matrix(x)){
    return((x + t(x)) / 2)
  } else return(x)
}

countPars <- function(x,tol=sqrt(.Machine$double.eps)){
  Matrices <- x@matrices
  # Pars per matrix:
  counts <- sapply(Matrices,function(mat){
    symm <- "SymmMatrix" %in% class(mat)
    # Index (either full or UT incl diag):
    if (symm){
      ix <- upper.tri(mat@values,diag=TRUE)
    } else {
      ix <- matrix(TRUE,nrow(mat@values),ncol(mat@values))
    }
    free <- mat@free

    # sum(abs(mat@values[free & ix]) > tol)    
    # Collect labels:
    length(unique(c(mat$labels[abs(mat@values) > tol & free & ix])))
  })
  
  nPar <- sum(counts)
  nVar <- ncol(x@data@observed)
  nObs <- nVar * (nVar+1) / 2
  
  return(list(nPar=nPar,DF=nObs - nPar,nObs=nObs,nVar=nVar))
}

lvnet <- function(
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
  startValues=list(), # Named list of starting values. CAN ALSO BE lvnet OBJECT!
  scale = FALSE, # Standardize cov to cor before estimation
  nLatents, # allows for quick specification of fully populated lambda matrix.
  
  # Experimental!!!
  lasso = 0, # IF NOT 0, use lasso penalty
  lassoMatrix, # Vector of character-string of matrices to apply LASSO penalty on
  # optimizer = c("default","SLSQP","NPSOL","CSOLNP")
  lassoTol = 1e-4,
  ebicTuning = 0.5,
  mimic = c("lavaan","lvnet"),
  fitFunction = c("penalizedML","ML"),
  exogenous # Vector of exogenous variables
  
  # Optimizer:
  # nCores = 1
){
  fitFunction <- match.arg(fitFunction)
  Nvar <- ncol(data)
  mimic <- match.arg(mimic)
  
  # Check args:
  if (!missing(lambda) & !missing(nLatents)){
    warning("'nLatents' ignored if 'lambda' is also assigned.")
  }
  
#   if (nCores > 1){
#     mxOption(NULL, "Number of Threads", nCores - 1)
#   } else {
#     mxOption(NULL, "Number of Threads", NULL)
#   }
  
  #   optimizer <- match.arg(optimizer)
  #   if (optimizer=="default"){
  #     if (lasso != 0)
  #       optimizer <- "NPSOL"
  #   } else {
  #     optimizer <- "SLSQP"
  #   }
  #   
  #   # Set optimizer:
  #   mxOption(NULL,"Default optimizer",optimizer)
  
  # Check for lasso:
  if (lasso != 0){
    if (missing(lassoMatrix)){
      stop ("'lassoMatrix' must not be missing if lasso != 0")
    }
    if (any(!lassoMatrix %in% c("lambda","psi","omega_psi","theta","omega_theta","beta"))){
      stop("LASSO only supported for 'lambda', 'beta', 'psi', 'omega_psi', 'theta' and 'omega_theta'")
    }
  }
  
  # If startvalues is lvnet, add to startvalues:
  if (is(startValues,"lvnet") || lasso != 0){
  # if (is(startValues,"lvnet")){
    if (is(startValues,"lvnet")){
      initRes <- startValues
      startValues <- list()
    } else {

      initRes <- lvnet(
        data=data, # Raw data or a covariance matrix
        lambda=lambda, # Lambda design matrix. NA indicates free parameters. If missing and psi is missing, defaults to identity matrix with warning
        beta=beta, # Structural matrix. If missing, defaults to zero.
        omega_theta=omega_theta, # Observed residual network. If missing, defaults to matrix of zeroes
        delta_theta=delta_theta, # Scaling matrix, can be missing
        omega_psi=omega_psi, # Latent residual network. If missing, defaults to matrix of zeroes
        delta_psi=delta_psi, # Scaling matrix, can be missing
        psi=psi, # Latent variance-covariance matrix. If missing, defaults to free
        theta=theta, # Used if model = "sem". Defaults to diagonal
        sampleSize=sampleSize,
        fitInd=fitInd,
        fitSat=fitSat,
        startValues=startValues,
        nLatents=nLatents
      )
    }

    if (is.null(startValues[['lambda']]) && !is.null(initRes$matrices$lambda) &&  ncol(initRes$matrices$lambda)>0 && nrow(initRes$matrices$lambda) > 0 && all(is.finite(initRes$matrices$lambda))){
      startValues[['lambda']] <- initRes$matrices$lambda
      if (!missing(lambda)){
        startValues[['lambda']] <- startValues[['lambda']] * is.na(lambda)
      }
    }
    if (is.null(startValues[['beta']])&& !is.null(initRes$matrices$beta) && ncol(initRes$matrices$beta)>0 && nrow(initRes$matrices$beta) > 0 && all(is.finite(initRes$matrices$beta))){
      startValues[['beta']] <- initRes$matrices$beta
      if (!missing(beta)){
        startValues[['beta']] <- startValues[['beta']] * is.na(beta)
      }
    }
    if (is.null(startValues[['omega_theta']]) &&  !is.null(initRes$matrices$omega_theta) && ncol(initRes$matrices$omega_theta)>0 && nrow(initRes$matrices$omega_theta) > 0 && all(is.finite(initRes$matrices$omega_theta))){
      startValues[['omega_theta']] <- setSym(initRes$matrices$omega_theta)
      if (!missing(omega_theta)){
        startValues[['omega_theta']] <- startValues[['omega_theta']] * is.na(omega_theta)
      }
    }
    if (is.null(startValues[['delta_theta']])  && !is.null(initRes$matrices$delta_theta) && ncol(initRes$matrices$delta_theta)>0 && nrow(initRes$matrices$delta_theta) > 0 && all(is.finite(initRes$matrices$delta_theta))){
      startValues[['delta_theta']] <- setSym(initRes$matrices$delta_theta)
      if (!missing(delta_theta)){
        startValues[['delta_theta']] <- startValues[['delta_theta']] * is.na(delta_theta)
      }
    }
    if (is.null(startValues[['omega_psi']])  &&!is.null(initRes$matrices$omega_psi) && ncol(initRes$matrices$omega_psi)>0 && nrow(initRes$matrices$omega_psi) > 0 && all(is.finite(initRes$matrices$omega_psi))){
      startValues[['omega_psi']] <- setSym(initRes$matrices$omega_psi)
      if (!missing(omega_psi)){
        startValues[['omega_psi']] <- startValues[['omega_psi']] * is.na(omega_psi)
      }
    }
    if (is.null(startValues[['delta_psi']])  && !is.null(initRes$matrices$delta_psi) && ncol(initRes$matrices$delta_psi)>0 && nrow(initRes$matrices$delta_psi) > 0 && all(is.finite(initRes$matrices$delta_psi))){
      startValues[['delta_psi']] <- setSym(initRes$matrices$delta_psi)
      if (!missing(delta_psi)){
        startValues[['delta_psi']] <- startValues[['delta_psi']] * is.na(delta_psi)
      }
    }
    if (is.null(startValues[['psi']])   && !is.null(initRes$matrices$psi) &&  ncol(initRes$matrices$psi)>0 && nrow(initRes$matrices$psi) > 0 && all(is.finite(initRes$matrices$psi))){
      startValues[['psi']] <- setSym(initRes$matrices$psi)
      if (!missing(psi)){
        startValues[['psi']] <- startValues[['psi']] * is.na(psi)
      }
    }
    if (is.null(startValues[['theta']])   && !is.null(initRes$matrices$theta) &&  ncol(initRes$matrices$theta)>0 && nrow(initRes$matrices$theta) > 0 && all(is.finite(initRes$matrices$theta))){
      startValues[['theta']] <- setSym(initRes$matrices$theta)
      if (!missing(theta)){
        startValues[['theta']] <- startValues[['theta']] * is.na(theta)
      }
    }
    if (missing(fitInd)){
      fitInd <- initRes$mxResults$independence
    }
    if (missing(fitSat)){
      fitInd <- initRes$mxResults$saturated
    }
  }
  
  ### Generate model:
  mod <- generatelvnetmodel(
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
    startValues=startValues,
    lasso = lasso,
    lassoMatrix=lassoMatrix,
    scale=scale,
    nLatents=nLatents,
    mimic=mimic,
    fitFunction=fitFunction)
  
  
  #   capture.output(fitMod <- OpenMx::mxRun(mod, silent = TRUE,
  #                   suppressWarnings = TRUE),type="message")
  fitMod <- OpenMx::mxRun(mod, silent = TRUE, suppressWarnings = TRUE)
  
  
  if (missing(fitSat)){
    # Saturated model:
    satMod <- generatelvnetmodel(
      data = data, 
      lambda = diag(Nvar), 
      psi = matrix(NA,Nvar,Nvar), 
      theta = matrix(0, Nvar,Nvar), 
      name = "saturated",
      sampleSize = sampleSize,
      mimic=mimic,
      fitFunction=fitFunction
    )
    
    capture.output(fitSat <- mxRun(satMod, silent = TRUE,
                                   suppressWarnings = TRUE)  ,type="message")
  }
  
  if (missing(fitInd)){
    
    # Construct Psi
    psiInd <- diag(NA, Nvar, Nvar)
    if (!missing(exogenous)){
      if (is.numeric(exogenous)){
        psiInd[exogenous,exogenous] <- NA
      } else {
        inds <- which(colnames(data) %in% exogenous)
        psiInd[exogenous,exogenous] <- NA
      }
    } 
 
    
    # Independence model:
    indMod <- generatelvnetmodel(
      data = data, 
      lambda = diag(Nvar), 
      psi = psiInd, 
      theta = matrix(0, Nvar, Nvar), 
      name = "independence",
      sampleSize = sampleSize,
      mimic=mimic,
      fitFunction=fitFunction
    )
    
    capture.output(fitInd <- mxRun(indMod, silent = TRUE,
                                   suppressWarnings = TRUE),type="message")
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
  
  sigma <- Results$matrices$sigma_positive
  S <- Results$sampleStats$covMat
  
  #   fitMod@matrices$omega_theta@values
  #   fitMod@algebras$penalty
  #   sum(abs(fitMod@matrices$omega_theta@values[upper.tri(fitMod@matrices$omega_theta@values,FALSE)]))
  #   
  # Compute chi-square:
  # if (lasso != 0){
    # Compute DF from non-zero elements
    # Count number of non-zero parameters:
    # Start with means:

  Pars <- countPars(fitMod, ifelse(lasso==0,sqrt(.Machine$double.eps),lassoTol))

  # Number of variables:
  Results$fitMeasures$nvar <- Pars$nVar
  
  # Number of observations:
  Results$fitMeasures$nobs <- Pars$nObs

    Results$fitMeasures$npar <- Pars$nPar
    Results$fitMeasures$df <-  Pars$DF
#   } else {
#     Results$fitMeasures$npar <- summary(fitMod)$estimatedParameters
#     Results$fitMeasures$df <- summary(fitMod)$degreesOfFreedom  
#   }
    
    #  Ncons = samplesize constant. Set to N if mimic = lavaan:
    if (mimic == "lavaan"){
      Ncons <- sampleSize
    } else {
      Ncons <- sampleSize - 1
    }

  Results$fitMeasures$fmin <- (sum(diag(S %*% corpcor::pseudoinverse(sigma)))- log(det(S %*% corpcor::pseudoinverse(sigma))) - Nvar)/2
  Results$fitMeasures$chisq <- 2 * Ncons * Results$fitMeasures$fmin
  Results$fitMeasures$pvalue <- pchisq(Results$fitMeasures$chisq, Results$fitMeasures$df, lower.tail = FALSE)
  
  # Baseline model:
  sigmaBase <- fitInd$algebras$sigma$result
  Results$fitMeasures$baseline.chisq <- Ncons * (sum(diag(S %*% corpcor::pseudoinverse(sigmaBase)))- log(det(S %*% corpcor::pseudoinverse(sigmaBase))) - Nvar)
  Results$fitMeasures$baseline.df <- countPars(fitInd)$DF
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
  Results$fitMeasures$rmsea <- sqrt( max(Tm - dfm,0) / (Ncons * dfm))
  
  # Codes for rmsea confidence interval taken from lavaan:
  lower.lambda <- function(lambda) {
    (pchisq(Tm, df=dfm, ncp=lambda) - 0.95)
  }
  if(is.na(Tm) || is.na(dfm)) {
    Results$fitMeasures$rmsea.ci.lower <- NA
  } else if(dfm < 1 || lower.lambda(0) < 0.0) {
    Results$fitMeasures$rmsea.ci.lower <- 0
  } else {
    if (lower.lambda(0) * lower.lambda(Tm) > 0){
      lambda.l <- NA
    } else {
      lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=Tm)$root,
                      silent=TRUE)      
    }
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
    
    if (upper.lambda(0) * upper.lambda(N.RMSEA) > 0){
      lambda.u <- NA
    } else {
      
      lambda.u <- try(uniroot(f=upper.lambda, lower=0,upper=N.RMSEA)$root,
                      silent=TRUE)  
    }
    
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
  
  # Add extended bic:
  Results$fitMeasures$ebic <-  -2*LL + Results$fitMeasures$npar * log(sampleSize) + 4 *  Results$fitMeasures$npar * ebicTuning * log(sampleSize)  
  
  Results$fitMeasures$ebicTuning <- ebicTuning
  
  class(Results) <- "lvnet"
  
  return(Results)
}
