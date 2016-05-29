start <- function(mat,list,alt){
  if (!is.null(list[[mat]])){
    return(list[[mat]])
  } else {
    return(alt)
  }
}

# Model matrices should contain NA for free elements and a value for fixed elements.
generatelvnetmodel <- function(
  data, # Raw data or a covariance matrix
  lambda, # Lambda design matrix. NA indicates free parameters. If missing and psi is missing, defaults to identity matrix with warning
  omega_theta, # Observed residual network. If missing, defaults to matrix of zeroes
  delta_theta, # Scaling matrix, can be missing
  omega_psi, # Latent residual network. If missing, defaults to matrix of zeroes
  delta_psi, # Scaling matrix, can be missing
  beta, # Structural matrix. If missing, defaults to zero.
  psi, # Latent variance-covariance matrix. If missing, defaults to free
  theta, # Used if model = "sem". Defaults to diagonal
  sampleSize,
  name = "mod",
  startValues = list(),
  lasso = 0,
  lassoMatrix = "",
  scale = TRUE,
  nLatents # allows for quick specification of fully populated lambda matrix.
  ){
  
  # Silly things to fool R check:
  I_lat <- NULL
  I_obs <- NULL
  PsiPlus <- NULL
  vec2diag <- NULL
  diag2vec <- NULL
  theta_inverse <- NULL
  eigenval <- NULL
  sigma_positive <- NULL
  P <- NULL
  penalty <- NULL
  
  # Check for input:
  stopifnot(is.matrix(data)|is.data.frame(data))
  
  # Number of variables:
  Nvar <- ncol(data)
  
  #   stopifnot(model[[1]] %in% c("lvnet","sem"))
  
  # Check matrices:
  # Lambda (Default to iden if psi is missing or full if psi is not)
  if (missing(lambda)){
    if (!missing(nLatents)){
      lambda <- matrix(NA,Nvar,nLatents)
    } else if (!missing(psi)){
      lambda <- matrix(NA, Nvar, ncol(psi))
    } else if (!missing(omega_psi)){
      lambda <- matrix(NA, Nvar, ncol(omega_psi))
    } else {
      lambda <- matrix(,Nvar,0)
      
      if (missing(theta) && missing(omega_theta)){
        if (lasso !=0){
          if (!(any(c("omega_theta","theta") %in% lassoMatrix))){
            theta <- matrix(NA, Nvar, Nvar)
          }
        } else {
          theta <- matrix(NA, Nvar, Nvar)
        }
      }
    }
  }
  
  Nlat <- ncol(lambda)
  
  if (nrow(lambda) != ncol(data)){
    stop("Number of rows in 'lambda' does not equal number of variables in 'data'")
  }
  
  # Check lasso matrix and set missing if needed:
  if (!missing(lassoMatrix)){
    if ("omega_psi" %in% lassoMatrix && missing(omega_psi)){
      omega_psi <- matrix(NA,Nlat,Nlat)
      diag(omega_psi) <- 0
    }
    if ("psi" %in% lassoMatrix && missing(psi)){
      psi <- matrix(NA,Nlat,Nlat)
      diag(psi) <- 1
    }
    if ("omega_theta" %in% lassoMatrix && missing(omega_theta)){
      omega_theta <- matrix(NA,Nvar,Nvar)
      diag(omega_theta) <- 0
    }
    if ("theta" %in% lassoMatrix && missing(theta)){
      theta <- matrix(NA,Nvar,Nvar)
    }
  }

  # psi and omega_psi may not both be missing:
  if(!missing(psi) & !missing(omega_psi)){
    stop("Both 'psi' and 'omega_psi' modeled.")
  }

  # Both theta and omega_theta may not be missing:
  if(!missing(theta) & !missing(omega_theta)){
    stop("Both 'theta' and 'omega_theta' modeled.")
  }
  
  # Identify if estimation should be on psi or omega_psi:
  estPsi <- (!missing(psi)) | (missing(psi) & missing(omega_psi))
  
  # Identify if estimation should be on theta or omega_theta:
  estTheta <- (!missing(theta)) | (missing(theta) & missing(omega_theta))
  
  # psi (default to freely estimable; default scaling in psi)
  if (missing(psi)){
    if (missing(beta)){
      psi <- matrix(NA, Nlat, Nlat)
      diag(psi) <- 1
    } else {
      psi <- diag(1, Nlat)
    }
  } 
  
  
  # Theta (defaults to diagonal)
  if (missing(theta)){
    theta <- diag(NA, Nvar)
  }
  
  stopifnot(isSymmetric(psi))
  stopifnot(isSymmetric(theta))    
  
  
  # Omega psi:
  if (missing(omega_psi)){
    if (missing(beta)){
      omega_psi <- matrix(NA, Nlat, Nlat)
      diag(omega_psi) <- 0
    } else {
      omega_psi <- matrix(0, Nlat, Nlat)
    }
  }
  
  # omega_theta (default to null)
  if (missing(omega_theta)){
    omega_theta <- matrix(0, Nvar, Nvar)
  }
  
  # Check omegas:
  if (any(is.na(diag(omega_theta))) || any(diag(omega_theta)!=0)){
    warning("'omega_theta' must have zero diagonal. Set to zero.")
    diag(omega_theta) <- 0
  }
  
  if (any(is.na(diag(omega_psi))) || any(diag(omega_psi)!=0)){
    warning("'omega_psi' must have zero diagonal. Set to zero.")
    diag(omega_psi) <- 0
  }
  
  # Delta (defaults to diagonal)
  if (missing(delta_theta)){
    delta_theta <- diag(NA, Nvar)
  }
  if (missing(delta_psi)){
    delta_psi <- diag(1, Nlat)
  }
  
  stopifnot(isSymmetric(omega_psi))
  stopifnot(isSymmetric(omega_theta))
  
  # Beta (defaults to null)
  if (missing(beta)){
    beta <- matrix(0, Nlat, Nlat)
  }
  
  if (is.null(colnames(data))){
    colnames(data) <- paste0("y",seq_len(ncol(data)))
  }

  # Mx data
  if (ncol(data) == nrow(data) && isSymmetric(unname(data))){
    if (missing(sampleSize)){
      stop("sampleSize needs to be assigned if input is covariance matrix.")
    }
    
    #     Mx_data <- OpenMx::mxData(observed = data, type = "cov", numObs = sampleSize)
    #     # means:
    #     Mx_means <- OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = ncol(data), values=0, 
    #                          free=FALSE, name = "means", dimnames = list("mean",colnames(data))
    #     )

    covMat <- data * (sampleSize - 1)/sampleSize
    rownames(covMat) <- colnames(covMat)
    
  } else {
    sampleSize <- nrow(data)   
    #     Mx_data <- OpenMx::mxData(observed = data, type = "raw")
    #     # means:
    #     Mx_means <- OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = ncol(data), values=colMeans(data), 
    #                          free=FALSE, name = "means", dimnames = list("mean",colnames(data))
    #     )
    
    data <- as.matrix(data)
    covMat <- cov(data, use = "pairwise.complete.obs")* (sampleSize - 1)/sampleSize
  }
  
  if (scale){
    covMat <- setSym(cov2cor(covMat))
  } else {
    if (lasso != 0){
      warning("It is advised to set 'scale = TRUE' when using LASSO estimation.")
    }
  }
  
  # 
  Mx_data <- OpenMx::mxData(observed = covMat, type = "cov", numObs = sampleSize)
  # means:
  Mx_means <- OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = ncol(data), values=0, 
                               free=FALSE, name = "means", dimnames = list("mean",colnames(data))
  )
  
  ### lvnet AND SEM ###
  # Lambda:
  if (Nlat > 0){

    Mx_lambda <- OpenMx::mxMatrix(
      type = "Full",
      nrow = nrow(lambda),
      ncol = ncol(lambda),
      free = is.na(lambda),
      values = start("lambda",startValues,ifelse(is.na(lambda),1,lambda)),
      name = "lambda"
    )
  } else {
    Mx_lambda <- OpenMx::mxMatrix(
      type = "Full",
      nrow = nrow(lambda),
      ncol = ncol(lambda),
      name = "lambda"
    )  
  } 
  
  # Beta:
  if (Nlat > 0){
    Mx_beta <- OpenMx::mxMatrix(
      type = "Full",
      nrow = nrow(beta),
      ncol = ncol(beta),
      free = is.na(beta),
      values = start("beta",startValues,0),
      name = "beta"
    )
  } else {
    Mx_beta <- OpenMx::mxMatrix(
      type = "Full",
      nrow = nrow(beta),
      ncol = ncol(beta),
      name = "beta"
    )
  }
  
  
  Mx_identity_lat <- OpenMx::mxMatrix(
    type = "Iden",
    nrow = Nlat,
    ncol = Nlat,
    name = "I_lat"
  )
  
  
  Mx_identity_obs <- OpenMx::mxMatrix(
    type = "Iden",
    nrow = Nvar,
    ncol = Nvar,
    name = "I_obs"
  )
  
  
  
  # Psi and omega_psi:
  if (estPsi){
    if (Nlat > 0){
      
      # bounds:
      ubound <- sqrt(diag(psi)) %o% sqrt(diag(psi))
      
      # Force no cors of 1:
      ubound <- 0.99*ubound
      diag(ubound) <- diag(psi)
      
      ubound <- ifelse(is.na(ubound),Inf,ubound)
      lbound <- -ubound
      diag(lbound) <- 0
      

      Mx_psi <- OpenMx::mxMatrix(
        type = "Symm",
        nrow = nrow(psi),
        ncol = ncol(psi),
        free = is.na(psi),
        values = start("psi",startValues,ifelse(is.na(psi),diag(ncol(psi)),psi)),
        lbound = lbound,
        ubound = ubound,
        name = "psi"
      )
    } else {
      Mx_psi <- OpenMx::mxMatrix(
        type = "Symm",
        nrow = nrow(psi),
        ncol = ncol(psi),
        name = "psi"
      )
    }

    Mx_delta_psi <- OpenMx::mxAlgebra(
      vec2diag(1/sqrt(diag2vec(solve(psi)))),
      name = "delta_psi"
    )
    
    Mx_omega_psi <- OpenMx::mxAlgebra(
      I_lat - delta_psi %*% solve(psi) %*% delta_psi,
      name = "omega_psi"
    )
    
#     Mx_psi_inverse <- OpenMx::mxAlgebra(
#       solve(psi),
#       name = "psi_inverse"
#     )
    
  } else {
    
    Mx_delta_psi <- OpenMx::mxMatrix(
      type = "Diag",
      nrow = nrow(delta_psi),
      ncol = ncol(delta_psi),
      free = is.na(delta_psi),
      values = start("delta_psi",startValues,ifelse(is.na(delta_psi),1,delta_psi)),
      lbound = 0,
      name = "delta_psi"
    )
#     
#     # Omega:
    if (Nlat > 0){
      Mx_omega_psi <- OpenMx::mxMatrix(
        type = "Symm",
        nrow = nrow(omega_psi),
        ncol = ncol(omega_psi),
        free = is.na(omega_psi),
        values = start("omega_psi",startValues,0),
        lbound = ifelse(diag(nrow(omega_psi)) == 1,0, -0.99),
        ubound = ifelse(diag(nrow(omega_psi)) == 1,0, 0.99),
        name = "omega_psi",
      )      
    } else {
      Mx_omega_psi <- OpenMx::mxMatrix(
        type = "Symm",
        nrow = nrow(omega_psi),
        ncol = ncol(omega_psi),
        name = "omega_psi"
      )      
    }
#         if (Nlat > 0){
#           Mx_psi_inverse <- OpenMx::mxMatrix(
#             type = "Symm",
#             nrow = Nlat,
#             ncol = Nlat,
#             free = is.na(omega_psi) | is.na(delta_psi),
#             values = diag(1, Nlat, Nlat),
#             lbound = ifelse(diag(nrow(omega_psi)) == 1,0, -1),
#             ubound = 1,
#             name = "psi_inverse"
#           )      
#         } else {
#           Mx_psi_inverse <- OpenMx::mxMatrix(
#             type = "Symm",
#             nrow = nrow(omega_psi),
#             ncol = ncol(omega_psi),
#             name = "psi_inverse"
#           )
#         }
    
#     Mx_omega_psi <-  OpenMx::mxAlgebra(
#       I_obs - delta_psi %*% psi_inverse %*% delta_psi, 
#       name = "omega_psi"
#     )
#     
    Mx_psi <- OpenMx::mxAlgebra(
      delta_psi %*% solve(I_lat - omega_psi) %*% delta_psi,
      name = "psi"
    )
    
#     Mx_delta_psi <- OpenMx::mxAlgebra(
#       vec2diag(1/sqrt(diag2vec(psi_inverse))),
#       name = "delta_psi"
#     )
#     
#     Mx_omega_theta <- OpenMx::mxAlgebra(
#       I_obs - delta_theta %*% theta_inverse %*% delta_theta,
#       name = "omega_theta"
#     )
#     
#     Mx_theta <- OpenMx::mxAlgebra(
#       solve(theta_inverse),
#       name = "theta"
#     )

    

  }
  
  # Theta and omega_theta:
  if (estTheta){
    
    Mx_theta <- OpenMx::mxMatrix(
      type = "Symm",
      nrow = nrow(theta),
      ncol = ncol(theta),
      free = is.na(theta),
      values = start("theta",startValues,ifelse(is.na(theta),diag(nrow(theta)),theta)),
      name = "theta"
    )
    
    Mx_delta_theta <- OpenMx::mxAlgebra(
      vec2diag(1/sqrt(diag2vec(solve(theta)))),
      name = "delta_theta"
    )
    
    Mx_omega_theta <- OpenMx::mxAlgebra(
      I_obs - delta_theta %*% solve(theta) %*% delta_theta,
      name = "omega_theta"
    )
    
    
    Mx_theta_inverse <- OpenMx::mxAlgebra(
      solve(theta),
      name = "theta_inverse"
    )
    
  } else {

    if (is.null(startValues[["delta_theta"]])){
      startValues[["delta_theta"]] <- diag(1/sqrt(diag(corpcor::pseudoinverse(covMat))))
    }
    
    
    ##### TEST CODES ####
    
    # Only estimate inverse of theta, obtain omega_theta and delta_theta afterwards:
    # Use

    Mx_theta_inverse <- OpenMx::mxMatrix(
      type = "Symm",
      nrow = Nvar,
      ncol = Nvar,
      free = is.na(omega_theta) | is.na(delta_theta),
      values = diag(1,Nvar),
      lbound = ifelse(diag(Nvar) == 1,0, -Inf),
      ubound = Inf,
      name = "theta_inverse"
    )
    
    Mx_delta_theta <- OpenMx::mxAlgebra(
      vec2diag(1/sqrt(diag2vec(theta_inverse))),
      name = "delta_theta"
    )
    
    Mx_omega_theta <- OpenMx::mxAlgebra(
      I_obs - delta_theta %*% theta_inverse %*% delta_theta,
      name = "omega_theta"
    )
    
    Mx_theta <- OpenMx::mxAlgebra(
      solve(theta_inverse),
      name = "theta"
    )
    
    #
    
    
#     
#     # Delta:
#     Mx_delta_theta <- OpenMx::mxMatrix(
#       type = "Diag",
#       nrow = nrow(delta_theta),
#       ncol = ncol(delta_theta),
#       free = is.na(delta_theta),
#       values = start("delta_theta",startValues,ifelse(is.na(delta_theta),1,delta_theta)),
#       lbound = 0,
#       name = "delta_theta"
#     )
# 
#     # Omega:
#     Mx_omega_theta <- OpenMx::mxMatrix(
#       type = "Symm",
#       nrow = nrow(omega_theta),
#       ncol = ncol(omega_theta),
#       free = is.na(omega_theta),
#       values = start("omega_theta",startValues,0),
#       lbound = ifelse(diag(nrow(omega_theta)) == 1,0, -0.99),
#       ubound = ifelse(diag(nrow(omega_theta)) == 1,0, 0.99),
#       name = "omega_theta"
#     )
#     
#     Mx_theta <- OpenMx::mxAlgebra(
#       delta_theta %*% solve(I_obs - omega_theta) %*% delta_theta,
#       name = "theta"
#     )
  }

  
  # Fake psi with shifted eigenvalues if needed:
#   Mx_Psi_Positive <- OpenMx::mxAlgebra(
#     psi - min(0,(min(eigenval(psi))-.00001)) * I_lat, name = "psi_positive"
#   )
#   Mx_Psi_Positive <- OpenMx::mxAlgebra(
#     psi -0* I_lat, name = "psi_positive"
#   )
  # Constraint on psi:
  # Small values for diagonal:
#   Mx_PsiDiagplus <- mxMatrix(
#     "Diag",
#     nrow(psi),
#     ncol(psi),
#     FALSE,
#     values = 1e-5,
#     name = "PsiPlus"
#   )
#   Mx_PsiCon <- OpenMx::mxConstraint(psi < sqrt(diag2vec(psi)) %*% sqrt(t(diag2vec(psi))) + PsiPlus)
#   
  

  
  ### FIT FUNCTIONS ###
#   if (lasso ==0){
#     # Expectation:
#     expFunction <- OpenMx::mxExpectationNormal(covariance = "sigma")
#     
#     # Fit function:
#     fitFunction <- OpenMx::mxFitFunctionML()
#     
#     # Implied covariance:
#     if (Nlat > 0){
#       Mx_sigma <- OpenMx::mxAlgebra(
#         lambda %*% solve(I_lat - beta) %*% psi %*% t(solve(I_lat - beta)) %*% t(lambda) + theta, 
#         name = "sigma",
#         dimnames = list(colnames(data),colnames(data)))    
#       
#       ### Model:
#       Mx_model <- OpenMx::mxModel(
#         name = name,
#         Mx_data,
#         Mx_means,
#         Mx_lambda,
#         Mx_theta,
#         Mx_psi,
#         Mx_delta_theta,
#         Mx_omega_theta,
#         Mx_delta_psi,
#         Mx_omega_psi,
#         Mx_identity_obs,
#         Mx_identity_lat,
#         Mx_beta,
#         Mx_sigma,
#         expFunction,
#         fitFunction
#         
#         # Mx_PsiCon,
#         # Mx_PsiDiagplus
#       )
#     } else {
#       Mx_sigma <- OpenMx::mxAlgebra(
#         theta, 
#         name = "sigma",
#         dimnames = list(colnames(data),colnames(data)))
#       
#       
#       ### Model:
#       Mx_model <- OpenMx::mxModel(
#         name = name,
#         Mx_data,
#         Mx_means,
#         Mx_theta,
#         Mx_delta_theta,
#         Mx_omega_theta,
#         Mx_identity_obs,
#         Mx_sigma,
#         expFunction,
#         fitFunction
#       )
#     }
#     
#   } else {
    # Observed covariance matrix (used in mxAlgebra for LASSO):
    mx_observedCovs <- OpenMx::mxMatrix(
      "Symm",
      nrow = nrow(covMat),
      ncol = ncol(covMat),
      free = FALSE,
      values = covMat,
      dimnames = dimnames(covMat),
      name = "C"
    )
    
    # Positive definite shifted sigma:
    Mx_Sigma_positive <- OpenMx::mxAlgebra(
      sigma - min(0,(min(eigenval(sigma))-.00001)) * I_obs, name = "sigma_positive"
    )
    
    # Tuning parameter:
    mx_Tuning <- OpenMx::mxMatrix(nrow=1,ncol=1,free=FALSE,values=lasso,name="tuning")
      
    # Number of obsrved variables:
    mx_P <- OpenMx::mxMatrix(nrow=1,ncol=1,free=FALSE,values=nrow(covMat),name="P")
    
    # LASSO penalty:
    # Construct the penalty:

    if (!missing(lassoMatrix) && length(lassoMatrix) > 0){
      
      # Put vechs around any matrix that is not lambda or beta:
      penString <- ifelse(lassoMatrix %in% c("lambda","beta"), lassoMatrix,
                          paste0("vech(",lassoMatrix,")")
                          )
      
      # sum absolute values and plus::
      penString <- paste0("sum(abs(",penString,"))", collapse = " + ")
      
      # Multiply with tuning:
      penString <- paste0("tuning * (",penString,")")
      
      # Add to penalty:
      Penalty <- OpenMx::mxAlgebraFromString( penString,name = "penalty")
 
    } else {
      Penalty <- OpenMx::mxAlgebra(0,name = "penalty")
    }

    
    
    # LASSO fit function:
    # logLik <- OpenMx::mxAlgebra(log(det(sigma)) + tr(C %*% solve(sigma)) - log(det(C)) - P + penalty,name = "logLik")
    
    # logLik <- OpenMx::mxAlgebra(log(det(sigma)) + tr(C %*% solve(sigma)) - log(det(C)) - P,name = "logLik")
    logLik <- OpenMx::mxAlgebra(log(det(sigma_positive)) + tr(C %*% solve(sigma_positive)) - log(det(C)) - P + penalty,name = "logLik")
    
    
    # Fit function:
    fitFunction <- OpenMx::mxFitFunctionAlgebra("logLik", numObs = sampleSize, numStats = ncol(covMat)*(ncol(covMat)+1)/2)
  
    if (Nlat > 0){
      Mx_sigma <- OpenMx::mxAlgebra(
        lambda %*% solve(I_lat - beta) %*% psi %*% t(solve(I_lat - beta)) %*% t(lambda) + theta, 
        name = "sigma",
        dimnames = list(colnames(data),colnames(data)))    
      
      ### Model:
      Mx_model <- OpenMx::mxModel(
        name = name,
        Mx_data,
        Mx_means,
        Mx_lambda,
        Mx_theta,
        Mx_psi,
        Mx_delta_theta,
        Mx_omega_theta,
        Mx_delta_psi,
        Mx_omega_psi,
        Mx_identity_obs,
        Mx_identity_lat,
        Mx_beta,
        Mx_sigma,
        Mx_theta_inverse,
        # Mx_psi_inverse,

        # LASSO stuff:
        fitFunction,
        mx_Tuning,
        mx_P,
        logLik,
        mx_observedCovs,
        Penalty,
        Mx_Sigma_positive 
      )
      
    } else {
      Mx_sigma <- OpenMx::mxAlgebra(
        theta, 
        name = "sigma",
        dimnames = list(colnames(data),colnames(data)))
      
      
      ### Model:
      Mx_model <- OpenMx::mxModel(
        name = name,
        Mx_data,
        Mx_means,
        Mx_theta,
        Mx_delta_theta,
        Mx_omega_theta,
        Mx_identity_obs,
        Mx_sigma,
        Mx_theta_inverse,
        # Mx_psi_inverse,

        # LASSO Stuff:
        fitFunction,
        mx_Tuning,
        mx_P,
        logLik,
        mx_observedCovs,
        Penalty,
        Mx_Sigma_positive 
      )

    }
    
  # }
  


  
  
  
  
  return(Mx_model)
}
