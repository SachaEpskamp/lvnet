
# Model matrices should contain NA for free elements and a value for fixed elements.
generateRIMmodel <- function(
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
  model = c("rim","sem")){
  
  # Check for input:
  stopifnot(is.matrix(data)|is.data.frame(data))
  
  # Number of variables:
  Nvar <- ncol(data)
  
  stopifnot(model[[1]] %in% c("rim","sem"))
  
  # Check matrices:
  # Lambda (Default to iden if psi is missing or full if psi is not)
  if (missing(lambda)){
    if (model[[1]] == "sem" && !missing(psi)){
      lambda <- matrix(NA, Nvar, ncol(psi))
    } else if (model[[1]] == "rim" && !missing(omega_psi)){
      lambda <- matrix(NA, Nvar, ncol(omega_psi))
    } else {
      lambda <- diag(Nvar)
    }
  }
  
  Nlat <- ncol(lambda)
  
  # psi (default to freely estimable; default scaling in psi)
  if (model[[1]] == "sem"){
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
    
  } else {
    
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
    
  }
  
  
  
  
  # Beta (defaults to null)
  if (missing(beta)){
    beta <- matrix(0, Nlat, Nlat)
  }
  
  if (is.null(colnames(data))){
    colnames(data) <- paste0("y",seq_len(ncol(data)))
  }
  
  # Mx data
  if (ncol(data) == nrow(data) && isSymmetric(data)){
    if (missing(sampleSize)){
      stop("sampleSize needs to be assigned if input is covariance matrix.")
    }
    
    #     Mx_data <- OpenMx::mxData(observed = data, type = "cov", numObs = sampleSize)
    #     # means:
    #     Mx_means <- OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = ncol(data), values=0, 
    #                          free=FALSE, name = "means", dimnames = list("mean",colnames(data))
    #     )
    covMat <- data * (sampleSize - 1)/sampleSize
    
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
  # 
  Mx_data <- OpenMx::mxData(observed = covMat, type = "cov", numObs = sampleSize)
  # means:
  Mx_means <- OpenMx::mxMatrix(type = "Full", nrow = 1, ncol = ncol(data), values=0, 
                               free=FALSE, name = "means", dimnames = list("mean",colnames(data))
  )
  
  #     # Sample covariances:
  #     Mx_S <- OpenMx::mxMatrix(
  #       type = "Symm",
  #       nrow = nrow(covMat),
  #       ncol = ncol(covMat),
  #       free = FALSE,
  #       values = covMat,
  #       name = "S"
  #       )
  #   
  #     # Number of parameters:
  #     if (model[[1]] == "rim"){
  #       nPar <- 
  #         sum(is.na(lambda)) + 
  #           sum(is.na(psi)&upper.tri(psi,diag=TRUE)) + 
  #           sum(is.na(delta)) + 
  #           sum(is.na(omega)&upper.tri(omega,diag=TRUE)) + 
  #           sum(is.na(beta))
  #   
  #     } else {
  #       nPar <- 
  #         sum(is.na(lambda)) + 
  #           sum(is.na(psi)&upper.tri(psi,diag=TRUE)) + 
  #           sum(is.na(theta)&upper.tri(theta,diag=TRUE)) +
  #           sum(is.na(beta))
  #   
  #     }
  #     Mx_p <- OpenMx::mxMatrix(
  #       type = "Full",
  #       nrow = 1,
  #       ncol = 1,
  #       free = FALSE,
  #       values = Nvar,
  #       name = "p"
  #       )
  # 
  #     Mx_n <- OpenMx::mxMatrix(
  #         type = "Full",
  #         nrow=1,
  #         ncol=1,
  #         free=FALSE,
  #         value = sampleSize,
  #         name = "n"
  #       )
  #   
  ### MODEL MATRICES ###
  
  ### RIM AND SEM ###
  
  # Lambda:
  Mx_lambda <- OpenMx::mxMatrix(
    type = "Full",
    nrow = nrow(lambda),
    ncol = ncol(lambda),
    free = is.na(lambda),
    values = ifelse(is.na(lambda),1,lambda),
    name = "lambda"
  )
  
  
  # Beta:
  Mx_beta <- OpenMx::mxMatrix(
    type = "Full",
    nrow = nrow(beta),
    ncol = ncol(beta),
    free = is.na(beta),
    values = 0,
    name = "beta"
  )
  
  
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
  
  # Expectation:
  expFunction <- OpenMx::mxExpectationNormal(covariance = "sigma")
  
  
  # Fit function:
  fitFunction <- OpenMx::mxFitFunctionML()
  
  ### SEM MODEL ####
  if (model[[1]] == "sem"){
    # Psi:
    Mx_psi <- OpenMx::mxMatrix(
      type = "Symm",
      nrow = nrow(psi),
      ncol = ncol(psi),
      free = is.na(psi),
      values = ifelse(is.na(psi),diag(ncol(psi)),psi),
      lbound = ifelse(diag(nrow(psi)) == 1, 0, NA),
      name = "psi"
    )
    
    
    # Theta:
    Mx_theta <- OpenMx::mxMatrix(
      type = "Symm",
      nrow = nrow(theta),
      ncol = ncol(theta),
      free = is.na(theta),
      values = diag(nrow(theta)),
      name = "theta"
    )
    
    # SEM model:
    Mx_sigma <- OpenMx::mxAlgebra(
      lambda %*% solve(I_lat - beta) %*% psi %*% t(solve(I_lat - beta)) %*% t(lambda) + theta, 
      name = "sigma",
      dimnames = list(colnames(data),colnames(data)))
    
    
    Mx_model <- OpenMx::mxModel(
      name = name,
      Mx_data,
      Mx_means,
      Mx_lambda,
      Mx_psi,
      Mx_identity_lat,
      Mx_theta,
      Mx_beta,
      Mx_sigma,
      expFunction,
      fitFunction
    )
    
  } else {
    ### RIM MODEL
    # Delta:
    Mx_delta_theta <- OpenMx::mxMatrix(
      type = "Diag",
      nrow = nrow(delta_theta),
      ncol = ncol(delta_theta),
      free = is.na(delta_theta),
      values = ifelse(is.na(delta_theta),1,delta_theta),
      lbound = 0,
      name = "delta_theta"
    )
    
    # Omega:
    Mx_omega_theta <- OpenMx::mxMatrix(
      type = "Symm",
      nrow = nrow(omega_theta),
      ncol = ncol(omega_theta),
      free = is.na(omega_theta),
      values = 0,
      lbound = ifelse(diag(nrow(omega_theta)) == 1,0, -1),
      ubound = ifelse(diag(nrow(omega_theta)) == 1,0, 1),
      name = "omega_theta"
    )
    
    
    Mx_delta_psi <- OpenMx::mxMatrix(
      type = "Diag",
      nrow = nrow(delta_psi),
      ncol = ncol(delta_psi),
      free = is.na(delta_psi),
      values = ifelse(is.na(delta_psi),1,delta_psi),
      lbound = 0,
      name = "delta_psi"
    )
    
    # Omega:
    Mx_omega_psi <- OpenMx::mxMatrix(
      type = "Symm",
      nrow = nrow(omega_psi),
      ncol = ncol(omega_psi),
      free = is.na(omega_psi),
      values = 0,
      lbound = ifelse(diag(nrow(omega_psi)) == 1,0, -1),
      ubound = ifelse(diag(nrow(omega_psi)) == 1,0, 1),
      name = "omega_psi"
    )
    
    
    # RIM model:
    Mx_sigma <- OpenMx::mxAlgebra(
      lambda %*% solve(I_lat - beta) %*% delta_psi %*% solve(I_lat - omega_psi) %*% delta_psi %*% t(solve(I_lat - beta)) %*% t(lambda) + delta_theta %*% solve(I_obs - omega_theta) %*% delta_theta, 
      name = "sigma",
      dimnames = list(colnames(data),colnames(data)))
    
    Mx_model <- OpenMx::mxModel(
      name = name,
      Mx_data,
      Mx_means,
      Mx_lambda,
      Mx_delta_theta,
      Mx_omega_theta,
      Mx_delta_psi,
      Mx_omega_psi,
      Mx_identity_obs,
      Mx_identity_lat,
      Mx_beta,
      Mx_sigma,
      expFunction,
      fitFunction
    )
    
  }
  
 
  
  return(Mx_model)
}