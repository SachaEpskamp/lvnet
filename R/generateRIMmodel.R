
# Model matrices should contain NA for free elements and a value for fixed elements.
generateRIMmodel <- function(
  data, # Raw data or a covariance matrix
  lambda, # Lambda design matrix. NA indicates free parameters. If missing and psi is missing, defaults to identity matrix with warning
  omega, # Network. If missing, defaults to matrix of zeroes
  psi, # Latent variance-covariance matrix. If missing, defaults to free
  beta, # Structural matrix. If missing, defaults to zero.
  delta, # Scaling matrix, can be missing
  theta, # Used if model = "sem". Defaults to diagonal
  sampleSize,
  name = "mod",
  model = c("rim","sem")){
  
  # Check for input:
  stopifnot(is.matrix(data)|is.data.frame(data))
  
  # Number of variables:
  Nvar <- ncol(data)
  
  # Check matrices:
  # Lambda (Default to iden if psi is missing or full if psi is not)
  if (missing(lambda)){
    if (missing(psi)){
      lambda <- diag(Nvar)
    } else {
      lambda <- matrix(NA, Nvar, ncol(psi))
    }
  }
  
  Nlat <- ncol(lambda)
  
  # psi (default to freely estimable; default scaling in psi)
  if (missing(psi)){
    if (missing(beta)){
      psi <- matrix(NA, Nlat, Nlat)
      diag(psi) <- 1
    } else {
      psi <- diag(1, Nlat)
    }
  }

  
  # omega (default to null)
  if (missing(omega)){
    omega <- matrix(0, Nvar, Nvar)
  }
  if (any(is.na(diag(omega))) || any(diag(omega)!=0)){
    warning("'omega' must have zero diagonal. Set to zero.")
    diag(omega) <- 0
  }
  
  # Delta (defaults to diagonal)
  if (missing(delta)){
    delta <- diag(NA, Nvar)
  }
  
  # Theta (defaults to diagonal)
  if (missing(theta)){
    theta <- diag(NA, Nvar)
  }
  
  # Beta (defaults to null)
  if (missing(beta)){
    beta <- matrix(0, Nlat, Nlat)
  }

  stopifnot(isSymmetric(psi))
  stopifnot(isSymmetric(omega))
  
  
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

  # Lambda:
  Mx_lambda <- OpenMx::mxMatrix(
    type = "Full",
    nrow = nrow(lambda),
    ncol = ncol(lambda),
    free = is.na(lambda),
    values = ifelse(is.na(lambda),1,lambda),
    name = "lambda"
  )
  
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
  
  # Delta:
  Mx_delta <- OpenMx::mxMatrix(
    type = "Diag",
    nrow = nrow(delta),
    ncol = ncol(delta),
    free = is.na(delta),
    values = ifelse(is.na(delta),1,delta),
    lbound = 0,
    name = "delta"
  )
  
  # Omega:
  Mx_omega <- OpenMx::mxMatrix(
    type = "Symm",
    nrow = nrow(omega),
    ncol = ncol(omega),
    free = is.na(omega),
    values = 0,
    lbound = ifelse(diag(nrow(omega)) == 1,0, -1),
    ubound = ifelse(diag(nrow(omega)) == 1,0, 1),
    name = "omega"
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
  
  # Theta:
  Mx_theta <- OpenMx::mxMatrix(
    type = "Symm",
    nrow = nrow(theta),
    ncol = ncol(theta),
    free = is.na(theta),
    values = diag(nrow(theta)),
    name = "theta"
  )
  
  Mx_identity_omega <- OpenMx::mxMatrix(
    type = "Iden",
    nrow = nrow(omega),
    ncol = ncol(omega),
    name = "I_omega"
  )
  
  
  Mx_identity_beta <- OpenMx::mxMatrix(
    type = "Iden",
    nrow = nrow(beta),
    ncol = ncol(beta),
    name = "I_beta"
  )
  
  # Sigma:
  if (model[[1]]=="rim"){
    # RIM model:
    Mx_sigma <- OpenMx::mxAlgebra(lambda %*% solve(I_beta - beta) %*% psi %*% t(solve(I_beta - beta)) %*% t(lambda) + delta %*% solve(I_omega - omega) %*% delta, 
                          name = "sigma",
                          dimnames = list(colnames(data),colnames(data)))
  } else if (model[[1]] == "sem"){
    # SEM model:
    Mx_sigma <- OpenMx::mxAlgebra(lambda %*% solve(I_beta - beta) %*% psi %*% t(solve(I_beta - beta)) %*% t(lambda) + theta, 
                          name = "sigma",
                          dimnames = list(colnames(data),colnames(data)))
    
    
  } else stop(paste("Model",model,"unknown."))
  
  # Expectation:
#   expFunction <- OpenMx::mxExpectationNormal(covariance = "sigma", means = "means")
expFunction <- OpenMx::mxExpectationNormal(covariance = "sigma")


  # Fit function:
  fitFunction <- OpenMx::mxFitFunctionML()
#     Mx_lik <- mxAlgebra(log(det(sigma)) + tr(S %*% solve(sigma)) - log(det(S)) - p,
#                          name = "loglik")
# 
# 
#     Mx_fit <- mxFitFunctionAlgebra("loglik", numObs = sampleSize, numStats = Nvar*(Nvar+1)/2)
#                                
  # Define model:
  if (model[[1]]=="rim"){
    Mx_model <- OpenMx::mxModel(
      name = name,
      Mx_data,
      Mx_means,
      Mx_lambda,
      Mx_psi,
      Mx_delta,
      Mx_omega,
      Mx_identity_beta,
      Mx_identity_omega,
      Mx_beta,
      Mx_sigma,
      expFunction,
      fitFunction
#       Mx_S,
#       Mx_p,
#       Mx_n,
#       Mx_lik,
#       Mx_fit
      )
  } else {
    Mx_model <- OpenMx::mxModel(
      name = name,
      Mx_data,
      Mx_means,
      Mx_lambda,
      Mx_psi,
      Mx_identity_beta,
      Mx_theta,
      Mx_beta,
      Mx_sigma,
      expFunction,
      fitFunction
#       Mx_S,
#       Mx_p,
#       Mx_n,
#       Mx_lik,
#       Mx_fit
      )
  }

  
  return(Mx_model)
}