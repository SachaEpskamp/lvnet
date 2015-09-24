wi2net <- function(x)
{
  x <- -cov2cor(x)
  diag(x) <- 0
  x <- forceSymmetric(x)
  return(x)
}

forcePositive <- function(x){
  if (any(eigen(x)$values<0)){
    cov2cor(x - (min(eigen(x)$values)-.1) * diag(nrow(x)))
  } else {
    x
  }
}
# E-step in the optimization algorithm:
Estep <- function(
  S, # Sample covariance
  Kcur, # Current estimate for K
  obs # Logical indicating observed
)
{
  stopifnot(Matrix::isSymmetric(S))
  if (missing(Kcur))
  {
    if (missing(obs)) 
    {
      Kcur <- diag(nrow(S)) 
    } else 
    {
      Kcur <- diag(length(obs))
    }
  }
  stopifnot(Matrix::isSymmetric(Kcur))
  if (missing(obs)) obs <- 1:nrow(Kcur) %in% 1:nrow(S)
  
  # To make life easier:
  H <- !obs
  O <- obs
  
  
  # Current estimate of S:
  Scur <- corpcor::pseudoinverse(Kcur)
  
  # Expected Sigma_OH:
  Sigma_OH <- S %*% corpcor::pseudoinverse(Scur[O,O]) %*% Scur[O, H]
  #   Sigma_OH
  
  # Expected Sigma_H:
  Sigma_H <- Scur[H, H] - Scur[H,O] %*% corpcor::pseudoinverse(Scur[O,O]) %*% Scur[O,H] + Scur[H,O] %*% corpcor::pseudoinverse(Scur[O,O]) %*% S %*% corpcor::pseudoinverse(Scur[O,O]) %*% Scur[O, H]
  
  # Construct expected sigma:
  Sigma_Exp <- rbind(cbind(S,Sigma_OH),cbind(t(Sigma_OH), Sigma_H))  
  
  return(Sigma_Exp)
}

# M-step in the optimization algorithm:
Mstep <- function(
  Sexp, # Expected full S
  obs, # Logica indiating oberved variables
  rho = 0,
  lambda
)
{
  if (!is.positive.definite(Sexp))
  {
    Sexp <- as.matrix(nearPD(Sexp)$mat)
    warning("Expected covariance matrix is not positive definite")
  }
  # Rho matrix:
  n <- nrow(Sexp)
  RhoMat <- matrix(rho, n, n)
  RhoMat[!obs,] <- 0
  RhoMat[,!obs] <- 0
  
  # Lambda:
  zeroes <- which(!lambda, arr.ind = TRUE)
  zeroes[,2] <- which(!obs)[zeroes[,2]]
  
  if (nrow(zeroes) > 0){
    K <- glasso(Sexp, RhoMat, penalize.diagonal=FALSE, zero = zeroes)$wi 
  } else {
    K <- glasso(Sexp, RhoMat, penalize.diagonal=FALSE)$wi
  }
  
  return(K)  
}

### Main lvglasso function
lvglasso <- function(
  S, # Sample cov
  nLatents, # Number of latents
  rho = 0, # Penalty
  thr = 1.0e-4, # Threshold for convergence (sum absolute diff)
  maxit = 1e4, # Maximum number of iterations
  lambda # Logical matrix indicating the free factor loadings. Defaults to full TRUE matrix.
)
{
  if (missing(nLatents)){
    stop("'nLatents' must be specified")
  }

  nobs <- nrow(S)
  ntot <- nobs + nLatents
  
  if (missing(lambda) || is.null(lambda)){
    lambda <- matrix(TRUE, nobs, nLatents)
  }
  
  if (nrow(lambda) != nobs | ncol(lambda) != nLatents) stop("Dimensions of 'lambda' are wrong.")
  
  # PCA prior for K:
  efaRes <- principal(S, nfactors = nLatents)
  
#   # If sampleSizeis missing, set to 1000. Is only used for prior anyway.
# #   if (missing(sampleSize)){
#     sampleSize <- 1000
# #   }
# 
#   #   Get prior for K:
#   efaRes <- fa(S, n.obs= sampleSize, nfactors=nLatents)
  resid <- residuals(efaRes)
  class(resid) <- "matrix"
  load <- loadings(efaRes)
  class(load) <- "matrix"
  r <- efaRes$r.scores
  class(r) <- "matrix"
  r <- diag(diag(r))
#   
#   # Stupid nonanalytic way to get prior:
#   # Simulate N random variales:
#   #   eta <- rmvnorm(10000, rep(0, nLatents), r)
#   # 
#   #   # Simulate oserved scores:
#   #   Y <- eta %*% t(load) + rmvnorm(10000, rep(0, nobs), diag(diag(resid)
  # RAM FRAMEWORK:
  
  Sym <- rbind(cbind(diag( efaRes$uniquenesses),matrix(0,nobs,nLatents) ),cbind(t(matrix(0,nobs,nLatents) ),r))
  As <- matrix(0, ntot, ntot)
  if (nLatents > 0) As[1:nobs, (nobs+1):ntot] <- load
  
  Sigma <- solve(diag(ntot) - As) %*% Sym %*% t(solve(diag(ntot) - As))
  
  # Compute K:
# browser()
  K <- solve(cor2cov(cov2cor(Sigma),c(sqrt(diag(S)), rep(1, nLatents) )))
  #  K <- K  
#   K <- cov2cor(K)
  if (!is.positive.definite(K))
  {
    K <- as.matrix(nearPD(K)$mat)
#     warning("Expected covariance matrix is not positive definite")
  }

# browser()

  #   K <- matrix(-0.5,ntot,ntot)
  #   K[1:nobs,1:nobs] <- 0
  #   diag(K) <- 1
  
  #   K <- round(K,2)
  
  #   K <- EBICglasso(cov(cbind(Y,eta)), sampleSize)
  #   K <- as.matrix(forceSymmetric(cbind(rbind(pseudoinverse(resid),t(-load)), rbind(-load,pseudoinverse(r)))))
  
  # K <- as.matrix(forceSymmetric(cbind(rbind(pseudoinverse(resid),-0.5), rbind(rep(-0.5,nobs),1))))
  #   K <- as.matrix(forceSymmetric(cbind(rbind(diag(nrow(S)),t(-load)), rbind(-load,pseudoinverse(r)))))
  # browser()
  #   K <- matrix(-0.5,ntot,ntot)
  #   K[1:nobs,1:nobs] <- 0
  #   diag(K) <- 1
  #   obs <- c(rep(TRUE,nrow(S)), rep(FALSE,nLatents))
  # 
  #   K <- rbind(cbind(2*diag(nobs),-load),cbind(t(-load), 2*diag(ntot - nobs)))
  # rownames(K) <- colnames(K) <- NULL
  # 
  # diag(K) <- diag(K) - min(eigen(pseudoinverse(K))$values)
  
  # #   K[1:nobs,1:nobs] <- 0
  #   diag(K) <- 1
  obs <- c(rep(TRUE,nrow(S)), rep(FALSE,nLatents))
# 
# # Stupid prior:
# K <- matrix(0, ntot, ntot)
# diag(K) <- 1
# if (nLatents > 0) K[1:nobs, (nobs+1):ntot] <- K[ (nobs+1):ntot, 1:nobs] <- -1/nobs
  
  #   is.positive.definite(as.matrix(forceSymmetric(Estep(S, K, obs))))
  ### EM ###
it <- 1
Kold <- K
  repeat
  {
    Sexp <- as.matrix(forceSymmetric(Estep(S, K, obs)))
    #     Sexp <- cor2cov(cov2cor(Sexp),ifelse(obs,sqrt(diag(Sexp)),1))
    K <- as.matrix(forceSymmetric(Mstep(Sexp, obs, rho, lambda)))
    #     qgraph(wi2net(K), layout = "spring")
    
    # If not pos def, shift eigenvalues:
    K <- forcePositive(K)
    
    # Check for convergence:
    if (sum(abs(cov2cor(corpcor::pseudoinverse(Kold)[obs,obs]) - cov2cor(corpcor::pseudoinverse(K)[obs,obs]))) < thr){
      break
    } else {
      it <- it + 1
      if (it > maxit){
        warning("Algorithm did not converge!")
        break
      } else {
        Kold <- K
      }
    }
  }
  
  if (!is.null(colnames(S)))
  {
    colnames(K) <- c(colnames(S),rep(paste0("F",seq_len(nLatents)))) 
  }
  
  if (is.null(rownames(S)))
  {
    rownames(K) <-c(paste0("x",seq_len(ncol(S))),paste0("F",seq_len(nLatents)))
  }

  # Partial correlations:
  pc <- wi2net(K)
  diag(pc) <- 1
  
rownames(pc) <- colnames(pc) <- rownames(K) <- colnames(K)

# Compute psychometric matrices:
Theta <- solve(K[obs, obs])
Lambda <- -Theta %*% K[obs,!obs]
Psi <- solve(K[!obs, !obs] - t(Lambda) %*% K[obs, obs] %*% Lambda)

# Return list mimics glasso:
Res <- list(
  w = corpcor::pseudoinverse(K), # Estimated covariance matrix
  wi = K, # Estimated precision matrix
  pcor = pc, # Estimated partial correlation matrix
  observed = obs, # observed and latents indicator
  niter = it, # Number of iterations used in the algorithm
  lambda = Lambda,
  theta = Theta,
  omega_theta = K[obs, obs],
  psi = Psi
    )

  class(Res) <- "lvglasso"

  return(Res)
}
