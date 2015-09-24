
# Computes optimal glasso network based on EBIC:
EBIClvglasso <- function(
  S, # Sample cov
  n, # Sample size
  nLatents, # Number of latents
  gamma = 0.5, # EBIC parameter
  nRho = 100,
  lambda, 
  ... # lvglasso arguments
){
  if (missing(lambda)) lambda <- NULL
  
  # If nLatents is vector, do this function for every latent:
  if (length(nLatents) > 1){
    Resses <- lapply(nLatents,function(nl)EBIClvglasso(S, n, nl, gamma, lambda, ...))
    opt <- which.max(sapply(Resses,'[[','ebic'))
    return(Resses[[opt]])
  }
  
  rho.max = max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  rho.min = rho.max/100
  rho = exp(seq(log(rho.min), log(rho.max), length = nRho))
  
  lvglas_res <- lapply(rho, function(r)try(lvglasso(S, nLatents, r,lambda =  lambda, ...)))

  failed <- sapply(lvglas_res,is,"try-res")
  # Likelihoods:
  EBICs <- sapply(lvglas_res[!failed],function(res){
    C <- solve(res$w[res$observed,res$observed])
    qgraph:::EBIC(S, C, n, gamma, E = sum(res$wi[lower.tri(res$wi, diag = TRUE)] != 0))
  })
  
  # Smalles EBIC:
  opt <- which.min(EBICs)
  
  Res <- lvglas_res[!failed][[opt]]
  Res$rho <- rho[!failed][opt]
  Res$ebic <- EBICs[opt]
  
  # Return 
  return(Res)
}


