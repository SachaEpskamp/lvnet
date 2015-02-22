# Plotting function:
plot.lvglasso <- function(
  object, # lvglasso object
  plot = c("network","loadings","residcors","residpcors"), # "full" the full network, S" will plot the sparse network between items and "L" the latent loadings
  ask,
  rotation = promax, # A rotation function to be used.
  ...
  ){
    
  if (missing(ask)) ask <- length(plot) > 1
  parOrig <- par()
  par(ask = ask)
  obs <- object$observed
  pcor <- object$pcor
  Res <- list()
  labs <- colnames(pcor)

  if ("network" %in% plot){
    Res$network <- qgraph(pcor, ..., title = "Estimated network", shape = ifelse(obs,"square", "circle"), layout = "spring")
  }
  
  if ("residpcors" %in% plot){
    Res$residpcors <- qgraph(object$pcor[object$observed,object$observed], ..., title = "Estimated residual partial correlations", shape = "square", layout = "spring", repulsion = 0.9)
  }
  
  if ("residcors" %in% plot){
    Res$residcors <- qgraph(cov2cor(object$theta), ..., title = "Estimated residual correlations", shape = "square", layout = "spring", repulsion = 0.9)
  }
  
  if ("loadings" %in% plot){
    fCovs <- object$psi
    load <- object$lambda
    
    # Rotate:
    rot <- rotation(load)
    if (is.matrix(rot)){
      load <- rot
      rotmat <- matrix(1,1,1)
    } else {
      load <- rot$loadings
      rotmat <- rot$rotmat
    }
    fCovs <- solve(rotmat) %*% fCovs %*% t(solve(rotmat))
    
    rownames(load) <- colnames(object$wi)[obs]
    Res$loadings <- qgraph.loadings(load, factorCors = fCovs, ..., title = "Estimated factor loadings", labels =  colnames(object$wi)[obs], model = "reflective")
  }
  
  
  #   
#   if ("full" %in% plot){
#     Res$S <- qgraph(pcor, ..., title = "Full structure", shape = ifelse(obs,"square", "circle"), layout = "spring")
#   }
#   
#   if ("S" %in% plot){
#     Res$S <- qgraph(pcor[obs,obs], ..., title = "Sparse structure", shape = "square", layout = "spring")
#   }
#   
#   if ("L" %in% plot){
#     fCors <- as.matrix(pcor[!obs,!obs])
#     load <- as.matrix(pcor[obs,!obs])
#     
#     rownames(load) <- colnames(object$wi)[obs]
#     
#     Res$L <- qgraph.loadings(load, factorCors = fCors, arrows = FALSE,..., title = "Low-rank structure", labels =  colnames(object$wi)[obs])
#   }
#   
  
  
  
  invisible(Res)
}
