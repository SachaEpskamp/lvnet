# Plotting function:
plot.lvglasso <- function(
  x, # lvglasso x
  plot = c("network","loadings","residcors","residpcors"), # "full" the full network, S" will plot the sparse network between items and "L" the latent loadings
  ask,
  rotation = promax, # A rotation function to be used.
  ...
  ){
    
  if (missing(ask)) ask <- length(plot) > 1
  parOrig <- par()
  par(ask = ask)
  obs <- x$observed
  pcor <- x$pcor
  Res <- list()
  labs <- colnames(pcor)

  if ("network" %in% plot){
    Res$network <- qgraph(pcor, ..., title = "Estimated network", shape = ifelse(obs,"square", "circle"), layout = "spring")
  }
  
  if ("residpcors" %in% plot){
    Res$residpcors <- qgraph(x$pcor[x$observed,x$observed], ..., title = "Estimated residual partial correlations", shape = "square", layout = "spring", repulsion = 0.9)
  }
  
  if ("residcors" %in% plot){
    Res$residcors <- qgraph(cov2cor(x$theta), ..., title = "Estimated residual correlations", shape = "square", layout = "spring", repulsion = 0.9)
  }
  
  if ("loadings" %in% plot){
    fCovs <- x$psi
    load <- x$lambda
    
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
    
    rownames(load) <- colnames(x$wi)[obs]
    Res$loadings <- qgraph.loadings(load, factorCors = fCovs, ..., title = "Estimated factor loadings", labels =  colnames(x$wi)[obs], model = "reflective")
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
#     rownames(load) <- colnames(x$wi)[obs]
#     
#     Res$L <- qgraph.loadings(load, factorCors = fCors, arrows = FALSE,..., title = "Low-rank structure", labels =  colnames(x$wi)[obs])
#   }
#   
  
  
  
  invisible(Res)
}
