matConverter <- function(x){
  mat <- ifelse(x[[1]]$par==0,x[[1]]$est,NA)
  mat
}

lav2lvnet <- function(
  model, # Lavaan model
  data, # data needed for ordering of lambda
  std.lv = TRUE, # DEFAULT IN lvnet!
  lavaanifyOps = list(auto=TRUE,std.lv=std.lv) # lavaanify options
){
  # Variable names:
  varNames <- colnames(data)
  
  # Lavaanify:
  lavMod <- do.call(lavaan::lavaanify,c(list(model), lavaanifyOps))
  
  # Test multiple groups:
  if (length(unique(lavMod$group)) > 1) stop("lvnet only supports single group analysis")
  
  # To model matrices:
  modelMatrices <- semPlot::modelMatrices(lavMod, "mplus")

  # Output list:
  output <- list(
    lambda = matConverter(modelMatrices$Lambda)[varNames,],
    beta = matConverter(modelMatrices$Beta),
    psi = matConverter(modelMatrices$Psi),
    theta = matConverter(modelMatrices$Theta)[varNames,varNames]
  )
  
  return(output)
}