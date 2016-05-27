lassoSelect <- function(
  object, # lvnetLasso object
  select, # an R expression.
  minimize = TRUE,
  refit = TRUE,
  lassoTol = 1e-4
){
  stopifnot(is(object,"lvnetLasso"))

  # Table of fit measures:
  fitTable <- as.data.frame(do.call(rbind,lapply(object$modList,function(x)unlist(x$fitMeasures))))

  # Eval selection:
  fits <- unlist(eval(substitute(select), envir = fitTable))

  if (minimize){
    best <- which.min(fits)
  } else {
    best <- which.max(fits)
  }
 
  if (refit){
    newMod <- lapply(object$lassoMatrix, function(m){
      mat <- object$modList[[best]]$matrices[[m]]
      ifelse(abs(mat) > lassoTol,NA,0)
    })
    names(newMod) <- object$lassoMatrix
    
    args <- object$args
    for (i in seq_along(object$lassoMatrix)){
      args[[object$lassoMatrix[[i]]]] <- newMod[[object$lassoMatrix[[i]]]] 
    }
    
    bestModel <- do.call(lvnet,c(object$args))
  } else {
    bestModel <- object$modList[[best]]
  }
  
  object$best <- bestModel  
  return(object)
}