lvnetCompare <- function(...){
 
  dots <- list(...)
  
  # Combine lvnet objects:
  lvnetObjects <- dots[sapply(dots,is,"lvnet")]
  if (is.null(names(lvnetObjects))){
    names(lvnetObjects) <- paste("Model",seq_along(lvnetObjects))
  }
  
  # Create rows:
  DF <- rbind(data.frame(Df = 0, AIC = NA, 
                         BIC = NA,
                         Chisq = 0),
              do.call(rbind,   lapply(lvnetObjects,function(x){
                data.frame(Df = x$fitMeasures$df, AIC = x$fitMeasures$aic, 
                           BIC = x$fitMeasures$bic,
                           Chisq = x$fitMeasures$chisq)
              }))
              )
  
  rownames(DF) <- c("Saturated",names(lvnetObjects))
  
  DF[['Chisq diff']] <- c(NA,abs(diff(DF[['Chisq']])))
  DF[['Df diff']] <- c(NA,abs(diff(DF[['Df']])))
  DF[['Pr(>Chisq)']] <- pchisq( DF[['Chisq diff']], DF[['Df diff']], lower.tail=FALSE)

  return(DF)
}

anova.lvnet <- function(object,...) lvnetCompare(object,...)
