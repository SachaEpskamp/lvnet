lvnetCompare <- function(...){
 
  dots <- list(...)

  # Combine lvnet objects:
  # LASSO objectS:
  lvnetLasObj <- dots[sapply(dots,is,"lvnetLasso")]
  lvnetObjects <- dots[sapply(dots,is,"lvnet")]
  
  if (length(lvnetLasObj) > 0){
    for (i in seq_along(lvnetLasObj)){
      lvnetObjects <- c(lvnetObjects,lvnetLasObj[[i]]$modList)
    }

  }

  if (length(lvnetObjects) == 0){
    stop("No 'lvnet' models in input.")
  }
  
  if (is.null(names(lvnetObjects))){
    names(lvnetObjects) <- paste("Model",seq_along(lvnetObjects))
  }
  
  # Create rows:
  DF <- rbind(data.frame(Df = 0, AIC = NA, 
                         BIC = NA,
                         EBIC = NA,
                         Chisq = 0),
              do.call(rbind,   lapply(lvnetObjects,function(x){
                data.frame(Df = x$fitMeasures$df, AIC = x$fitMeasures$aic, 
                           BIC = x$fitMeasures$bic,
                           EBIC = x$fitMeasures$ebic,
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
