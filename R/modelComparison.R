rimCompare <- function(...){
 
  dots <- list(...)
  
  # Combine rim objects:
  rimObjects <- dots[sapply(dots,is,"rim")]
  if (is.null(names(rimObjects))){
    names(rimObjects) <- paste("Model",seq_along(rimObjects))
  }
  
  # Create rows:
  DF <- rbind(data.frame(Df = 0, AIC = NA, 
                         BIC = NA,
                         Chisq = 0),
              do.call(rbind,   lapply(rimObjects,function(x){
                data.frame(Df = x$fitMeasures$df, AIC = x$fitMeasures$aic, 
                           BIC = x$fitMeasures$bic,
                           Chisq = x$fitMeasures$chisq)
              }))
              )
  
  rownames(DF) <- c("Saturated",names(rimObjects))
  
  DF[['Chisq diff']] <- c(NA,diff(DF[['Chisq']]))
  DF[['Df diff']] <- c(NA,diff(DF[['Df']]))
  DF[['Pr(>Chisq)']] <- pchisq( DF[['Chisq diff']], DF[['Df diff']], lower.tail=FALSE)

  return(DF)
}

anova.rim <- function(object,...) rimCompare(object,...)