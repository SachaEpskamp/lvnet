summary.lvnet <- function(object, include = c('input','chisq','infcrit','fitindices','rmsea','parests'), digits = 3,...){
  
  cat("========== lvnet ANALYSIS RESULTS ========== ")
  
  if ('input' %in% include){
    cat(
      "\n\nInput:",
      "\n\tModel:\t\t\t",object$model,
      "\n\tNumber of manifests:\t",ncol(object$sampleStats$covMat),
      "\n\tNumber of latents:\t",ncol(object$matrices$lambda),
      "\n\tNumber of parameters:\t",object$fitMeasures$npar,
      "\n\tNumber of observations\t",object$sampleStats$sampleSize
    )
  }

  if ('chisq' %in% include){
    cat(
      "\n\nTest for exact fit:",
      "\n\tChi-square:\t\t",round(object$fitMeasures$chisq,digits),
      "\n\tDF:\t\t\t",round(object$fitMeasures$df,digits),
      "\n\tp-value:\t\t",round(object$fitMeasures$pvalue,digits)
    )
  }

  if ('infcrit' %in% include){
    cat(
      "\n\nInformation criteria:",
      "\n\tAIC:\t\t\t",round(object$fitMeasures$aic,digits),
      "\n\tBIC:\t\t\t",round(object$fitMeasures$bic,digits),
      "\n\tAdjusted BIC:\t\t",round(object$fitMeasures$bic2,digits),
      "\n\tExtended BIC:\t\t",round(object$fitMeasures$ebic,digits)
    )
  }
  
  if ('fitindices' %in% include){
    cat(
      "\n\nFit indices:",
      "\n\tCFI:\t\t\t",round(object$fitMeasures$cfi,digits),
      "\n\tNFI:\t\t\t",round(object$fitMeasures$nfi,digits),
      "\n\tTLI:\t\t\t",round(object$fitMeasures$tli,digits),
      "\n\tRFI:\t\t\t",round(object$fitMeasures$rfi,digits),
      "\n\tIFI:\t\t\t",round(object$fitMeasures$ifi,digits),
      "\n\tRNI:\t\t\t",round(object$fitMeasures$rni,digits),
      "\n\tRMR:\t\t\t",round(object$fitMeasures$rmr,digits),
      "\n\tSRMR:\t\t\t",round(object$fitMeasures$srmr,digits)
    )
    
  }
  
  if ('rmsea' %in% include){
    cat(
      "\n\nRMSEA:",
      "\n\tRMSEA:\t\t\t",round(object$fitMeasures$rmsea,digits),
      "\n\t90% CI lower bound:\t",round(object$fitMeasures$rmsea.ci.lower,digits),
      "\n\t90% CI upper bound:\t",round(object$fitMeasures$rmsea.ci.upper,digits),
      "\n\tp-value:\t\t",round(object$fitMeasures$rmsea.pvalue,digits)
    )
  }

  if ('parests' %in% include){
    cat("\n\nParameter estimates:\n")
    sum <- summary(object$mxResults$model)
    parMat <- sum$parameters
    parMat[['Estimate']] <- round(parMat[['Estimate']],digits)
    parMat[['Std.Error']] <- round(parMat[['Std.Error']],digits)
    print(parMat[,c('matrix','row','col','Estimate','Std.Error')], row.names=FALSE)
  }
  
}

print.lvnet <- function(x,...) summary(x,...)
