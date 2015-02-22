print.rim <- summary.rim <- function(object, include = c('input','chisq','infcrit','fitindices','rmsea','parests')){
  
  cat("========== RIM ANALYSIS RESULTS ========== ")
  
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
      "\n\tChi-square:\t\t",object$fitMeasures$chisq,
      "\n\tDF:\t\t\t",object$fitMeasures$df,
      "\n\tp-value:\t\t",object$fitMeasures$pvalue
    )
  }

  if ('infcrit' %in% include){
    cat(
      "\n\nInformation criteria:",
      "\n\tAIC:\t\t\t",object$fitMeasures$aic,
      "\n\tBIC:\t\t\t",object$fitMeasures$bic,
      "\n\tAdjusted BIC:\t\t",object$fitMeasures$bic2
    )
  }
  
  if ('fitindices' %in% include){
    cat(
      "\n\nFit indices:",
      "\n\tCFI:\t\t\t",object$fitMeasures$cfi,
      "\n\tNFI:\t\t\t",object$fitMeasures$nfi,
      "\n\tTLI:\t\t\t",object$fitMeasures$tli,
      "\n\tRFI:\t\t\t",object$fitMeasures$rfi,
      "\n\tIFI:\t\t\t",object$fitMeasures$ifi,
      "\n\tRNI:\t\t\t",object$fitMeasures$rni,
      "\n\tRMR:\t\t\t",object$fitMeasures$rmr,
      "\n\tSRMR:\t\t\t",object$fitMeasures$srmr
    )
    
  }
  
  if ('rmsea' %in% include){
    cat(
      "\n\nRMSEA:",
      "\n\tRMSEA:\t\t\t",object$fitMeasures$rmsea,
      "\n\t90% CI lower bound:\t",object$fitMeasures$rmsea.ci.lower,
      "\n\t90% CI upper bound:\t",object$fitMeasures$rmsea.ci.upper,
      "\n\tp-value:\t\t",object$fitMeasures$rmsea.pvalue
    )
  }

  if ('parests' %in% include){
    cat("\n\nParameter estimates:\n")
    sum <- summary(object$mxResults$model)
    parMat <- sum$parameters
    print(parMat[,c('matrix','row','col','Estimate','Std.Error')])
  }
  
}