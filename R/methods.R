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
    
    # if theta_inverse in name, compute omega_theta and delta_theta:
    if (any(grepl("theta_inverse",parMat$name))){
      
      theta_inverse <- parMat %>% filter_(~matrix == "theta_inverse")
      # Compute variance:
      theta_inverse <- theta_inverse %>% mutate_(var = ~Std.Error^2) 
      # Obtain diagonals:
      diag_theta_inverse <- theta_inverse %>% filter(row==col) %>% select_(~row,diag=~Estimate)
      # Add to matrix:
      theta_inverse <- theta_inverse %>% left_join(diag_theta_inverse %>% rename_(diagI = ~diag), by = "row") %>%
        left_join(diag_theta_inverse %>% rename_(diagJ = ~diag), by = c("col" = "row"))
      
      # Compute pcors and standard deviations:
      theta_inverse <- theta_inverse %>% mutate_(
        Estimate = ~-Estimate / (sqrt(diagI)*sqrt(diagJ)),
        Std.Error = ~(-1/ (sqrt(diagI)*sqrt(diagJ)))^2 * Std.Error
        )
      theta_inverse$matrix <- "omega_theta"
      theta_inverse <- theta_inverse %>% filter(row != col)
      
      # compute delta_theta:
      delta_theta <- parMat %>% filter_(~matrix == "theta_inverse") %>%
        filter_(~row==col) %>% mutate_(Estimate = ~1/sqrt(Estimate))
    delta_theta$matrix <- "delta_theta"
     
    
    parMat <- dplyr::bind_rows(parMat[parMat$matrix != "theta_inverse",],theta_inverse,delta_theta) 
    }
    
    parMat[['Estimate']] <- round(parMat[['Estimate']],digits)
    parMat[['Std.Error']] <- round(parMat[['Std.Error']],digits)
    print.data.frame(parMat[,c('matrix','row','col','name','Estimate')], row.names=FALSE)
  }
  
}

print.lvnet <- function(x,...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  if (name=="x") name <- "object"
  
  cat("\nlvnet estimation completed.:\n",
      paste0("\t- Chi-square (",x$fitMeasures$df,") = ",round(x$fitMeasures$chisq,2),", p = ",round(x$fitMeasures$pvalue,2),"\n"),
      paste0("\t- RMSEA = ",round(x$fitMeasures$rmsea,2)," (95% CI: ",round(x$fitMeasures$rmsea.ci.lower,2)," - ",round(x$fitMeasures$rmsea.ci.upper,2),")\n")
    )
  
  ### Tips
  cat("\n",
      paste0("Use summary(",name,") to inspect more fitmeasures and parameter estimates (see ?summary.lvnet)"),
      "\n",
      paste0("Use plot(",name,") to plot estimated networks and factor structures (see ?plot.lvnet)"),
      "\n",
      paste0("Use lvnetCompare(object1, object2) to compare lvnet models (see ?lvnetCompare)")
  )
}



# lasso and search print and summary:
summary.lvnetSearch <- summary.lvnetLasso <- function(object,...){
  summary(object$best,...)
}

print.lvnetLasso <- function(x,...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  if (name=="x") name <- "object"

  cat("\nlvnetLasso completed:\n",
      paste0("\t- Criterion used: ",x$criterion,"\n"),
      paste0("\t- # tuning parameters: ",length(x$tuning),"\n"),
      paste0("\t- # tuning parameter range: ",min(x$tuning)," - ",max(x$tuning), "\n")
  )

  cat("\nBest model:\n",
      paste0("\t- Tuning: ",round(x$tuning[x$bestID],2)),"\n",
      paste0("\t- Chi-square (",x$best$fitMeasures$df,") = ",round(x$best$fitMeasures$chisq,2),", p = ",round(x$best$fitMeasures$pvalue,2),"\n"),
      paste0("\t- RMSEA = ",round(x$best$fitMeasures$rmsea,2)," (95% CI: ",round(x$best$fitMeasures$rmsea.ci.lower,2)," - ",round(x$best$fitMeasures$rmsea.ci.upper,2),")\n")
  )
  
  ### Tips
  cat("\n",
      paste0("Best model is stored under ",name,"$best"),
      "\n",
      paste0("Use summary(",name,") to inspect best model (see ?summary.lvnet)"),
      "\n",
      paste0("Use plot(",name,") to plot best model (see ?plot.lvnet)"),
      "\n",
      paste0("Use lvnetCompare(object1, object2) to compare lvnet models (see ?lvnetCompare)")
  )
}



print.lvnetSearch <- function(x,...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  if (name=="x") name <- "object"
  
  cat("\nlvnetSearch completed:\n",
      paste0("\t- Criterion used: ",x$criterion,"\n"),
      paste0("\t- # iterations: ",length(x$nIter),"\n")
  )
  
  cat("\nBest model:\n",
      paste0("\t- Chi-square (",x$best$fitMeasures$df,") = ",round(x$best$fitMeasures$chisq,2),", p = ",round(x$best$fitMeasures$pvalue,2),"\n"),
      paste0("\t- RMSEA = ",round(x$best$fitMeasures$rmsea,2)," (95% CI: ",round(x$best$fitMeasures$rmsea.ci.lower,2)," - ",round(x$best$fitMeasures$rmsea.ci.upper,2),")\n")
  )
  
  ### Tips
  cat("\n",
      paste0("Best model is stored under ",name,"$best"),
      "\n",
      paste0("Use summary(",name,") to inspect best model (see ?summary.lvnet)"),
      "\n",
      paste0("Use plot(",name,") to plot best model (see ?plot.lvnet)"),
      "\n",
      paste0("Use lvnetCompare(object1, object2) to compare lvnet models (see ?lvnetCompare)")
  )
}




# Plot method:
plot.lvnetSearch <- plot.lvnetLasso <- function(x,...) plot.lvnet(x$best,...)

plot.lvnet <- function(x, 
        what = c("factorStructure","residual","latent"),
        partial, # Should partial correlations be plotted? defaults to TRUE if partial correlations are modeled
        layout = "circle",
        ... # Arguments sent to qgraph
        ){
  
  what <- match.arg(what)
  
  if (what == "factorStructure"){
     # get Lambda:
    Lambda <- x$matrices$lambda
    
    # Get residual variances:
    residVar <- diag(x$matrices$theta)
    
    # Is omega_psi estimated?
    if (missing(partial)){
      partial <- any(grepl("omega_psi",names(x$mxResults$model@output$matrices)))
    }
    
    if (!partial){
      graph <- qgraph::qgraph.loadings(Lambda, model = "reflective", factorCors = x$matrices$psi, resid=residVar,layout = layout,  ...)
    } else {
      graph <- qgraph::qgraph.loadings(Lambda, model = "reflective", factorCors = x$matrices$omega_psi, DoNotPlot=TRUE, layout = layout, ...)
      graph$Edgelist$bidirectional[] <- FALSE
      graph$Edgelist$directed[graph$Edgelist$to > nrow(Lambda)] <- FALSE
      plot(graph)
    }
  } else  if (what == "latent"){
    
    # Partial?
    if (missing(partial)){
      partial <- any(grepl("omega_psi",names(x$mxResults$model@output$matrices)))
    }
    
    if (partial){
      graph <- qgraph::qgraph(x$matrices$omega_psi,layout = layout,...)
    } else {
      graph <- qgraph::qgraph(cov2cor(x$matrices$psi),directed=TRUE,bidirectional=TRUE,layout = layout,...)
    }
    
  } else if (what == "residual"){
    
    # Partial?
    if (missing(partial)){
      partial <- any(grepl("theta_inverse",names(x$mxResults$model@output$matrices)))
    }
    
    if (partial){
      graph <- qgraph::qgraph(x$matrices$omega_theta,layout = layout,...)
    } else {
      graph <- qgraph::qgraph(cov2cor(x$matrices$theta),directed=TRUE,bidirectional=TRUE,layout = layout,...)
    }
    
  } else stop("Graph not supported.")
  
  invisible(graph)
}





