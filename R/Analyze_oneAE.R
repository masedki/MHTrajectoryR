Analyze_oneAE <- function(ae,
                          drug,
                          maxit,
                          alpha,
                          nbinit)
{
  if (sum(ae)>1){
    datau <- uniquecombs(cbind(ae,drug[,-which(colSums(drug[which(ae==1),])==0)]))
    w <- table(attr(datau,"index"))
    x <- datau[,-1]
    y <- datau[,1]
    
    drugs.names <- colnames(drug[,-which(colSums(drug[which(ae==1),])==0)])
    
    if (ncol(datau) > 12){
      cat("Metropolis-Hastings is used since many drugs are considered\n")
      
      maxitlist <- as.list(rep(maxit,nbinit))
      nb.cpus <- min(detectCores(all.tests = FALSE, logical = FALSE) , nbinit)
      if(Sys.info()["sysname"] == "Windows")
      {
        print("TODO")
      }
      else
      {
        results <- mclapply(X = maxitlist,
                            FUN = MHLogisticw,
                            y = y,
                            x = x,
                            w = w,
                            drugs.names = drugs.names,
                            alpha = alpha,
                            mc.cores = nb.cpus,
                            mc.preschedule = TRUE,
                            mc.cleanup = TRUE
        )
        tmpres <- rep(NA, length(results))
        for (k in 1:length(results))
          tmpres[k] <- results[[k]]$best.model.bic
        if (any(is.na(tmpres)))
          tmpres[which(is.na(tmpres))] <- -Inf

        results <- list(best=results[[which(tmpres==max(tmpres))[1]]], all=results)
      }      
      results$signals <- FindSignals(results)
    }else{
      cat("Exhaustive approach is used since few drugs are considered\n")
      results <- list(best=ExhaustiveLogisticw(y = y, x = as.matrix(x), w=w, drugs.names = drugs.names))
      results$signals <- FindSignals(results)
    }
  }else{
    cat("The adverse event is notified only one time, so the estimation is not doable!\n")
    results <- list(best=list(model=NULL))
  }
  
  return(results)  
}


FindSignals <- function(results){
  output <- data.frame(ATC=as.character(results$best$drugs.names[ which(results$best$best.model == 1) ]),
                       beta=results$best$best.fit.coef[-1], 
                       stringsAsFactors = F, 
                       row.names = NULL)
  output <- output[order(output$beta, decreasing = T),]
  output <- output[which(output$beta>0),]
  return(output)
}
  
