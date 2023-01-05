#################
# Qprof_process #
#################

# function to process the results after a SIM or CIM computation

Qprof_process <- function(mppData, Q.eff, log.pval, nEnv){
  
  Qprof <- data.frame(mppData$map, log.pval)
  
  if(Q.eff == "cr"){
    
    Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.cr)
    Qeff_names <- paste0(rep(unique(mppData$cross.ind), nEnv), Env_name)
    
    colnames(Qprof)[5] <- "log10pval"
    colnames(Qprof)[(6 + nEnv):dim(Qprof)[2]] <- Qeff_names
    
  } else if (Q.eff == "anc") {
    
    Env_name <- rep(paste0("_E", 1:nEnv), each = mppData$n.par)
    Qeff_names <- paste0(rep(mppData$parents, nEnv), Env_name)
    
    colnames(Qprof)[5] <- "log10pval"
    colnames(Qprof)[(6 + nEnv):dim(Qprof)[2]] <- Qeff_names
    
  } else { # par and bi-allelic no modif for gen effects names
    
    colnames(Qprof)[5] <- "log10pval"
    
  }
  
  class(Qprof) <- c("QTLprof", "data.frame")
  
  # 5.1: Verify the positions for which model could not be computed
  
  if(sum(Qprof$log10pval == 0) > 0) {
    
    if (sum(Qprof$log10pval) == 0){
      
      message("The computation of the QTL models failled for all positions.")
      
    } else {
      
      list.pos <- mppData$map[(Qprof$log10pval == 0), 1]
      
      text <- paste("The computation of the QTL model failed for the following",
                    "positions: ", paste(list.pos, collapse = ", "),
                    ". This could be due to singularities or function issues.")
      
      message(text)
      
    }
    
  }
  
  return(Qprof)
  
}