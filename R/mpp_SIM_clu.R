###############
# mpp_SIM_clu #
###############

# SIM function for general procedure and CV. This function is an exact copy from
# mpp_SIM function. The only difference is that it keep the possibility to
# pass to the function a cluster object that is already defined.

mpp_SIM_clu <- function(mppData, trait = 1, Q.eff = "cr", 
                    plot.gen.eff = FALSE, parallel = FALSE, cluster = NULL) {
  
  # 1. Check data format and arguments
  ####################################
  
  check.model.comp(mppData = mppData, trait = trait, Q.eff = Q.eff,
                   VCOV = 'h.err', plot.gen.eff = plot.gen.eff, fct = "SIM")
  
  # 2. Form required elements for the analysis
  ############################################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  ### 2.3 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  vect.pos <- 1:dim(mppData$map)[1]
  
  # 3. computation of the SIM profile (genome scan)
  #################################################
  
  if (parallel) {
    
    log.pval <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelSIM,
                          mppData = mppData, trait = t_val,
                          cross.mat = cross.mat,
                          Q.eff = Q.eff, VCOV = 'h.err',
                          plot.gen.eff = plot.gen.eff)
    
  } else {
    
    log.pval <- lapply(X = vect.pos, FUN = QTLModelSIM,
                       mppData = mppData, trait = t_val, cross.mat = cross.mat,
                       Q.eff = Q.eff, VCOV = 'h.err', plot.gen.eff = plot.gen.eff)
    
  }
  
  log.pval <- t(data.frame(log.pval))
  if(plot.gen.eff){log.pval[is.na(log.pval)] <- 1}
  log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
  log.pval[is.na(log.pval[, 1]), 1] <- 0
  
  
  # 4. form the results
  #####################
  
  SIM <- data.frame(mppData$map, log.pval)
  
  if(plot.gen.eff){
    
    if(Q.eff == "cr"){ Qeff_names <- unique(mppData$cross.ind)
    
    } else { Qeff_names <- mppData$parents }
    
    colnames(SIM)[5:dim(SIM)[2]] <- c("log10pval", Qeff_names)
    
  } else {colnames(SIM)[5] <- "log10pval"}
  
  
  class(SIM) <- c("QTLprof", "data.frame")
  
  ### 4.1: Verify the positions for which model could not be computed
  
  if(sum(SIM$log10pval == 0) > 0) {
    
    if (sum(SIM$log10pval) == 0){
      
      warning("the computation of the QTL model failled for all positions")
      
    } else {
      
      list.pos <- mppData$map[(SIM$log10pval == 0), 1]
      
      prob_pos <- paste(list.pos, collapse = ", ")
      
      message("the computation of the QTL model failed for the following ",
              "positions: ", prob_pos,
              ". This could be due to singularities or function issues")
      
    }
    
  }
  
  return(SIM)
  
}