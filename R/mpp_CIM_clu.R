###############
# mpp_CIM_clu #
###############

# CIM function for general procedure and CV. This function is an exact copy from
# mpp_CIM function. The only difference is that it keep the possibility to
# pass to the function a cluster object that is already defined.

mpp_CIM_clu <- function(mppData, trait = 1, Q.eff = "cr",
                    cofactors = NULL,  window = 20, plot.gen.eff = FALSE,
                    parallel = FALSE, cluster = NULL)
{
  
  # 1. Check data format and arguments
  ####################################
  
  check.model.comp(mppData = mppData, trait = trait, Q.eff = Q.eff,
                   VCOV = 'h.err', plot.gen.eff = plot.gen.eff,
                   cofactors = cofactors, fct = "CIM")
  
  # 2. Form required elements for the analysis
  ############################################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  ### 2.2 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.3 Formation of the list of cofactors
  
  if(is.character(cofactors)){
    
    cof.pos <- which(mppData$map[, 1] %in% cofactors)
    
  } else {
    
    cof.pos <- which(mppData$map[, 1] %in% cofactors[, 1])
    
  }
  
  cof.list <- lapply(X = cof.pos, FUN = inc_mat_QTL, mppData = mppData,
                     Q.eff = Q.eff, order.MAF = TRUE)
  
  ### 2.4 Formation of the genome-wide and cofactors partition
  
  vect.pos <- 1:dim(mppData$map)[1]
  
  # cofactor partition tells if the cofactor should be included or
  # not in the model at each position.
  
  if (is.character(cofactors)){
    
    cofactors2 <- mppData$map[mppData$map[, 1] %in% cofactors, c(2, 4)]
    
  } else { cofactors2 <- cofactors[, c(2, 4)] }
  
  test.cof <- function(x, map, window) {
    
    t1 <- map$chr == as.numeric(x[1])
    t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
    !(t1 & t2)
    
  }
  
  cof.part <- apply(X = cofactors2, MARGIN = 1, FUN = test.cof,
                    map = mppData$map, window = window)
  
  
  # 3. computation of the CIM profile (genome scan)
  #################################################
  
  if (parallel) {
    
    log.pval <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelCIM,
                          mppData = mppData, trait = t_val, cross.mat = cross.mat,
                          Q.eff = Q.eff, VCOV = 'h.err', cof.list = cof.list,
                          cof.part = cof.part, plot.gen.eff = plot.gen.eff)
    
  } else {
    
    log.pval <- lapply(X = vect.pos, FUN = QTLModelCIM,
                       mppData = mppData, trait = t_val, cross.mat = cross.mat,
                       Q.eff = Q.eff, VCOV = 'h.err', cof.list = cof.list,
                       cof.part = cof.part, plot.gen.eff = plot.gen.eff)
    
  }
  
  
  log.pval <- t(data.frame(log.pval))
  if(plot.gen.eff ){log.pval[is.na(log.pval)] <- 1}
  log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
  log.pval[is.na(log.pval[, 1]), 1] <- 0
  
  # 4. form the results
  #####################
  
  CIM <- data.frame(mppData$map, log.pval)
  
  
  if(plot.gen.eff){
    
    if(Q.eff == "cr"){ Qeff_names <- unique(mppData$cross.ind)
    
    } else { Qeff_names <- mppData$parents }
    
    colnames(CIM)[5:dim(CIM)[2]] <- c("log10pval", Qeff_names)
    
  } else {colnames(CIM)[5] <- "log10pval"}
  
  class(CIM) <- c("QTLprof", "data.frame")
  
  ### 4.1: Verify the positions for which model could not be computed
  
  if(sum(CIM$log10pval == 0) > 0) {
    
    if (sum(CIM$log10pval) == 0){
      
      warning("the computation of the QTL model failled for all positions")
      
    } else {
      
      list.pos <- mppData$map[(CIM$log10pval == 0), 1]
      
      prob_pos <- paste(list.pos, collapse = ", ")
      
      message("the computation of the QTL model failed for the following ",
              "positions: ", prob_pos,
              ". This could be due to singularities or function issues")
      
    }
    
  }
  
  return(CIM)
  
}