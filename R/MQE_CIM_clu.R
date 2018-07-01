###############
# mpp_CIM_clu #
###############

# CIM function for general procedure and forward. This function is an exact
# copy from MQE_CIM function. The only difference is that it keep the
# possibility to pass to the function a cluster object that is already defined.

MQE_CIM_clu <- function(mppData = NULL, trait = 1, Q.eff = "cr", VCOV = "h.err",
                    cofactors = NULL, cof.Qeff, chg.Qeff = FALSE, window = 20,
                    parallel = FALSE, cluster = NULL){
  
  # 1. Check data format and arguments
  ####################################
  
  check.MQE(mppData = mppData, Q.eff = Q.eff, trait = trait,
            VCOV = VCOV, cofactors = cofactors, cof.Qeff = cof.Qeff,
            n.cores = 1, fct = "CIM")
  
  # 2. Form required elements for the analysis
  ############################################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  ### 2.3 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.4 Formation of the list of cofactors
  
  # order the list of cofactors and the corresponding Q.eff
  
  cof.pos <- vapply(X = cofactors,
                    FUN = function(x, mppData) which(mppData$map[, 1] == x),
                    FUN.VALUE = numeric(1), mppData = mppData)
  
  cof.ord <- data.frame(cofactors, cof.Qeff, cof.pos,
                        stringsAsFactors = FALSE)
  
  cof.ord <- cof.ord[order(cof.pos), ]
  
  cofactors <- cof.ord[, 1]; cof.Qeff <- cof.ord[, 2]
  
  cof.pos <- which(mppData$map[, 1] %in% cofactors)
  
  
  cof.list <- mapply(FUN = inc_mat_QTL, x = cof.pos, Q.eff = cof.Qeff,
                     MoreArgs = list(mppData = mppData,order.MAF = TRUE),
                     SIMPLIFY = FALSE)
  
  ### 2.5 Formation of the genome-wide and cofactors partition
  
  vect.pos <- 1:dim(mppData$map)[1]
  
  # 2.5.1 cofactor partition tells if the cofactor should be included or
  # not in the model at each position.
  
  cofactors2 <- mppData$map[cof.pos, c(2,4)]
  
  test.cof <- function(x, map, window) {
    
    t1 <- map$chr == as.numeric(x[1])
    t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
    !(t1 & t2)
    
  }
  
  cof.part <- apply(X = cofactors2, MARGIN = 1, FUN = test.cof,
                    map = mppData$map, window = window)
  
  # make a QTL effect partition to change type of QTL effect of the tested pos.
  
  if(chg.Qeff) {
    
    cof.part2 <- (!cof.part)*1
    
    Qeff.partition <- function(x, Q.eff, cof.Qeff){
      if(sum(x) == 0){Q.eff } else { cof.Qeff[max(which(x == 1))] }
    }
    
    Qeff.part <- apply(X = cof.part2, MARGIN = 1, FUN = Qeff.partition,
                       Q.eff = Q.eff, cof.Qeff = cof.Qeff)
    
  } else {Qeff.part <- rep(Q.eff, length(vect.pos))}
  
  
  # 3. computation of the CIM profile (genome scan)
  #################################################
  
  if (parallel) {
    
    log.pval <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelCIM_MQE,
                          mppData = mppData, trait = t_val, cross.mat = cross.mat, 
                          Qeff.part = Qeff.part, VCOV = VCOV,
                          cof.list = cof.list, cof.part = cof.part)
    
    
  } else {
    
    log.pval <- lapply(X = vect.pos, FUN = QTLModelCIM_MQE,
                       mppData = mppData, trait = t_val, cross.mat = cross.mat,
                       Qeff.part = Qeff.part, VCOV = VCOV,
                       cof.list = cof.list, cof.part = cof.part)
  }
  
  log.pval <- t(data.frame(log.pval))
  log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
  log.pval[is.na(log.pval[, 1]), 1] <- 0
  
  # 4. form the results
  #####################
  
  CIM <- data.frame(mppData$map, log.pval)
  colnames(CIM)[5] <- "log10pval"
  class(CIM) <- c("QTLprof", "data.frame")
  
  ### 4.1: Verify the positions for which model could not be computed
  
  if(sum(CIM$log10pval == 0) > 0) {
    
    if (sum(CIM$log10pval) == 0){
      
      warning("the computation of the QTL models failled for all positions")
      
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