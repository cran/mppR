###########
# MQE_CIM #
###########

# Composite interval maping for MQE model
# 
# Compute a QTL profile with cofactors position having different type of QTL
# effects.
# 
# \strong{WARNING!} The computation of random pedigree models
# (\code{VCOV = "pedigree" and "ped_cr.err"}) can sometimes fail. This could be
# due to singularities due to a strong correlation between the QTL term(s) and 
# the polygenic term. This situation can appear in the parental model.
# the error can also sometimes come from the \code{asreml()} function. From
# our experience, in that case, trying to re-run the function one or two times
# allow to obtain a result.
# 
# @param mppData An object of class \code{mppData}.
# 
# @param trait \code{Numerical} or \code{character} indicator to specify which
# trait of the \code{mppData} object should be used. Default = 1.
# 
# @param Q.eff \code{Character} expression indicating the type of QTL effect at
# the tested position. Possibility to choose among: "cr", "par", "anc" or
# "biall". For details look at \code{\link{mpp_SIM}}. Default = "cr".
#
# @param VCOV \code{Character} expression defining the type of variance
# covariance structure used: 1) "h.err" for an homogeneous variance residual term
# (HRT) linear model; 2) "h.err.as" for a HRT model fitted by REML using
# \code{ASReml-R}; 3) "cr.err" for a cross-specific variance residual terms
# (CSRT) model; 4) "pedigree" for a random pedigree term and HRT model;
# and 5) "ped_cr.err" for random pedigree and CSRT model.
# For more details see \code{\link{mpp_SIM}}. Default = "h.err".
# 
# @param cofactors Vector of \code{character} marker or in between marker
# positions names. Default = NULL.
# 
# @param cof.Qeff \code{Character} vector indicating for each cofactor position
# the type of QTL effect from: "cr", "par", "anc" and "biall".
# 
# @param chg.Qeff \code{Logical} value. If \code{chg.Qeff = TRUE}.
# The type of QTL effect of the tested position will change for the one
# specified in \code{cof.Qeff} when the function enter the region of a cofactor
# delimited by \code{window}. Default = FALSE.
# 
# @param window \code{Numeric} distance (cM) on the left and the right of a
# cofactor position where it is not included in the model. Default = 20.
# 
# @param n.cores \code{Numeric}. Specify here the number of cores you like to
# use. Default = 1.
# 
#   
# @return Return:
# 
# \item{CIM }{\code{Data.frame} of class \code{QTLprof} with five columns :
# 1) QTL marker or in between position names; 2) chromosomes;
# 3) interger position indicators on the chromosome;
# 4) positions in centi-Morgan; and 5) -log10(p-values)}
# 
# @author Vincent Garin
# 
# 
# @seealso \code{\link{mpp_SIM}}, \code{\link{QTL_select}}
# 
# 
# @examples
# 
# data(mppData)
# 
# SIM <- mpp_SIM(mppData = mppData)
# cofactors <- QTL_select(SIM)[, 1]
# 
# CIM <- MQE_CIM(mppData = mppData, Q.eff = "cr", cofactors = cofactors,
#                cof.Qeff = c("anc", "par", "biall"))
#
# plot(x = CIM)
#                                
# @export
# 

MQE_CIM <- function(mppData = NULL, trait = 1, Q.eff = "cr", VCOV = "h.err",
                    cofactors = NULL, cof.Qeff, chg.Qeff = FALSE, window = 20,
                    n.cores = 1){
  
  # 1. Check data format and arguments
  ####################################
  
  check.MQE(mppData = mppData, Q.eff = Q.eff, trait = trait,
            VCOV = VCOV, cofactors = cofactors, cof.Qeff = cof.Qeff,
            n.cores = n.cores, fct = "CIM")
  
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
  
  ### 2.6 Optional cluster
  
  if(n.cores > 1){
    
    parallel <- TRUE
    cluster <- makeCluster(n.cores)
    
  } else {
    
    parallel <- FALSE
    cluster <- NULL
    
  }
  
  
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
  
  if(n.cores > 1){stopCluster(cluster)}
  
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