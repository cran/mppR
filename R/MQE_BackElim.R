################
# MQE_BackElim #
################

# Backward elimination on multi-QTL effect candidates
# 
# Performs a backward elimination using a list of given QTLs positions. These
# position can have different type of QTL effects (cross-specific, parental,
# ancestral or bi-allelic).
# 
# The function starts with all QTL positions in the model and test the inclusion
# of each position as the last in the model. If all position p-values are below
# \code{alpha} the procedure stop. If not the position with the highest p-value
# is remove and the procedure continue until there is no more unsignificant
# position.
# 
# \strong{WARNING!} The computation of random pedigree models
# (\code{VCOV = "pedigree" and "ped_cr.err"}) can sometimes fail. This could be
# due to singularities due to a strong correlation between the QTL term(s) and 
# the polygenic term. This situation can appear in the parental model.
# the error can also sometimes come from the \code{asreml()} function. From
# our experience, in that case, trying to re-run the function one or two times
# allow to obtain a result.
# 
# @param mppData An object of class \code{mppData}
#
# @param trait \code{Numerical} or \code{character} indicator to specify which
# trait of the \code{mppData} object should be used. Default = 1.
# 
# @param QTL Vector of \code{character} markers or in between marker positions
# names. Default = NULL.
# 
# @param Q.eff \code{Character} vector indicating for each QTL position the
# type of QTL effect among: "cr", "par", "anc" and "biall". For details look
# at \code{\link{mpp_SIM}}.
#
# @param VCOV \code{Character} expression defining the type of variance
# covariance structure used: 1) "h.err" for an homogeneous variance residual term
# (HRT) linear model; 2) "h.err.as" for a HRT model fitted by REML using
# \code{ASReml-R}; 3) "cr.err" for a cross-specific variance residual terms
# (CSRT) model; 4) "pedigree" for a random pedigree term and HRT model;
# and 5) "ped_cr.err" for random pedigree and CSRT model.
# For more details see \code{\link{mpp_SIM}}. Default = "h.err".
# 
# @param alpha \code{Numeric} value indicating the level of significance for
# the backward elimination. Default = 0.05.
# 
# @return Return:
# 
# \item{QTL }{\code{Data.frame} of class \code{QTLlist} with six columns :
# 1) QTL marker or in between position names; 2) chromosomes;
# 3) interger position indicators on the chromosome;
# 4) positions in centi-Morgan; 5) -log10(p-values), and 6) type of QTL effect.}
#  
# @author Vincent Garin
# 
# @examples
# 
# data(mppData)
# 
# SIM <- mpp_SIM(mppData = mppData)
# QTL <- QTL_select(SIM)
# 
# QTL <-  MQE_BackElim(mppData = mppData, QTL = QTL[, 1],
#                      Q.eff = c("anc", "par", "biall"))
# 
# 
# @export
# 


MQE_BackElim <- function (mppData = NULL, trait = 1, QTL = NULL, Q.eff,
                          VCOV = "h.err", alpha = 0.05) {
  
  # 1. Check data format
  ######################
  
  check.MQE(mppData = mppData, trait = trait, Q.eff = Q.eff,
            VCOV = VCOV, QTL = QTL, fct = "QTLeffects")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  ### 2.3 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.4 Formation of the list of QTL
  
  # order list of QTL positions
  
  Q.pos <- vapply(X = QTL,
                  FUN = function(x, mppData) which(mppData$map[,1] == x),
                  FUN.VALUE = numeric(1), mppData = mppData)
  
  Q.ord <- data.frame(QTL, Q.eff, Q.pos, stringsAsFactors = FALSE)
  
  Q.ord <- Q.ord[order(Q.pos),]
  
  QTL <- Q.ord[, 1]; Q.eff <- Q.ord[, 2]
  
  Q.pos <- which(mppData$map[, 1] %in% QTL)
  
  Q.list <- mapply(FUN = inc_mat_QTL, x = Q.pos, Q.eff = Q.eff,
                   MoreArgs = list(mppData = mppData, order.MAF = TRUE),
                   SIMPLIFY = FALSE)
  
  names(Q.list) <- paste0("Q", 1:length(Q.list))
  
  
  
  # 3. Compute the models
  #######################
  
  ind <- TRUE
  
  while(ind) {
    
    ### 3.1 elaboration of model formulas
    
    model.formulas <- formula_backward(Q.names = names(Q.list), VCOV = VCOV)
    
    ### 3.2 computation of the models
    
    pvals <- lapply(X = model.formulas, FUN = QTLModelBack, mppData = mppData,
                    trait = t_val, Q.list = Q.list, cross.mat = cross.mat,
                    VCOV = VCOV)
    
    
    pvals <- unlist(pvals)
    
    
    ### 3.4 test the p-values
    
    
    if(sum(pvals > alpha) > 0) {
      
      # remove the QTL position from the list
      
      Q.list <- Q.list[!(pvals==max(pvals))]
      
      # test if there is no more positions
      
      if(length(Q.list) == 0){ind <- FALSE}
      
    } else {
      
      # stop the procedure
      
      ind <- FALSE
      
    }
    
  }
  
  # 4. return the final list of QTLs
  ##################################
  
  index.rem <- as.numeric(substr(names(Q.list), 2, nchar(names(Q.list))))
  
  QTL.rem <- QTL[index.rem] ; Qeff.rem <- Q.eff[index.rem]
  
  QTL <- mppData$map[mppData$map[, 1] %in% QTL.rem, ]
  
  QTL <- data.frame(QTL, Qeff.rem, stringsAsFactors = FALSE)
  
  colnames(QTL)[5] <- "QTL.eff"
  
  if(dim(QTL)[1] == 0){
    
    QTL <- NULL
    
  }
  
  return(QTL)
  
}