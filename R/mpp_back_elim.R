################
# mpp_back_elim #
################

#' Backward elimination on QTL candidates
#' 
#' Performs a backward elimination using a list of given QTLs positions. The
#' positions with a p-value above the significance level \code{alpha}, are
#' successively removed.
#' 
#' The function starts with all QTL positions in the model and test the inclusion
#' of each position as the last in the model. If all position p-values are below
#' \code{alpha} the procedure stop. If not the position with the highest p-value
#' is remove and the procedure continue until there is no more unsignificant
#' position.
#' 
#' @param mppData An object of class \code{mppData}.
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function \code{\link{QTL_select}} or
#' vector of \code{character} marker positions names.
#' Default = NULL.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. For more details see
#' \code{\link{mpp_SIM}}. Default = "cr".
#' 
#' @param alpha \code{Numeric} value indicating the level of significance for
#' the backward elimination. Default = 0.05.
#' 
#' @return Return:
#' 
#' \item{QTL }{\code{Data.frame} of class \code{QTLlist} with five columns :
#' 1) QTL marker names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-values).}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mpp_SIM}}
#' 
#' @examples
#' 
#' data(mppData)
#' 
#' SIM <- mpp_SIM(mppData)
#' 
#' QTL <- QTL_select(SIM)
#' 
#' QTL.sel <- mpp_back_elim(mppData = mppData, QTL = QTL)
#' 
#' @export
#' 


mpp_back_elim <- function (mppData, trait = 1, QTL = NULL, Q.eff = "cr",
                         alpha = 0.05) {
  
  # 1. Check data format
  ######################
  
  check.model.comp(mppData = mppData, trait = trait, Q.eff = Q.eff,
                   VCOV = 'h.err', QTL = QTL, fct = "back")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  
  ### 2.3 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.4 Formation of the list of QTL
  
  if(is.character(QTL)){
    
    Q.pos <- which(mppData$map[, 1] %in% QTL)
    
    QTL <- mppData$map[mppData$map[, 1] %in% QTL, ]
    
  } else {
    
    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])
    
  }
  
  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                    Q.eff = Q.eff, order.MAF = TRUE)
  
  names(Q.list) <- paste0("Q", 1:length(Q.list))
  
  
  # 3. Compute the models
  #######################
  
  ind <- TRUE
  
  while(ind) {
    
    ### 3.1 elaboration of model formulas
    
    model.formulas <- formula_backward(Q.names = names(Q.list), VCOV = 'h.err')
    
    ### 3.2 computation of the models
    
    pvals <- lapply(X = model.formulas, FUN = QTLModelBack, mppData = mppData,
                    trait = t_val, Q.list = Q.list, cross.mat = cross.mat,
                    VCOV = 'h.err')
    
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
  
  QTL <- QTL[as.numeric(substr(names(Q.list), 2, nchar(names(Q.list)))), ,
             drop = FALSE]
  
  if(dim(QTL)[1] == 0){
    
    QTL <- NULL
    
    return(QTL)
    
  } else{
    
    class(QTL) <- c("QTLlist", "data.frame")
    
    return(QTL)
    
  }
  
}