##########
# QTL_CI #
##########

# -log10(p-value) drop-off confidence interval
# 
# Computes a -log10(p-value) drop-off confidence interval for QTL positions.
# 
# @param QTL Object of class \code{QTLlist} representing a list of
# selected position obtained with the function \code{\link{QTL_select}} or
# vector of \code{character} marker or in between marker positions names.
# Default = NULL.
# 
# @param Qprof Object of class \code{QTLprof} returned by the function
# \code{\link{mpp_SIM}} or \code{\link{mpp_CIM}}. The CIM profile that is
# used should be calculated excluding cofactors from the tested chromosome.
# Such a profile can be obtained using a \code{window} larger than the longer
# chromosome.
# 
# @param drop \code{numeric} -log10(p-value) drop value at the limits of the
# interval. Default = 1.5.
# 
# @return Return:
# 
# \item{CI }{\code{data.frame} with the following columns: 1) QTL marker or
# in between position names; 2) chromosomes; 3) QTL positions in cM;
# 4) inferior marker or in between position names; 5) inferior positions in cM;
# 6) superior marker or in between position names; 7) superior positions in cM;
# 8) ranges in cM.}
# 
# @author Vincent Garin
# 
# @seealso \code{\link{mpp_SIM}}, \code{\link{mpp_CIM}},
# \code{\link{QTL_select}}
# 
# @examples
# 
# data(mppData)
# 
# SIM <- mpp_SIM(mppData = mppData)
# 
# QTL <- QTL_select(Qprof = SIM, threshold = 3, window = 20)
# 
# QTL_CI(QTL = QTL, Qprof = SIM, drop = 1.5)
# 
# @export
#


QTL_CI <- function(QTL = NULL, Qprof, drop = 1.5) {
  
  # test data format
  ##################
  
  stopifnot(inherits(Qprof, "QTLprof"))
  
  if(is.null(QTL)){
    
    stop("no QTL position has been specified in 'QTL'")
    
  }
  
  
  
  ### end test data format
  
  Qprof <- Qprof[, 1:5]
  
  if(is.character(QTL)){
    
    QTL <- Qprof[Qprof[, 1] %in% QTL, ]
    
  } else {
    
    QTL <- Qprof[Qprof[, 1] %in% QTL[, 1], ]
    
  }

  # function to determine the two limits
  ######################################
  
  lim.fct <- function(x, QTL, Qprof, drop){
    
    Qi <- QTL[x, ]
    Qprofi <- Qprof[Qprof[, 2] == Qi[, 2], ]
    
    ### inferior position
    
    Qprofi.inf <- Qprofi[Qprofi[, 4] <= Qi[, 4], ]
    inf.diff <- Qi[, 5] - Qprofi.inf[, 5] # -log10(pval) difference
    
    # three situations:
    
    # A) position is the first of the chromosome
    
    if(length(inf.diff) == 1){ inf.pos <- Qi[, ] } else {
      
      # B) position is not a border position but l drop is not enough
      
      if(sum((inf.diff > drop)) == 0) { # all position are above the drop difference
        
        inf.pos <- Qprofi.inf[1, ]
        
      } else {
        
        # C) positin is not a border positin and l drop is sufficient
        
        inf.pos <- Qprofi.inf[max(which(inf.diff >= drop)), 1:5, ]
        
      }
      
    } 
    
    ### superior position
    
    Qprofi.sup <- Qprofi[Qprofi[, 4] >= Qi[, 4], ]
    sup.diff <- Qi[, 5] - Qprofi.sup[, 5] # -log10(pval) difference
    
    # three situations:
    
    # A) position is the last of the chromosome
    
    if(length(sup.diff) == 1){ sup.pos <- Qi[, ] } else {
      
      # B) position is not a border position but l drop is not enough
      
      if(sum((sup.diff > drop)) == 0) { # all position are above the drop difference
        
        sup.pos <- Qprofi.sup[dim(Qprofi.sup)[1], ]
        
      } else {
        
        # C) positin is not a border positin and l drop is sufficient
        
        sup.pos <- Qprofi.sup[min(which(sup.diff >= drop)), 1:5, ]
        
      }
      
    }
    
    range <- sup.pos[, 4] - inf.pos[, 4] 
    res <- data.frame(Qi[, -c(3, 5)], inf.pos[, -c(2, 3, 5)],
                      sup.pos[, -c(2, 3, 5)], range,
                      stringsAsFactors = FALSE)
    
    colnames(res) <- c("QTL.mk", "chr", "QTL.pos[cM]", "inf.mk",
                       "inf.lim[cM]", "sup.mk", "sup.lim[cM]", "range[cM]")
    
    return(res)
    
  }
  
  n.QTL <- dim(QTL)[1]
  
  res <- lapply(X = 1:n.QTL, FUN = lim.fct, QTL = QTL,
         Qprof = Qprof, drop = drop)
  
  CI <- c()
  for(i in seq_along(res)){ CI <- rbind.data.frame(CI, res[[i]]) }
  
  return(CI)
  
}
