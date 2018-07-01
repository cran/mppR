##############
# QTL_select #
##############

#' QTL candidates selection
#' 
#' Selection of QTL candidate positions.
#' 
#' The function select QTL positions that are above the given
#' \code{threshold} per chromosome. Once a position has been selected, and
#' exclusion \code{window} is set around that position. Positions falling into
#' that region will not be candidate anymore. The search continue until there
#' is no more candidate position above the \code{threshold}.
#' 
#' @param Qprof Object of class \code{QTLprof} returned by the function 
#' \code{\link{mpp_SIM}} or \code{\link{mpp_CIM}}. 
#' 
#' @param threshold \code{Numeric} value representing -log10(p-value) threshold
#' above which a position can be considered as a QTL candidate. Default = 3.
#' 
#' @param window \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected positions. Default = 50.
#' 
#' @param verbose \code{Logical} value specifying if the detection of no QTL
#' should be printed. Default = TRUE.
#' 
#' @return Return:
#' 
#' \item{QTL }{\code{Data.frame} of class \code{QTLlist} with five columns :
#' 1) QTL marker names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-values).}
#' 
#' @seealso \code{\link{mpp_SIM}}, \code{\link{mpp_CIM}}, \code{\link{mpp_perm}}
#' 
#' @references
#' 
#' This function is a modification of the QTL.reduce function
#' coming from the Biometris pipepline.
#' 
#' RAP (R Analytical Pipeline) (V0.9.1) May 2011
#' 
#' Authors: Paul Eilers (1), Gerrit Gort (1), Sabine Schnabel (1), Lucia
#' Gutierrez(1, 2), Marcos Malosetti(1), Joost van Heerwaarden, and Fred van
#' Eeuwijk(1)
#' 
#' (1) Wageningen University and Research Center, Netherlands (2) Facultad de
#' Agronomia, UDELAR, Uruguay
#' 
#' @examples
#' 
#' data(mppData)
#' 
#' SIM <- mpp_SIM(mppData)
#' 
#' QTL <- QTL_select(Qprof = SIM, threshold = 3)
#' 
#' @export
#' 


QTL_select <- function(Qprof, threshold = 3, window = 50, verbose = TRUE) {
  
  # 1. verify the format of the data
  ##################################
  
  stopifnot(inherits(Qprof, "QTLprof"))
  
  
  # 2. QTL selection procedure
  ############################
  
  pot.qtl <- Qprof[which(Qprof$log10pval > threshold), ]  
  res.qtl <- c()  
  
  if (nrow(pot.qtl) > 0) {
    
    pot.qtl$select <- 1  #select flag to be modified in loop
    pot.qtl$eval <- 0  #evaluated flag to be modified in loop
    
    for (Chr in unique(pot.qtl$chr)) {
      
      # Select the potential on the chromosome we are inspecting
      t.pot.qtl <- pot.qtl[pot.qtl$chr == Chr, ]
      
      while (sum(t.pot.qtl$eval) < nrow(t.pot.qtl)) {
        
        # Search the maximal score of the position for which the evaluation
        # flag is 0
        max.p <- max(t.pot.qtl$log10pval[which(t.pot.qtl$eval == 0)])
        # Selection of the row for which we have found the maximum score in the
        # previous test
        sel.row <- which(t.pot.qtl$log10pval == max.p & t.pot.qtl$eval == 0)[1]
        # We calculate the distance between the the highest peak and the all
        # the selected one
        d <- abs(t.pot.qtl$pos.cM - t.pot.qtl$pos.cM[sel.row])
        # Then we flag the potential qtl positions:
        
        t.pot.qtl$select[d <= window] <- 0  # the Marker positionned inside the
        # confidence range of the highest marker receive a 0 for selected
        t.pot.qtl$eval[d <= window] <- 1  # and a 1 for evaluation
        t.pot.qtl$select[sel.row] <- 1  # the marker which had the highest
        # position receive a 1 for selection
        
      } 
      
      t.pot.qtl <- t.pot.qtl[which(t.pot.qtl$select == 1), ]
      res.qtl <- rbind(res.qtl, t.pot.qtl)
    }
    
    res.qtl$select <- NULL
    res.qtl$eval <- NULL
    
  }
  
  QTL <- res.qtl
  
  if(is.null(QTL)){
    
    if(verbose){
      message("no position has been selected as QTL candidate")
    }
    
    return(QTL)
    
    
  } else {
    
    QTL <- QTL[, 1:5]
    
    class(QTL) <- c("QTLlist", "data.frame")
    
    return(QTL)
    
  }
  
}