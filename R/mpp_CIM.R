###########
# mpp_CIM #
###########

#' MPP Composite Interval Mapping
#' 
#' Compute QTL models along the genome using cofactors representing other
#' genetic positions for control.
#' 
#' For more details about the different models, see documentation of the
#' function \code{\link{mpp_SIM}}. The function returns a -log10(p-value) QTL
#' profile.
#' 
#' @param mppData An object of class \code{mppData}.
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. For more details see
#' \code{\link{mpp_SIM}}. Default = "cr".
#' 
#' @param cofactors Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function \code{\link{QTL_select}} or
#' vector of \code{character} marker positions names.
#' Default = NULL.
#' 
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 20.
#' 
#' @param plot.gen.eff \code{Logical} value. If \code{plot.gen.eff = TRUE},
#' the function will save the decomposed genetic effects per cross/parent.
#' These results can be plotted with the function \code{\link{plot.QTLprof}}
#' to visualize a genome-wide decomposition of the genetic effects.
#' \strong{This functionality is only available for the cross-specific,
#' parental and ancestral models.}
#' Default value = FALSE.
#' 
#' @param n.cores \code{Numeric}. Specify here the number of cores you like to
#' use. Default = 1.
#' 
#'   
#' @return Return:
#' 
#' \item{CIM }{\code{Data.frame} of class \code{QTLprof}. with five columns :
#' 1) QTL marker names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-val). And if
#' \code{plot.gen.eff = TRUE}, p-values of the cross or parental QTL effects.}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mpp_SIM}}, \code{\link{QTL_select}}
#' 
#' @examples
#' 
#' # Cross-specific effect model
#' #############################
#' 
#' data(mppData)
#' 
#' SIM <- mpp_SIM(mppData = mppData, Q.eff = "cr")
#' 
#' cofactors <- QTL_select(Qprof = SIM, threshold = 3, window = 20)
#' 
#' CIM <- mpp_CIM(mppData = mppData, Q.eff = "cr", cofactors = cofactors,
#'                window = 20, plot.gen.eff = TRUE)
#' 
#' plot(x = CIM)
#' plot(x = CIM, gen.eff = TRUE, mppData = mppData, Q.eff = "cr")
#' 
#' # Bi-allelic model
#' ##################
#' 
#' cofactors <- mppData$map[c(15, 63), 1]
#' 
#' CIM <- mpp_CIM(mppData = mppData, Q.eff = "biall", cofactors = cofactors,
#'                window = 20)
#' 
#' plot(x = CIM, type = "h")
#'                                
#' @export
#' 


mpp_CIM <- function(mppData, trait = 1, Q.eff = "cr", cofactors = NULL,
                    window = 20, plot.gen.eff = FALSE, n.cores = 1)
{
  
  # 1. Check data format and arguments
  ####################################
  
  check.model.comp(mppData = mppData, trait = trait, Q.eff = Q.eff,
                   VCOV = 'h.err', plot.gen.eff = plot.gen.eff,
                   n.cores = n.cores, cofactors = cofactors, fct = "CIM")
  
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
  
  ### 2.5 Optional cluster
  
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
  
  if(n.cores > 1){stopCluster(cluster)}
  
  log.pval <- t(data.frame(log.pval))
  if(plot.gen.eff){log.pval[is.na(log.pval)] <- 1}
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