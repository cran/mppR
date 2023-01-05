##############
# inc_mat_QTL #
##############

#' QTL incidence matrix
#'
#' Build a single position QTL incidences matrix.
#' 
#' @param x \code{Integer} value indicating the genetic position on the map
#' (\code{mppData$map}) of the QTL incidence matrix.
#' 
#' @param mppData An object of class \code{mppData}.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. For more details see
#' \code{\link{mpp_SIM}}. Default = "cr".
#' 
#' @param order.MAF \code{Logical} value specifying if the QTL incidence matrix
#' should be ordered by allele frequency for a parental and ancestral QTL
#' incidence matrix. The colum will be ordred from the least to the most frequent
#' allele. Default = FALSE.
#' 
#' @param ref_par Optional \code{Character} expression defining the parental
#' allele that will be used as reference for the parental model. Default = NULL
#' 
#' @return Return:
#' 
#' \item{QTL.mat}{QTL incidence matrix. For the
#' cross-specific model, it represents the difference between the
#' number of allele from parent 2 or B and parent 1 or A divided by two. For
#' parental (ancestral) model it represents the expected number of parental
#' (ancestral) allele copies. For the bi-allelic model, it represents the number
#' of copies of the least frequent allele.}
#' 
#' @author Vincent Garin
#'
#' @examples
#' 
#' data(mppData)
#' 
#' QTLmatCr <- inc_mat_QTL(x = 2, mppData = mppData, Q.eff = "cr")
#' 
#' QTLmatPar <- inc_mat_QTL(x = 2, mppData = mppData, Q.eff = "par")
#' 
#' QTLmatAnc <- inc_mat_QTL(x = 2, mppData = mppData, Q.eff = "anc")
#' 
#' QTLmatBi <- inc_mat_QTL(x = 2, mppData = mppData, Q.eff = "biall")
#' 
#' 
#' @export
#' 


inc_mat_QTL <- function(x, mppData, Q.eff, order.MAF = FALSE, ref_par = NULL) {
  
  pos <- unlist(mppData$map[x, c(2,3)])
  
  if(Q.eff == "cr"){
    
    cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
    
    alpha.pred <- mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], mppData$n.zigo] -
      mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 1]
    
    alpha.pred <- t(rep(1, dim(cross.mat)[2])) %x% alpha.pred
    QTL.mat <- cross.mat * alpha.pred
    
  } else if(Q.eff == "par") {
    
    par.mat <- IncMat_parent(mppData = mppData)
    
    # get number of alleles from parent A and B
    
    if (mppData$n.zigo == 3){
      
      alleleA <- (2*mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 1]) +
        mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 2]
      
      alleleB <- (2*mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 3]) +
        mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 2]
      
    } else if (mppData$n.zigo == 2){
      
      alleleA <- (2*mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 1])
      alleleB <- (2*mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 2])
      
    }
    
    # form parental QTL matrix
    
    PA_pos <- (t(rep(1, dim(par.mat$PA)[2])) %x% alleleA) * par.mat$PA
    PB_pos <- (t(rep(1, dim(par.mat$PA)[2])) %x% alleleB) * par.mat$PB
    
    QTL.mat <- PA_pos + PB_pos
    
    
  } else if (Q.eff == "anc") {
    
    par.mat <- IncMat_parent(mppData = mppData)
    
    # form a parental matrix. Same as for Q.eff == "par"
    
    if (mppData$n.zigo == 3){
      
      alleleA <- (2*mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 1]) +
        mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 2]
      
      alleleB <- (2*mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 3]) +
        mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 2]
      
    } else if (mppData$n.zigo == 2){
      
      alleleA <- (2*mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 1])
      alleleB <- (2*mppData$geno.IBD$geno[[pos[1]]]$prob[, pos[2], 2])
      
    }
    
    PA_pos <- (t(rep(1, dim(par.mat$PA)[2])) %x% alleleA) * par.mat$PA
    PB_pos <- (t(rep(1, dim(par.mat$PA)[2])) %x% alleleB) * par.mat$PB
    
    QTL.mat <- PA_pos + PB_pos
    
    # modify parental matrix according to ancestral matrix (par. clustering)
    
    A.allele <- as.factor(mppData$par.clu[x, ])
    
    A <- model.matrix(~ A.allele - 1)
    
    QTL.mat <- QTL.mat %*% A
    
    
  } else if (Q.eff == "biall"){
    
    QTL.mat <- subset(x = mppData$geno.IBS, select = x, drop = FALSE)
    
  }
  
  if(order.MAF){
    
    # this part can be later extended to each interconnected part.
    
    if((Q.eff == "par") || (Q.eff == "anc")){
      
      # determine the most frequent allele
      
      all.freq <- apply(X = QTL.mat, MARGIN = 2,
                        FUN = function(x) sum(round(x, 3) != 0))
      
      QTL.mat <- QTL.mat[, names(sort(all.freq))]
      
    }
    
  }
  
  if((Q.eff == "par") & !is.null(ref_par)){
    
    name.Qmat <- colnames(QTL.mat)
    name.Qmat <- name.Qmat[name.Qmat != ref_par]
    QTL.mat <- QTL.mat[, c(name.Qmat, ref_par)]
    
  }
  
  return(QTL.mat)
  
}