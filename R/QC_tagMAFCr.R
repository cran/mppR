###############
# QC_tagMAFCr #
###############

# Tag problematic MAF
# 
# Detect markers with problematic minor allele frequency (MAF)
# within crosses. 
# 
# If \code{tag.mono = FALSE}, the within cross monomorphic markers are not
# considered as problematic.
# 
# @param MAF Object of class \code{mafRes} obtained with the function
# \code{\link{QC_MAF}} using a non null \code{cross.ind} argument.
# 
# @param MAF.lim \code{Numeric} vector with length equal to the number of
# crosses specifying the within cross MAF value bellow which a marker is tagged.
# Default = 0.05 for each cross.
# 
# @param tag.mono \code{Logical} value. If tag.mono = TRUE. Monomorphic
# markers within cross are also tagged as problematic. Default = FALSE.
# 
# 
# @param parallel \code{Logical} value specifying if the function should be
# executed in parallel on multiple cores. Default = FALSE.
# 
# @param cluster Cluster object obtained with the function \code{makeCluster()}
# from the \code{parallel} package. Default = NULL.   
# 
# @return Return:
#  
# \item{MAF.cr.ind}{\code{Logical} vector specifying if for at least one cross
# the MAF is bellow \code{MAF.lim} (without beeing monomorphic if
# \code{tag.mono = FALSE}).}
# 
# 
# @author Vincent Garin
# 
# @seealso
# 
# \code{\link{QC_MAF}}
# 
# @examples
# 
# data('USNAM_geno')
# 
# # obtain MAF object
# 
# cross.ind <- substr(rownames(USNAM_geno)[7:dim(USNAM_geno)[1]], 1, 4)
# 
# MAF <- QC_MAF(mk.mat = USNAM_geno[7:dim(USNAM_geno)[1], ],
#               cross.ind = cross.ind)
# 
# # Analyse MAF within crosses
# 
# MAF.cr.ind <- QC_tagMAFCr(MAF)
# 
# MAF$MAF.cr[, MAF.cr.ind] # Only a single marker has a MAF that is lower than
# # 0.05 but that is not monomorphic
# 
# MAF.cr.ind <- QC_tagMAFCr(MAF, tag.mono = TRUE)
# 
# sum(MAF.cr.ind) # 72 marker are tagged because they are monomorphic in at least
# # 1 cross, 1 because it MAF is bellow 0.05 in at least one cross without
# # being monomorphic in any cross.
#
# @export
# 


QC_tagMAFCr <- function(MAF, MAF.lim = NULL, tag.mono = FALSE, parallel = FALSE,
                        cluster = NULL){
  
  
  # 1. check that the object if of class mafRes
  #############################################
  
  if(!is(MAF, "mafRes")) {
    
    stop("'MAF' must be of class ", dQuote("mafRes"))
    
    if(length(MAF) !=2 ){
      
      stop("'MAF' does not contain the within cross MAF. ",
           "Use QC_MAF() with 'cross.ind' to get such an object")
      
    }
    
  }
  
  # 2. check the within cross MAF according to the given criteria
  ###############################################################
  
  MAF.cr <- MAF$MAF.cr
  
  if(is.null(MAF.lim)) ub <- rep(0.05, dim(MAF.cr)[1]) else ub <- MAF.lim
  
  # default value for the MAF
  # limit (0.05) if the user did not enter any value.
  
  
  if(!tag.mono){ # bellow MAF.lim but not monomorphic
    
    tag.cr <- function(x, ub) sum((x != 0) & (x < ub)) > 0  
    
  } else { # only belloe MAF.lim
    
    tag.cr <- function(x, ub) sum(x < ub) > 0  
    
  }
  
    
    if(parallel) {
      
      MAF.cr.ind <- parApply(cl = cluster, X = MAF.cr, MARGIN = 2,
                             FUN = tag.cr, ub = ub)
      
      return(MAF.cr.ind)
      
    } else {
      
      MAF.cr.ind <- apply(X = MAF.cr, MARGIN = 2, FUN = tag.cr, ub = ub)
      
      return(MAF.cr.ind)
      
    }
  
}