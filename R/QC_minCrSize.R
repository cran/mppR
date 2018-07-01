################
# QC_minCrSize #
################

# Minimum cross size
# 
# Remove the crosses that do not have the minimum critical size.
# 
# @param mk.mat Marker score \code{matrix} or \code{data.frame} with genotype as
# row and markers as column.
# 
# @param cross.ind \code{Character} vector indicating to which cross each
# genotype belong.
# 
# @param n.lim \code{Numeric} value specifying the minimum cross size.
# Default = 15.   
# 
# @return Return:
# 
# \item{mk.mat.red}{Reduced marker \code{matrix} with too small crosses removed.}
# 
# @author Vincent Garin
# 
# @examples
# 
# data('USNAM_geno')
# 
# USNAM_geno.red <- USNAM_geno[-c(7:101), ] # remove 95 genotypes in the first cross
# 
# cross.ind <- substr(rownames(USNAM_geno.red), 1, 4)
# 
# mk.mat.red <- QC_minCrSize(mk.mat = USNAM_geno.red, cross.ind = cross.ind)
# 
# @export
# 



QC_minCrSize <- function(mk.mat, cross.ind, n.lim = 15) {
  
  ### 1. Check if cross.ind has the same length as the number of genotypes
  
  if(length(cross.ind) != dim(mk.mat)[1]){
    
    stop("'cross.ind' length is not equal to the number of genotypes in 'mk.mat'")
  }
  
  ### 2. remove the too small crosse
  
  n.cr <- table(cross.ind)
  
  cr.ind <- cross.ind %in% attr(n.cr, "dimnames")[[1]][n.cr >= n.lim]
  
  mk.mat.red <- subset(x = mk.mat, subset = cr.ind, drop = FALSE)
  
  return(mk.mat.red)
  
}
