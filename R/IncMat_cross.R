################
# IncMat_cross #
################

# Cross effect incidence matrix
# 
# Form a cross effect incidence matrix.
# 
# 
# @param cross.ind \code{Character} or \code{factor} vector specifying to
# which cross the genotypes belong.
# 
# @return Return:
# 
# \item{cross.mat}{Cross effect incidence matrix composed of 0 and 1 value
# indicating to which cross each genotype belongs.}
# 
# @author Vincent Garin
# 
# @examples
#
# data(mppData)
# 
# cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
#
# @export
#


IncMat_cross <- function(cross.ind){
  
  if(length(unique(cross.ind)) == 1){
    
    cross.mat <- matrix(1, length(cross.ind), 1)
    colnames(cross.mat) <- paste0("Cr", cross.ind[1])
    
  } else {
    
    Cr <- factor(x = cross.ind, levels = unique(cross.ind))
    cross.mat <- model.matrix( ~ Cr - 1, data = Cr)
    
  }
  
  return(cross.mat)
  
}