#################
# IncMat_parent #
#################

# Skeleton of parental incidence matrices
# 
# Form the skeleton matrices used to form parental incidence
# matrices in parental effect models.
# 
# @param mppData An object of class \code{mppData}.
# 
# @return Return:
# 
# \code{List} with two matrices
# 
# \item{PA}{Incidence matrix specifying for each genotype which parent
# is the parent 1 or A}
# 
# \item{PB}{Incidence matrix specifying for each genotype which parent
# is the parent 2 or B}
# 
# @author Vincent Garin
# 
# @examples
# 
# data(mppData)
# 
# par.mat <- IncMat_parent(mppData)
# 
# @export
#    


IncMat_parent <- function(mppData){
  
  stopifnot(inherits(mppData, "mppData"))
  
  # cross incidence matrix formation
  
  cross.matrix <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  # list of parents present in at least one of the crosses
  
  par.per.cross <- mppData$par.per.cross
  
  parents <- mppData$parents
  
  # Parent 1
  
  PA <- c()
  
  for(i in 1:mppData$n.cr){
    
    # which parent is the parent of the cross
    
    PA <- rbind(PA, (par.per.cross[i, 2] == parents)*1)
    
  }
  
  # Parent 2
  
  PB <- c()
  
  for(i in 1:mppData$n.cr){
    
    # which parent is the parent of the cross
    
    PB <- rbind(PB,(par.per.cross[i, 3] == parents)*1)
    
  }
  
  # multiplication of the cross matrix by the parent matrix
  
  PA <- cross.matrix %*% PA
  
  PB <- cross.matrix %*% PB
  
  colnames(PA) <- colnames(PB) <- parents
  
  return(list(PA = PA, PB = PB))
  
  
}