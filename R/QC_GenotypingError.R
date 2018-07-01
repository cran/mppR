######################
# QC_GenotypingError #
######################

# Detect genotyping errors
# 
# Check if markers have more than two alleles.
# 
# @param mk.mat \code{Character} marker scores \code{matrix} with genotypes as
# row and markers as column. \strong{Rows and columns names must be the genotype
# and marker identifiers respectively. Marker scores must be coded using one
# letter per allele. For example, AA, CC, GG, TT, AC, AG, AT, CA, CG, CT,
# GA, GC, GT, TA, TC, TG. Missing values must be coded \code{NA}.}
# 
# @param parallel \code{Logical} value specifying if the function should be
# executed in parallel on multiple cores. Default = FALSE.
# 
# @param cluster Cluster object obtained with the function \code{makeCluster()}
# from the parallel package. Default = NULL.
# 
# @return Return:
# 
# \item{prob.mk}{\code{Character} vector containing the name of marker with
# more than two alleles.}
# 
# @author Vincent Garin
# 
# @seealso
# 
# \code{\link{QC_MAF}}, \code{\link{QC_missing}}
# 
# @examples
# 
# data(USNAM_geno)
# 
# mk.mat <- USNAM_geno
# mk.mat[1, 10] <- "TC" 
# QC_GenotypingError(mk.mat)
#
# @export
# 
# 

QC_GenotypingError <- function(mk.mat, parallel = FALSE, cluster = NULL){
  
  # 1. check data format
  ######################
  
  stopifnot(is.matrix(mk.mat), is.character(mk.mat))
  
  
  # 2. define a single vector function
  ####################################

  test.vect <- function(x) {
    
    x <- x[!is.na(x)] # remove missing values
    
    if (length(x) == 0){ # Only missing values
      
      NA
      
    } else {
      
      # split all marker score into single allele scores
      
      alleles <- c(vapply(X = x, FUN = function(x) unlist(strsplit(x, split = "")),
                          FUN.VALUE = character(2)))
      
      if(length(table(alleles)) > 2){ TRUE } else { FALSE } 
      
    } 
    
  }
  
  # 3. apply the function on the marker matrix
  ############################################
  
  if(parallel){
    
    test <- parApply(cl = cluster, X = mk.mat, MARGIN = 2, FUN = test.vect)
    
  } else {
   
    test <- apply(X = mk.mat, MARGIN = 2, FUN = test.vect) 
    
  }
  
  
  
  # 3. outputs
  ############
  
  if(sum(test, na.rm = TRUE) > 0){
    
  prob.mk <- colnames(mk.mat)[test]
  return(prob.mk)
  
    
  } else {
    
    prob.mk <- NULL
    return(prob.mk)
  
    }
  
  
}