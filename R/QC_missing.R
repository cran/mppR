##############
# QC_missing #
##############

# Missing value for marker or genotype
# 
# Identify markers or genotypes with a missing rate higher than the
# specified \code{threshold}.
# 
# @param mk.mat \code{Character} marker score \code{matrix} with genotype as
# row and marker as column. \strong{Rows and columns names must be the genotype
# and marker identifiers respectively. Missing values must be coded \code{NA}.}
# 
# @param threshold \code{Numeric} value representing the missing values rate
# above which a marker or a genotype is considered as problematic.
# Default = 0.1.
# 
# @param MARGIN \code{Numeric} value giving the subscript which the function
# will be applied over. If MARGIN = 1,  the function looks at the genotypes
# (rows). If MARGIN = 2, the function looks at the markers (columns).
# Default = 2.
# 
# @return Return:
# 
# \item{mis.ind }{Three columns \code{data.frame} with: 1) marker or
# genotype names; 2) column or row positions in the matrix; and
# 3) The percentages of missing values.}
# 
# @author Vincent Garin
# 
# @examples
# 
# data('USNAM_geno')
# 
# # missingness of markers
# QC_missing(mk.mat = USNAM_geno)
# 
# # missingness of genotypes
# QC_missing(mk.mat = USNAM_geno, threshold = 0.2, MARGIN = 1)
# 
# @export
# 

QC_missing <- function(mk.mat, threshold = 0.1, MARGIN = 2) {
  
  # 1. Check data format
  ######################
  
  ### 1.1 check if it is a matrix or a data frame
  
  if (!(is.matrix(mk.mat))){
    
    stop("'mk.mat' is not a matrix")  
    
  }
  
  ### 1.2 check all columns are character
  
  if (!(is.character(mk.mat))){
    
    stop("'mk.mat' is not a character matrix")  
    
  }
  
  # 2. Compute missing rate
  #########################
  
  # function computing missing value rate
  
  missing.freq <- function(x) { sum(is.na(x))/length(x) }
  
  m.rate <- apply(X = mk.mat, MARGIN = MARGIN, FUN = missing.freq)
  
  ind <- m.rate > threshold # test marker or genotype problematic missing rate
  
  # 3. Format the results
  #######################
  
  # get the name of the marker/genotype with missing problem
  
  if (MARGIN == 1) { names <- rownames(mk.mat)[ind]
  col.names <- c("genotype", "position", "% missing")
    
  } else if (MARGIN == 2) { names <- colnames(mk.mat)[ind]
  col.names <- c("marker", "position", "% missing")}
  
  
  mis.ind <- data.frame(names, which(ind), m.rate[ind])
  colnames(mis.ind) <- col.names
 
  return(mis.ind)
  
}