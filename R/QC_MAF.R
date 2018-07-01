##########
# QC_MAF #
##########

# Minor allele frequency
# 
# Compute minor allele frequency (MAF) at the whole population and
# within crosses if a cross indicator vector is provided to \code{cross.ind}.
# 
# @param mk.mat \code{Character} marker score \code{matrix} with genotypes as
# row and markers as column. \strong{Rows and columns names must be the genotype
# and marker identifiers respectively. Marker scores must be coded using one
# letter per allele. For example, AA, CC, GG, TT, AC, AG, AT, CA, CG, CT,
# GA, GC, GT, TA, TC, TG. Missing values must be coded \code{NA}.}
# 
# @param cross.ind \code{Character} vector indicating to which cross each
# genotype belongs. If \code{cross.ind = NULL}, the function do not consider cross
# subdivisions and only calculate MAF at the population level. Default = NULL.
# 
# @param parallel \code{Logical} value specifying if the function should be
# executed in parallel on multiple cores. Default = FALSE.
# 
# @param cluster Cluster object obtained with the function \code{makeCluster()}
# from the parallel package. Default = NULL.   
# 
# @return Return:
# 
# Object of class mafRes MAF containing the following element(s).
# If cross.ind = NULL
# 
# \item{MAF.pop}{Vector of marker allele MAF at the
# population level. \code{NA} if all markers scores are missing.}
# 
# If cross.ind is not NULL
# 
# \item{MAF.pop}{Vector of marker allele MAF at the
# population level. \code{NA} if all markers scores are missing.}
# 
# \item{MAF.cr}{\code{Matrix} of marker allele MAF within crosses.
# \code{NA} if all markers scores are missing.}
# 
# @author Vincent Garin
# 
# @seealso
# 
# \code{\link{QC_tagMAFCr}}
# 
# @examples
# 
# data(USNAM_geno)
# cross.ind <- substr(rownames(USNAM_geno)[7:dim(USNAM_geno)[1]], 1, 4)
# 
# MAF <- QC_MAF(mk.mat = USNAM_geno[7:dim(USNAM_geno)[1], ],
#              cross.ind = cross.ind)
#
# @export
# 


QC_MAF <- function(mk.mat, cross.ind = NULL, parallel = FALSE,
                   cluster = NULL) {
  
  # 1. check data format
  ######################
  
  stopifnot(is.matrix(mk.mat), is.character(mk.mat))
  
  
  # 2. function to calculate MAF on a single vector
  #################################################
  
  MAF.vect <- function(x) {
    
    x <- x[!is.na(x)] # remove missing values
    
    if (length(x) == 0){ # Only missing values
      
      return(NA)
      
    } else {
      
      # split all marker score into single allele scores
      
      alleles <- c(vapply(X = x, FUN = function(x) unlist(strsplit(x, split = "")),
                          FUN.VALUE = character(2)))
      
      if(length(table(alleles)) == 1){ # monomorphic marker
        
        return(0)
        
      } else {
        
        return(min(table(alleles) / length(alleles)))
        
      } 
      
    } 
    
  } # end single vector MAF function
  
  if(is.null(cross.ind)){ # compute MAF only at population level
    
    ### 3.1 MAF at population level
    
    if(parallel) {
      
      MAF.pop <- parApply(cl = cluster, X = mk.mat, MARGIN = 2, FUN = MAF.vect)  
      
    } else {
      
      MAF.pop <- apply(X = mk.mat, MARGIN = 2, FUN = MAF.vect)  
      
    }
    
    MAF <- MAF.pop
    class(MAF) <- c("mafRes")
    
    return(MAF)
    
    
  } else { # compute MAF at population and cross level
    
    ### 3.1 check that the cross indicator vector has the same length
    
    if(length(cross.ind) != dim(mk.mat)[1]){
      
      stop("'cross.ind' length is not the same as the number of genotype in 'mk.mat'")
    }
    
    ### 3.2 MAF at population level
    
    if(parallel) {
      
      MAF.pop <- parApply(cl = cluster, X = mk.mat, MARGIN = 2, FUN = MAF.vect)  
      
    } else {
      
      MAF.pop <- apply(X = mk.mat, MARGIN = 2, FUN = MAF.vect)  
      
    }
    
    ### 3.3 MAF at the cross level
    
    # function to apply MAF computation on within cross
    
    cross.ind <- factor(cross.ind, levels = unique(cross.ind))
    
    MAF.vect.cr <- function(x) tapply(X = x, INDEX = cross.ind, FUN = MAF.vect)
    
    if(parallel) {
      
      MAF.cr <- parApply(cl = cluster, X = mk.mat, MARGIN = 2, FUN = MAF.vect.cr)
      
    } else {
      
      MAF.cr <- apply(X = mk.mat, MARGIN = 2, FUN = MAF.vect.cr)  
      
    }
    
    MAF <-  list(MAF.pop = MAF.pop, MAF.cr = MAF.cr)
    class(MAF) <- c("list", "mafRes")
    
    return(MAF)
    
    }
  
}