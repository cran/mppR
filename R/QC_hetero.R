#############
# QC_hetero #
#############

# Marker heterozigosity rate
# 
# Compute marker heterozigosity rate.
# 
# @param mk.mat \code{Character} marker score \code{matrix} with genotypes as
# row and markers as column. \strong{Rows and columns names must be the genotype
# and marker identifiers respectively. Marker scores must be coded using one
# letter per allele. For example, AA, CC, GG, TT, AC, AG, AT, CA, CG, CT,
# GA, GC, GT, TA, TC, TG. Missing values must be coded \code{NA}.}
# 
# @param parallel Logical value specifying if the function should be executed
# in parallel on multiple cores. Default value = FALSE.
# 
# @param cluster Cluster object obtained with the function \code{makeCluster()}
# from the parallel package. Default = NULL.
# 
# @return Return:
# 
# \item{het.rate}{\code{Numeric} vector specifying the heterozygosity rate
# for each marker.}
# 
# 
# @author Vincent Garin
# 
# @examples
# 
# data('USNAM_geno')
# 
# hetero <- QC_hetero(mk.mat = USNAM_geno[7:dim(USNAM_geno)[1], ])
#
# @export
# 


QC_hetero <- function(mk.mat, parallel = FALSE, cluster = NULL) {
  
  # 1. Check data format
  ######################
  
  # 1. check data format
  ######################
  
  stopifnot(is.matrix(mk.mat), is.character(mk.mat))
  
  # 2. function to detect heterozygosity
  ######################################
  
  het.test <- function(x){
    
    if(is.na(x)) { NA 
      
      } else {
      
        alleles <- strsplit(x, "")
      
      vapply(X = alleles, FUN = function(x) (x[1] != x[2])*1, numeric(1))
      
    }
    
  }
  
  # 3. determine heterozygosity rate
  ##################################
    
  par.het <- apply(X = mk.mat, MARGIN = c(1,2), FUN = het.test)
  
  if(parallel){
    
    het.rate <- parApply(cl = cluster, X = par.het, MARGIN = 2,
                  FUN = function(x) sum(x, na.rm = TRUE) / length(x[!is.na(x)]))
    
  } else {
    
    
    het.rate <- apply(X = par.het, MARGIN = 2,
                  FUN = function(x) sum(x, na.rm = TRUE) / length(x[!is.na(x)])) 
    
  }
  
  return(het.rate)
    
}