############
# geno_012 #
############

# Marker scores into 0, 1, 2 format
#
# Transform genotype marker score into 0, 1, 2 format. The score represents the
# number of copies of the allele with the lowest marker allele frequency (MAF).
# 
# @param mk.mat \code{Character} genotype  marker score \code{matrix}.
# \strong{Marker scores must be coded using one letter for each allele.
# For example, AA, CC, GG, TT, AC, AG, AT, CA, CG, CT, GA, GC, GT, TA, TC, TG.
# Missing values must be coded NA.}
# 
# @param parallel \code{Logical} value specifying if the function should be
# executed in parallel on multiple cores. To run function in parallel user must
# provide cluster in the \code{cluster} argument. Default = FALSE.
# 
# @param cluster Cluster object obtained with the function \code{makeCluster()}
# from the parallel package. Default = NULL.
#   
# @return Return:
# 
# \code{List} with two matrices
# 
# \item{geno012}{Marker \code{matrix} with marker scores coded as 0, 1, 2
# corresponding to the number of copies of the least frequent SNP allele.}
# 
# \item{all.ref}{\code{matrix} with reference allele scores. The first row
# represents the minor allele (lowest frequency), the second the one represent the
# major allele (largest frequency) and the two others the heterozygous scores.} 
# 
# @author Vincent Garin    
# 
# @examples
# 
# data(USNAM_geno)
# 
# data <- geno_012(mk.mat = USNAM_geno)
# 
# data_012 <- data[[1]]
# 
# @export
#           


geno_012 <- function(mk.mat, parallel = FALSE, cluster = NULL) {
  
  geno.names <- rownames(mk.mat)
  
  # 1. Establish the allele with the highest MAF
  ##############################################
  
  ### 1.1 function to apply along the matrix
  
  # The function return vector of four characters. The first one represent
  # the allele with the highest MAF. Then the second one and the two
  # heterozygous
  
  
  allele.MAF <- function(x) {
    
    x <- x[!is.na(x)] # remove missing values
    
    if (length(x) == 0){ # Only missing values
      
      return(c(NA, NA, NA, NA))
      
    } else {
      
      # split all marker score into single allele scores
      
      alleles <- c(vapply(X = x, FUN = function(x) unlist(strsplit(x, split = "")),
                          FUN.VALUE = character(2)))
      
      alleles.sum <- table(alleles)
      
      if (length(alleles.sum) == 1) { # monomorphic case
        
        return(c(NA, NA, NA, NA))
        
      } else {
        
        all.sc <- unlist(attr(sort(alleles.sum, decreasing = FALSE), "dimnames"))
        all.sc <- c(paste0(all.sc[1], all.sc[1]), paste0(all.sc[2], all.sc[2]),
                    paste0(all.sc[1], all.sc[2]), paste0(all.sc[2], all.sc[1]))
        
        return(all.sc)
        
      }
      
    } 
    
  }
  
  ### 1.2 computation of the allele with the highest MAF
  
  if(parallel){
    
    all.class <- parApply(cl = cluster, X = mk.mat, MARGIN = 2,
                          FUN = allele.MAF)
    
  } else {
    
    all.class <- apply(X = mk.mat, MARGIN = 2, FUN = allele.MAF)
    
  }
  
  
  
  
  # 2. transform the marker scores. Given the allele with the highest MAF
  #######################################################################
  
  trans012 <- function(x, ref) {
    
    if(is.na(ref[1])){
      
      x <- rep(NA, length(x))
      
    } else {
      
      x[is.na(x)] <- NA
      x[x == ref[1]] <- 2
      x[x == ref[2]] <- 0
      x[x == ref[3]] <- 1
      x[x == ref[4]] <- 1
      
    }
    
    return(as.numeric(x))
    
  }
  
  mk.mat.list <- data.frame(mk.mat, stringsAsFactors = FALSE)
  all.class.list <- data.frame(all.class, stringsAsFactors = FALSE)
  
  if(parallel){
    
    geno012 <- clusterMap(cl = cluster, fun = trans012, x = mk.mat.list,
                          ref = all.class.list, SIMPLIFY = TRUE)
    
  } else {
    
    geno012 <- mapply(FUN = trans012, x = mk.mat.list, ref = all.class.list)
    
  }
  
  rownames(geno012) <- geno.names
  
  
  results <- list(geno012 = geno012, all.ref = all.class)
  
  return(results)
  
}


