#############
# check.inf #
#############

# function to check and replace if some -log10(pval) are +/- Inf.

check.inf <- function(x){
  
  if (sum(is.infinite(x))>0){
    
    # find the maximum value after infinity
    
    max.sc <- max(x[-which(is.infinite(x))])
    
    # replace -infinity by 0
    
    x[((x<0) & is.infinite(x))] <- 0
    
    # replace +infinity by the maximum LOD score value
    
    x[((x>0) & is.infinite(x))] <- max.sc
    
  }
  
  return(x)
  
}
