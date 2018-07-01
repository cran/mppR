#######################
# parent_clusterCheck #
#######################

# Check the result of parent clustering
# 
# Detect marker or position for which all parents have been clustered into
# the same ancestral group (monomorphic). Rewrite these positions with one
# allele per parent.
# 
# @param par.clu \code{Interger matrix} representing the results of a parents
# genotypes clustering. The columns represent the parental lines and the rows
# the different markers or in between positions. At a particular
# position, parents with the same value are assumed to inherit from the same
# ancestor.
# 
# @return Return:
# 
# \code{List} containing the two following objects:
# 
# \item{par.clu}{Reorganised \code{par.clu} object with monomorphic position
# put back to the parental situation (1 allele per parent).}
# 
# \item{marker.prob}{Vector of marker(s) or in between positions where all
# parents where groupped into a single group}
# 
# @author Vincent Garin
# 
# 
# @export
# 


parent_clusterCheck <- function(par.clu){
  
  # 1. detect monomorphic positions
  #################################
  
  mono.pos <- apply(X = par.clu, MARGIN = 1,
                    FUN = function(x) length(unique(x)) == 1 )
  
  
  # 2. replace monomorphic positions with 1 allele per parents
  ############################################################
  
  n.mono <- sum(mono.pos)
  
  if(n.mono > 0){ # there are some monomorphic positions
    
    par.clu[mono.pos, ] <-  matrix(rep(1:dim(par.clu)[2], n.mono),
                                   nrow = n.mono, byrow = TRUE)
    
    marker.prob <- rownames(par.clu)[mono.pos]
    
  } else { # no monomorphic positions
    
    marker.prob <- NULL
    
  }
  
  return(list(par.clu = par.clu, marker.prob = marker.prob))
  
}