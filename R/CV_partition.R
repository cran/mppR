################
# CV_partition #
################

#' Cross validation partition
#' 
#' Partition the genotype indices into training and validation sets for
#' cross-validation (CV).
#' 
#' The genotype indices are randomly assigned within cross to k subsets (folds).
#' Then each subset is used once as validation set, the remaining data go in the
#' training set.
#' 
#' @param cross.ind \code{Character} vector with the same length as the number
#' of genotypes which specifies to which cross each genotype belongs.
#' 
#' @param k \code{Numeric} value representing the number of subsets (fold) into
#' which data are spread within cross. Default = 5.
#' 
#' @return Return:
#' 
#' \item{fold}{\code{List} of k lists (one for each fold). Each fold list
#' contains two vectors with genotypes indices of the training (\code{$train.set}) and
#' the validation set (\code{$val.set}).}
#' 
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mpp_CV}}
#' 
#' @examples
#' 
#' data(mppData)
#'
#' part.cv <- CV_partition(cross.ind = mppData$cross.ind, k = 5)
#'
#' part.cv[[1]]$train.set
#' part.cv[[1]]$val.set
#'
#' @export
#'


CV_partition <- function(cross.ind, k = 5) {
  
  # 1. define a general indices
  #############################
  
  
  gen.indices <- 1:length(cross.ind)
  
  
  # split the ith cross into k subsets
  ####################################
  
  cross.id <- unique(cross.ind)
  part.cross <- vector(mode = "list", length = length(cross.id))
  
  for(i in seq_along(cross.id)){
   
    cross.i <- gen.indices[cross.ind == cross.id[i]]
    
    int.div <- length(cross.i) %/% k
    rest <- length(cross.i) %% k
    
    within.part <- c(rep(int.div, (k-rest)), rep((int.div + 1), (rest)))
    
    part.cross[[i]] <- split(x = sample(cross.i),
                          f = as.factor(rep(letters[1:k], within.part))) 
    
  }
  
  
  # rearrange the subsets into training and validation set
  ########################################################
  
  partition <- vector(mode = "list", length = k)
  
  for(i in 1:k){
    
    val.set <- c()
    
    for(j in seq_along(part.cross)){
      
     val.set <- c(val.set, part.cross[[j]][[i]]) 
      
    }
    
    val.set <- sort(val.set)
    
    train.set <- gen.indices[-val.set]
    
    
    partition[[i]] <- list(train.set = train.set, val.set = val.set)
    
  }
 
  return(partition)
   
}