##########################
# parent_cluster.mppData #
##########################

#' Parent clustering for \code{mppData} objects
#' 
#' Integrate the parent clustering information to the mppData
#' object. The parent clustering is necessary to compute the ancestral model.
#' If the parent clustering step is skipped, the ancestral model can not be
#' used but the other models (cross-specific, parental, and bi-allelic) can
#' still be computed.
#' 
#' At a single marker position, two parents can be grouped into a similar
#' ancestral classes if we assume that they receive there allele from a common
#' ancestor. The parent clustering information (\code{par.clu}) describe parental
#' relatedness and which parent belong to which ancestral group. For example,
#' at marker i, we could have five parents (pA, pB, pC, pD, pE) and the following
#' clustering information (1, 2, 1, 2, 3). This means that pA and pC received
#' their allele from the same ancestor (A1). pB and pD also have a shared
#' ancestor (A2) who is different from (A1). And pE was not included in any
#' group and can be seen as an independent ancestral group (A3).
#' 
#' The parent clustering information is provided via \code{par.clu}. It is an
#' \code{interger matrix} with markers in row and parents in columns.
#' At a particular marker position, parents with the same value are assumed to
#' inherit from the same ancestor. for more details, see \code{\link{par_clu}}.
#' 
#' The marker positions that are considered as monomorphic given the parent
#' clustering information are set back to one allele per parent to still allow
#' the computation of the QTL allelic effect at those positions later.
#' 
#' The parent clustering can be performed using the R package
#' 'clusthaplo' that can be found there:
#' \url{https://cran.r-project.org/src/contrib/Archive/clusthaplo/}.
#' The 'clusthaplo' option is not integrated in this version of mppR. However,
#' a version of mppR with function calling clusthaplo can be found on github
#' \url{https://github.com/vincentgarin/mppR} (branch master).
#' 
#' @param mppData  An object of class \code{mppData}. the \code{mppData} must
#' have been processed using: \code{\link{create.mppData}},
#' \code{\link{QC.mppData}}, \code{\link{IBS.mppData}},
#' and \code{\link{IBD.mppData}}.
#' 
#' @param par.clu \code{Interger matrix} representing the results of a
#' parents genotypes clustering. The columns represent the parental lines and
#' the rows the markers. The columns names must be the same as the parents
#' list of the mppData object. The rownames must be the same as the map marker
#' list of the mppData object. At a particular position, parents with the same
#' value are assumed to inherit from the same ancestor. for more details,
#' see \code{\link{par_clu}}. Default = NULL.
#' 
#' @return
#' 
#' An increased \code{mppData} object containing the the same elements
#' as the \code{mppData} object provided as argument and the
#' following new elements:
#' 
#' \item{par.clu}{\code{Integer matrix} with rows representing markers and
#' columns corresponding to the parents. At a single marker position, parents
#' with the same value were clustered in the same ancestral group.}
#' 
#' \item{n.anc}{Average number of ancestral clusters along the genome.}
#' 
#' \item{mono.anc}{Positions for which the ancestral clustering was monomorphic.}
#' 
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{create.mppData}}, \code{\link{QC.mppData}},
#' \code{\link{IBS.mppData}}, \code{\link{IBD.mppData}}, \code{\link{par_clu}}
#' 
#' @examples
#' 
#' data(mppData_init)
#' data(par_clu)
#' 
#' mppData <- QC.mppData(mppData_init)
#' mppData <- IBS.mppData(mppData = mppData)
#' 
#' mppData <- IBD.mppData(mppData = mppData, type = 'RIL',
#'                        type.mating = 'selfing')
#'                        
#' mppData <- parent_cluster.mppData(mppData = mppData, par.clu  = par_clu)                         
#'                        
#' 
#' @export


parent_cluster.mppData <- function(mppData, par.clu = NULL){

  
  # check the format of the data
  #################################
  
  if(!is_mppData(mppData)){
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  # test if correct step in the mppData processing
  
  if(!(mppData$status %in% c('IBD', 'complete'))){
    
    stop("you have to process 'mppData' in a strict order: ",
         "create.mppData, QC.mppData, IBS.mppData, IBD.mppData, ",
         "parent_cluster.mppData. You can only use parent_cluster.mppData ",
         "after create.mppData, QC.mppData, IBS.mppData, and IBD.mppData")
    
  }
  
  
  # check the content of par.clu
    
    if(!is.matrix(par.clu)){
      
      stop("'par.clu' argument is not a matrix")
      
    }
    
    if(!is.integer(par.clu)){
      
      stop("'par.clu' is not integer")
      
    }
    
    # list parent
    
    new_par <- colnames(par.clu)
    
    if(!all(new_par %in% mppData$parents)) {
      
      wrong.par <- new_par[!(new_par %in% mppData$parents)]
      pbpar <- paste(wrong.par, collapse = ", ")
      
      message <- sprintf(ngettext(length(wrong.par),
                                  "the following parent %s is not present in 'mppData'",
                                  "the following parents %s are not present in 'mppData'"),
                         pbpar)
      
      stop(message)
      
    }
    
    # list markers
    
    if(!identical(rownames(par.clu), mppData$map[, 1])){
      
      stop("the markers of 'par.clu' and 'mppData' are not identical")
      
    }
    
    # Check monomorphism in par.clu
    ###############################
    
    par.clu <- par.clu[, mppData$parents]
    
    nb.cl <- apply(X = par.clu, MARGIN = 1,
                   FUN = function(x) length(unique(x)))
    
    av.cl <- mean(nb.cl)
    
    par_clu <- parent_clusterCheck(par.clu = par.clu)
    
    # Calculate the number of ancestral cluster
    ###########################################
    
    mppData$par.clu <- par_clu[[1]]
    
    mppData$n.anc <- av.cl
    
    mppData$mono.anc <- par_clu[[2]]
    
    mppData$status <- 'complete'
    
    mppData$geno.off <- NULL
    
    class(mppData) <- c("mppData", "list")
    
    return(mppData)
    
  
}