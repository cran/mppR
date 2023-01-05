#####################
# mppData_mdf_pheno #
#####################

#' Modify the phenotypic values of a mppData object
#'
#' Modify the phenotypic values of a \code{mppData} object.
#' 
#' @param mppData An object of class \code{mppData}.
#' 
#' @param pheno Two columns \code{data.frame} with : 1) \code{character}
#' genotypes identifiers; 2) \code{numeric} trait values. \strong{The genotypes
#' identifiers must be identical to \code{mppData$geno.id}. The trait value identifiers
#' must correspond to a trait already in the mppData object.}
#' 
#' @return Return:
#' 
#' \item{mppData}{New \code{mppData} object with modified phenotypic values added.}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_add_pheno}}, \code{\link{subset.mppData}},
#' 
#' @examples
#' 
#' data(mppData)
#' pheno_new <- data.frame(geno.id = mppData$geno.id, ULA = rnorm(498))
#'
#' mppData <- mppData_mdf_pheno(mppData = mppData, pheno = pheno_new)
#' 
#' @export
#' 

mppData_mdf_pheno <- function(mppData, pheno){
  
  if(!is_mppData(mppData)) {
    stop("'mppData' must be of class ", dQuote("mppData"))
  }
  
  if(!identical(pheno[, 1], mppData$geno.id)) {
    stop("the genotypes identifiers of 'mppData' and 'pheno' are not identical")
  }
  
  ph_ref <- mppData$pheno
  ph_new <- pheno[, 2:ncol(pheno), drop = FALSE]
  ph_ref_id <- colnames(ph_ref)
  ph_new_id <- colnames(ph_new)
  
  if(any(!(ph_new_id %in% ph_ref_id))){
    
    s_mess <- paste('pheno colnames do not match pheno colnames of the mppData object. It must be one or some of:',
                    paste(ph_ref_id, collapse = ', ' ))
    
    stop(s_mess)
    
  }
  
  # modify the trait(s) value(s)
  ph_new <- as.matrix(ph_new)
  ph_ref[, colnames(ph_ref) %in% colnames(ph_new)] <- ph_new
  
  mppData$pheno <- ph_ref
  
  return(mppData)
  
}