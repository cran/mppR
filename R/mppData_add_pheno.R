#####################
# mppData_add_pheno #
#####################

#' Add new phenotypic values to a mppData object
#'
#' Add the new phenotypic values contained in 'pheno' to a \code{mppData} object.
#' 
#' @param mppData An object of class \code{mppData}.
#' 
#' @param pheno \code{data.frame} with : 1) \code{character}
#' genotypes identifiers; 2) \code{numeric} trait values. \strong{The genotypes
#' identifiers must be identical to \code{mppData$geno.id}.}
#' 
#' @return Return:
#' 
#' \item{mppData}{New \code{mppData} object with new phenotypic values added.}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_mdf_pheno}}, \code{\link{subset.mppData}},
#' 
#' @examples
#' 
#' data(mppData)
#' pheno_new <- data.frame(geno.id = mppData$geno.id, ph1 = rnorm(498))
#'
#' mppData <- mppData_add_pheno(mppData = mppData, pheno = pheno_new)
#' 
#' @export
#' 

mppData_add_pheno <- function(mppData, pheno){
  
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
  
  if(any(ph_new_id %in% ph_ref_id)){
    
    prob_id <- ph_new_id[ph_new_id %in% ph_ref_id]
    s_mess <- paste('The following identifiers: ', paste(prob_id, collapse = ', '),
                    'are already used in mppData object pheno values.')
    
    stop(s_mess)
    
  }
  
  # add the new trait(s)
  mppData$pheno <- as.matrix(cbind(ph_ref, ph_new))
  
  return(mppData)
  
}