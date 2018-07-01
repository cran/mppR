#########################
# print.summary.mppData #
#########################

#' Print summary.mppData object
#' 
#' @param x object of class \code{summary.mppData}
#' 
#' @param ... Ignored. 
#' 
#' @examples 
#' 
#' data(mppData)
#' sum.mppData <- summary(mppData)
#' print(sum.mppData)
#' 
#' @export
#' 

print.summary.mppData <- function(x, ...){
  
  cat("object of class 'mppData'", "\n \n")
  cat("\t Type of population: ", x$typePop, "\n \n")
  cat("\t No. Genotypes: ", x$Ngeno, "\n")
  print(x$par.per.cross, quote = FALSE)
  cat("\n")
  cat("\t Phenotype(s): ", paste(x$phenoName), "\n")
  cat("\t Percent phenotyped: ", paste(x$phenoPer), "\n \n")
  cat("\t Total marker: ", x$mkNb, "\n")
  cat("\t No. markers:", x$mkChr)
  
}
