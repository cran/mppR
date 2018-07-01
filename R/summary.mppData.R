###################
# summary.mppData #
###################

#' Summary of \code{mppData} object
#' 
#' \code{summary} for object of class \code{mppData}.
#' 
#' @param object An object of class \code{mppData}.
#' 
#' @param ... Ignored.
#' 
#' @examples
#' 
#' data(mppData)
#' summary(mppData)
#' 
#' @export
#' 

summary.mppData <- function(object, ...) {
  
  stopifnot(inherits(object, "mppData"))
  ans <- list()
  
  par.per.cross <- t(object$par.per.cross)
  colnames(par.per.cross) <- attr(table(object$cross.ind), "dimnames")[[1]]
  par.per.cross <- rbind(par.per.cross, table(object$cross.ind))
  colnames(par.per.cross) <- rep(" ", dim(par.per.cross)[2])
  rownames(par.per.cross) <- c("Crosses", "Parent 1", "Parent 2", "N")
  
  # genotype info
  
  ans$typePop <- object$type
  ans$Ngeno <- dim(object$pheno)[1]
  ans$par.per.cross <- par.per.cross
  
  # phenotype info
  
  ans$phenoName <- colnames(object$pheno)
  
  miss_fct <- function(x){ round((1 - (sum(is.na(x))/ length(x)))*100, 1) }
  
  ans$phenoPer <- apply(X = object$pheno, MARGIN = 2, FUN = miss_fct)
  
  # map information

    ans$mkNb <- dim(object$map)[1]
    chr_ind <- factor(x = object$map[, 2], levels = unique(object$map[, 2]))
    ans$mkChr <- table(chr_ind)
  
  class(ans) <- "summary.mppData"
  ans
  
}