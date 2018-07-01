###################
# summary.QeffRes #
###################

#' Summary of \code{QeffRes} object
#' 
#' \code{summary} for object of class \code{QeffRes}.
#' 
#' @param object An object of class \code{QeffRes} obtained with
#' function \code{\link{QTL_gen_effects}}.
#' 
#' @param QTL \code{Numeric} vector indicating the QTL positions for which the
#' QTL effect must be printed. Default = NULL.
#' 
#' @param ... Ignored.
#' 
#' @seealso \code{\link{QTL_gen_effects}}
#' 
#' @examples
#' 
#' data(mppData)
#' SIM <- mpp_SIM(mppData)
#' QTL <- QTL_select(SIM)
#' QTL.effects <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "cr")
#' summary(QTL.effects)
#' 
#' @export
#' 

summary.QeffRes <- function(object, QTL = NULL, ...) {
  
  stopifnot(inherits(object, "QeffRes"))
  ans <- list()

  if(is.null(QTL)){
    
    ans$QTL.effects <- object[[1]]
    ans$QTL.info <- t(object[[2]][1:3, ])
    ans$QTL.info <- data.frame(ans$QTL.info[, 1],
                               as.numeric(ans$QTL.info[, 2]),
                               as.numeric(ans$QTL.info[, 3]),
                               stringsAsFactors = FALSE)
    colnames(ans$QTL.info) <- c("mk.names", "chr", "pos.cM")
    ans$Q_ind <- 1:length(object[[1]])
    
  } else {
    
    ans$QTL.effects <- object[[1]]
    ans$QTL.effects <- ans$QTL.effects[QTL]
    ans$QTL.info <- t(object[[2]][1:3, QTL])
    ans$QTL.info <- data.frame(ans$QTL.info[, 1],
                               as.numeric(ans$QTL.info[, 2]),
                               as.numeric(ans$QTL.info[, 3]),
                               stringsAsFactors = FALSE)
    colnames(ans$QTL.info) <- c("mk.names", "chr", "pos.cM")
    ans$Q_ind <- QTL
    
  }
  
  class(ans) <- "summary.QeffRes"
  ans
  
}