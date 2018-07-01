##################
# summary.QR2Res #
##################

#' Summary of \code{QR2Res} object
#' 
#' \code{summary} for object of class \code{QR2Res}.
#' 
#' @param object An object of class \code{QR2Res} obtained with
#' function \code{\link{QTL_R2}}.
#' 
#' @param ... Ignored.
#' 
#' @seealso \code{\link{QTL_R2}}
#' 
#' @examples
#' 
#' data(mppData)
#' SIM <- mpp_SIM(mppData)
#' QTL <- QTL_select(SIM)
#' Q_R2 <- QTL_R2(mppData = mppData, QTL = QTL, Q.eff = "cr")
#' summary(Q_R2)
#' 
#' @export
#' 

# data(mppData)
# 
# SIM <- mpp_SIM(mppData)
# QTL <- QTL_select(Qprof = SIM, threshold = 3, window = 20)
# Q_R2 <- QTL_R2(mppData = mppData, QTL = QTL, Q.eff = "cr")
# 
# object <- Q_R2
# QTL = c(1, 3)

summary.QR2Res <- function(object, ...) {
  
  stopifnot(inherits(object, "QR2Res"))
  ans <- list()
  
  if(length(object) == 2){
    
    ans$type <- 'glb.only'
    
  } else {
    
    ans$type <- 'complete'
    
    
  }
  
  ans$R2 <- object
  
  class(ans) <- "summary.QR2Res"
  ans
  
}