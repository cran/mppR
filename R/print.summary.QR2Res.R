#########################
# print.summary.QR2Res #
#########################

#' Print summary.QR2Res object
#' 
#' @param x object of class \code{summary.QR2Res}
#' 
#' @param ... Ignored. 
#' 
#' @examples 
#' 
#' data(mppData)
#' SIM <- mpp_SIM(mppData)
#' QTL <- QTL_select(SIM)
#' Q_R2 <- QTL_R2(mppData = mppData, QTL = QTL, Q.eff = "cr")
#' sum.QR2Res <- summary(Q_R2)
#' print(sum.QR2Res)
#' 
#' @export
#' 

print.summary.QR2Res <- function(x, ...){
  
  # write the title
  
  cat("QTL R2")
  cat("\n")
  cat(rep("*", 6), sep = "")
  cat("\n")
  cat("\n")
  
  if(x[[1]] == 'glb.only'){
    
    cat(paste("Global R2:", round(x$R2[[1]], 1)))
    cat("\n")
    cat("\n")
    cat(paste("Global adjusted R2:", round(x$R2[[2]], 1)))
    
  } else {
    
    cat(paste("Global R2:", round(x$R2[[1]], 1)))
    cat("\n")
    cat("\n")
    cat(paste("Global adjusted R2:", round(x$R2[[2]], 1)))
    cat("\n")
    cat("\n")
    
    # individual QTL results
    
    ind_Q_res <- data.frame(round(x$R2$part.adj.R2.sg, 1),
                            round(x$R2$part.adj.R2.diff, 1))
    
    colnames(ind_Q_res) <- c('part. adj. sg. R2', 'part. adj. diff. R2')
    
    print.data.frame(ind_Q_res)
    
    
  }
  
}
