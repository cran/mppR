#########################
# print.summary.QeffRes #
#########################

#' Print summary.QeffRes object
#' 
#' @param x object of class \code{summary.QeffRes}
#' 
#' @param ... Ignored. 
#' 
#' @examples 
#' 
#' data(mppData)
#' SIM <- mpp_SIM(mppData)
#' QTL <- QTL_select(SIM)
#' QTL.effects <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "cr")
#' sum.QeffRes <- summary(QTL.effects)
#' print(sum.QeffRes)
#' 
#' @export
#' 

print.summary.QeffRes <- function(x, ...){
  
  # write the title
  
  cat("QTL effects")
  cat("\n")
  cat(rep("*", 11), sep = "")
  cat("\n")
  cat("\n")
  
  # write the number of QTL and total R squared
  
  cat(paste("Number of QTL(s):", length(x$QTL.effects)))
  cat("\n")
  
  for (i in 1:length(x$QTL.effects)){
    
    cat("\n")
    cat("\n")
    cat(paste("QTL", x$Q_ind[i]))
    cat("\n")
    cat(rep("-",5), sep = "")
    cat("\n")
    cat("\n")
    
    # QTL position information (with or without CI information)
    
    print.data.frame(x$QTL.info[i, ], row.names = FALSE)
    
    
    cat("\n")
    cat("\n")
    cat("QTL effect per cross or parent:")
    cat("\n")
    cat("\n")
    
    print(x$QTL.effects[[i]])
    
    
  }
  
}
