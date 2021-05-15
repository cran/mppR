###############
# QTL_forward #
###############

#' Forward regression QTL model
#' 
#' Determines a multi-QTL model using a forward regression.
#' 
#' Forward regression to determine the a multi-QTL model. The function
#' selects successively QTL positions with -log10(pval) above the threshold.
#' Those positions are added as cofactors for following detection run.
#' The procedure stop when no more position has a -log10(pval) above the
#' threshold.
#'
#' @param mppData An object of class \code{mppData}
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#' 
#' @param Q.eff \code{Character} vector of possible QTL effects the user want to
#' test. Elements of Q.eff can be "cr", "par", "anc" or "biall". For details
#' look at \code{\link{mpp_SIM}}.
#' 
#' @param threshold \code{Numeric} value representing the -log10(p-value) threshold
#' above which a position can be considered as significant. Default = 4.
#' 
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 30.
#' 
#' @param n.cores \code{Numeric}. Specify here the number of cores you like to
#' use. Default = 1.
#' 
#' @param verbose \code{Logical} value indicating if the steps of the
#' forward regression should be printed. Default = TRUE.
#' 
#' 
#' @return Return:
#' 
#' \item{QTL }{\code{Data.frame} of class \code{QTLlist} with five columns :
#' 1) QTL marker names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-values).}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mpp_SIM}},
#'
#' @examples
#'
#' data(mppData)
#' 
#' QTL <- QTL_forward(mppData = mppData, Q.eff = "par")
#'
#' @export
#'

QTL_forward <- function(mppData = NULL, trait = 1, Q.eff,
                        threshold = 4, window = 30, n.cores = 1,
                        verbose = TRUE) {
  
  # 1. test data format
  #####################
  
  check.MQE(mppData = mppData, trait = trait, Q.eff = Q.eff,
            VCOV = 'h.err', n.cores = n.cores, fct = "QTL_forward")
  
  # initialize list of selected QTL and list of type of effects.
  
  QTL.list <- c()
  pos.ind <- 1
  
  # form cluster
  
  if(n.cores > 1){
    
    parallel <- TRUE
    cluster <- makeCluster(n.cores)
    
  } else {
    
    parallel <- FALSE
    cluster <- NULL
    
  }
  
  # 2. step 1: SIM
  ################
  
  
  SIM <- mpp_SIM_clu(mppData = mppData, trait = trait, Q.eff = Q.eff,
                     parallel = parallel, cluster = cluster)
  
  max.pval <- max(SIM$log10pval, na.rm = TRUE)
  
  
  if(max.pval > threshold){
    
    # store the detected QTL
    QTL.list <- rbind(QTL.list, SIM[which.max(SIM$log10pval), ])
    
    if(verbose){
      
      cat("\n")
      cat(paste("position", pos.ind))
      cat("\n")
      
    }
    
    pos.ind <- pos.ind + 1
    
    # 3. Step 2 : CIM with cofactors (QTL already selected)
    #######################################################
    
    rest.sign.pos <- TRUE
    
    cofactors <- QTL.list[, 1]
    
    while(rest.sign.pos) {
      
      CIM <- mpp_CIM_clu(mppData = mppData, trait = trait, Q.eff = Q.eff,
                        cofactors = cofactors, window = window,
                         plot.gen.eff = FALSE, parallel = parallel,
                         cluster = cluster)
      
      # remove candidate positions around the selected cofactors
      
      pos.QTL <- mppData$map[mppData$map[, 1] %in% QTL.list[, 1], c(2,4)]
      
      test.cof <- function(x, map, window) {
        
        t1 <- map$chr == as.numeric(x[1])
        t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
        (t1 & t2)
        
      }
      
      cof.part <- apply(X = pos.QTL, MARGIN = 1, FUN = test.cof,
                        map = mppData$map, window = window)
      
      CIM <- CIM[rowSums(cof.part) == 0, ]
      
      max.pval <- max(CIM$log10pval, na.rm = TRUE)
      
      if(max.pval > threshold){
        
        # store the detected QTL
        QTL.list <- rbind(QTL.list, CIM[which.max(CIM$log10pval), ])
        cofactors <- QTL.list[, 1]
        
        if(verbose){
          
          cat("\n")
          cat(paste("position", pos.ind))
          cat("\n")
          
        }
        
        pos.ind <- pos.ind + 1
        
      } else {
        
        rest.sign.pos <- FALSE
        
      }
      
    }
    
    
    if(n.cores > 1){stopCluster(cluster)}
    
    # return the final list of QTL
    
    QTL.list <- QTL.list[order(QTL.list$chr, QTL.list$pos.ind), ]
    
    class(QTL.list) <- c("QTLlist", "data.frame")
    
    return(QTL.list)
    
    
  } else {
    
    warning(paste("No position is above the threshold, the stepwise procedure",
                  "could not select any QTL position."))
    
    return(NULL)
    
  }
  
}