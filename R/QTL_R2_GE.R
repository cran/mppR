#############
# QTL_R2_GE #
#############

#' MPP GxE QTL R2
#'
#' Compute global and partial R2 statistics for MPP GxE QTL using a linear model.
#' The global R2 is the contribution of all QTL positions while the partial R2
#' is the specific contribution of an individual QTL position.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits (environments) should be used.
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function QTL_select() or
#' a vector of \code{character} marker positions names. Default = NULL.
#'
#' @param glb.only \code{Logical} value. If glb.only = TRUE, only the global and
#' global adjusted R squared will be returned. Default = FALSE.
#'
#' @return Return:
#'
#' \code{List} containing the global unadjusted R2, the global adjusted R2,
#' the partial unadjusted R2, and the partial adjusted R2.
#'
#' @author Vincent Garin
#'
#' @examples
#'
#' data(mppData_GE)
#'
#' Qpos <- c("PZE.105068880", "PZE.106098900")
#'
#' R2 <- QTL_R2_GE(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
#'                 QTL = Qpos)
#'
#' @export
#'

QTL_R2_GE <- function(mppData, trait, QTL = NULL, glb.only = FALSE){
  
  ##### Check data format and arguments #####
  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = "par", VCOV = 'ID',
                  QTL_ch = TRUE, QTL = QTL)
  
  ##### form the trait value ####
  nEnv <- length(trait)
  TraitEnv <- c(mppData$pheno[, trait])
  
  
  ##### form the list of QTLs #####
  if(is.character(QTL)){
    
    Q.pos <- which(mppData$map[, 1] %in% QTL)
    
    QTL <- mppData$map[mppData$map[, 1] %in% QTL, ]
    
  } else {
    
    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])
    
  }
  
  n.QTL <- length(Q.pos)
  
  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                   Q.eff = "par", order.MAF = TRUE)
  
  Q.list <- lapply(X = Q.list, FUN =  function(x, nEnv) diag(nEnv) %x% x,
                   nEnv = nEnv)
  
  names(Q.list) <- paste0("Q", 1:length(Q.list))
  
  #### Compute the global R2 ####
  R2.all <- R2_lin_GE(mppData = mppData, trait = TraitEnv, nEnv = nEnv,
                      QTL = do.call(cbind, Q.list))
    
    R2 <- R2.all[[1]]
    R2.adj <- R2.all[[2]]
    
  
  #### Partial R2 #####
  
  if(glb.only) {
    
    QR2Res <- list(glb.R2 = R2, glb.adj.R2 = R2.adj)
    
  } else {
    
    if(n.QTL > 1){
        
    R2.dif <- lapply(X = 1:n.QTL, FUN = part.R2.diff.lin, Q.list = Q.list,
                    mppData = mppData, trait = TraitEnv, nEnv = nEnv)
        
    R2_i.dif <- lapply(X = R2.dif, FUN = function(x) x[[1]])
    R2_i.dif.adj <- lapply(X = R2.dif, FUN = function(x) x[[2]])
        
    R2_i.dif <- R2 - unlist(R2_i.dif) # difference full model and model minus i
    R2_i.dif.adj <- R2.adj - unlist(R2_i.dif.adj)
        
    names(R2_i.dif) <- names(R2_i.dif.adj) <- paste0("Q", 1:n.QTL)
        
    QR2Res <- list(glb.R2 = R2,
                   glb.adj.R2 = R2.adj,
                   part.R2.diff = R2_i.dif,
                   part.adj.R2.diff = R2_i.dif.adj)
      
    } else {
      
      names(R2) <- names(R2.adj) <- "Q1"
      
      QR2Res <- list(glb.R2 = R2,
                     glb.adj.R2 = R2.adj,
                     part.R2.diff = R2,
                     part.adj.R2.diff = R2.adj)
      
    }
    
  }
  
  class(QR2Res) <- c("QR2Res", "list")
  
  return(QR2Res)
  

}
