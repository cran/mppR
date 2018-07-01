##########
# MQE_R2 #
##########


# Global and partial R squared for multi-QTL effects 
# 
# Computes the global and partial (adjusted) R squared of a list of QTLs
# including positions with different types of QTL effects. For example
# position one is a parental effect, position two a bi-allelic QTL, etc.
# The type of effect of the QTL position are specified in \code{Q.eff}.
# 
# The R squared computation is done using a linear model
# corresponding to \code{VCOV = 'h.err'}. For more details about R squared
# computation and adjustement look at function \code{\link{QTL_R2}}.
#
# @param mppData An object of class \code{mppData}.
# 
# @param trait \code{Numerical} or \code{character} indicator to specify which
# trait of the \code{mppData} object should be used. Default = 1.
# 
# @param QTL Vector of \code{character} markers or in between marker positions
# names. Default = NULL.
# 
# @param Q.eff \code{Character} vector indicating for each QTL position the
# type of QTL effect among: "cr", "par", "anc" and "biall". For details look at
# \code{\link{mpp_SIM}}.
#
# @param glb.only \code{Logical} value. If \code{glb.only = TRUE}, only
# the global and global adjusted R squared will be returned. Default = FALSE.
# 
#
# @return Return:
#
# List containing the following objects:
#
# \item{glb.R2 }{ Global R squared of all QTL terms.}
#
# \item{glb.adj.R2 }{ Global adjusted R squared of all QTL terms.}
#
# \item{part.R2.diff }{ Vector of partial R squared doing
# the difference between the full model and a model minus the ith QTL.}
# 
# \item{part.adj.R2.diff }{ Vector of partial adjusted R squared doing
# the difference between the full model and a model minus the ith QTL.}
# 
# \item{part.R2.sg }{ Vector of partial R squared using only the ith QTL.}
# 
# \item{part.adj.R2.sg }{ Vector of partial adjusted R squared using only the
# ith QTL.}
# 
# 
# @author Vincent Garin
# 
# @seealso \code{\link{mpp_SIM}}, \code{\link{QTL_R2}}
#
# @examples
#
# data(mppData)
# 
# SIM <- mpp_SIM(mppData = mppData)
# QTL <- QTL_select(SIM)
# 
# MQE_R2(mppData = mppData, QTL = QTL[, 1], Q.eff = c("par", "anc", "biall"))
#
# @export
#


MQE_R2 <- function(mppData = NULL, trait = 1, QTL = NULL, Q.eff,
                   glb.only = FALSE){
  
  # 1. check the data format
  ##########################
  
  check.MQE(mppData = mppData, trait = trait, Q.eff = Q.eff,
            QTL = QTL, fct = "R2")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  
  ### 2.2 Formation of the list of QTL incidence matrices
  
  # order list of QTL positions
  
  Q.pos <- vapply(X = QTL,
                  FUN = function(x, mppData) which(mppData$map[, 1] == x),
                  FUN.VALUE = numeric(1), mppData = mppData)
  
  Q.ord <- data.frame(QTL, Q.eff, Q.pos, stringsAsFactors = FALSE)
  
  Q.ord <- Q.ord[order(Q.pos), ]
  
  QTL <- Q.ord[, 1]; Q.eff <- Q.ord[, 2]
  
  Q.pos <- which(mppData$map[, 1] %in% QTL)
  
  # form a list of QTL incidence matrices with different type of QTL effect.
  
  # function to produce different type of QTL incidence matricdes
  
  Q.list <- mapply(FUN = inc_mat_QTL, x = Q.pos, Q.eff = Q.eff,
                   MoreArgs = list(mppData = mppData, order.MAF = TRUE),
                   SIMPLIFY = FALSE)
  
  n.QTL <- length(Q.list)
  
  # 3. Compute the R squared
  ##########################
  
  ### 3.1 Global adjusted and unadjusted linear R squared
  
  R2.all <- R2_lin(mppData = mppData, trait = t_val,
                   QTL = do.call(cbind, Q.list))
  
  R2 <- R2.all[[1]]
  R2.adj <- R2.all[[2]]
  
  if(glb.only) {
    
    return(list(glb.R2 = R2, glb.adj.R2 = R2.adj))
    
  } else {
    
    if(n.QTL > 1){
      
      ### 3.2 Compute the partial R squared
      
      # functions to compute the R squared or all QTL minus 1 or only 1 QTL position
      
      part.R2.diff <- function(x, QTL, mppData, trait) {
        R2_lin(mppData = mppData, trait = trait,
               QTL = do.call(cbind, Q.list[-x]))
      }
      
      part.R2.sg <- function(x, QTL, mppData, trait) {
        R2_lin(mppData = mppData, trait = trait,
               QTL = do.call(cbind, Q.list[x]))
      }
      
      R2.dif <- lapply(X = 1:n.QTL, FUN = part.R2.diff, QTL = Q.list,
                       mppData = mppData, trait = t_val)
      
      R2_i.dif <- lapply(X = R2.dif, FUN = function(x) x[[1]])
      R2_i.dif.adj <- lapply(X = R2.dif, FUN = function(x) x[[2]])
      
      R2_i.dif <- R2 - unlist(R2_i.dif) # difference full model and model minus i
      R2_i.dif.adj <- R2.adj - unlist(R2_i.dif.adj)
      
      R2.sg <- lapply(X = 1:n.QTL, FUN = part.R2.sg, QTL = Q.list,
                      mppData = mppData, trait = t_val)
      
      R2_i.sg <- unlist(lapply(X = R2.sg, FUN = function(x) x[[1]]))
      R2_i.sg.adj <- unlist(lapply(X = R2.sg, FUN = function(x) x[[2]]))
      
      
      names(R2_i.dif) <- names(R2_i.dif.adj) <- paste0("Q", 1:n.QTL)
      names(R2_i.sg) <- names(R2_i.sg.adj) <- paste0("Q", 1:n.QTL)
      
      return(list(glb.R2 = R2,
                  glb.adj.R2 = R2.adj,
                  part.R2.diff = R2_i.dif,
                  part.adj.R2.diff = R2_i.dif.adj,
                  part.R2.sg = R2_i.sg,
                  part.adj.R2.sg = R2_i.sg.adj))
      
    } else {
      
      names(R2) <- names(R2.adj) <- "Q1"
      
      return(list(glb.R2 = R2,
                  glb.adj.R2 = R2.adj,
                  part.R2.diff = R2,
                  part.adj.R2.diff = R2.adj,
                  part.R2.sg = R2,
                  part.adj.R2.sg = R2.adj))
    }
    
  }
  
}