##########
# QTL_R2 #
##########

#' QTL global and partial R squared
#' 
#' Computes the global and partial (adjusted) R squared of a  list of QTLs using
#' a linear model.
#' 
#' The function computes R squared statistics using a linear model. The extra
#' variance explained by a full model containing the QTL terms with respect
#' to a reduced model containing only the cross intercept terms and uses the
#' ratio between the residual sum of square of these two models:
#' R2 = 1-(RSS(f))/(RSS(r)).
#' 
#' Partial R squared for each individual QTL position can also be calculated.
#' Two types of partial R squared are returned. The first one
#' uses the difference between the R squared obtained with all QTL
#' positions and the R squared obtain with all position minus the ith one
#' (difference R squared). The second method used only the ith QTL position
#' in the model (single R squared).
#' 
#' For both global and partial R squared, it is possible to obtained adjusted
#' measurements taking the number of degrees of freedom into consideration using
#' an adaptation of the formula given by Utz et al. (2000):
#' R.adj = R-(z/(N-z-n.cr))*(1-R) where z is the total
#' number of estimated components of the genetic effect. N is the total number
#' of phenotypic information, and n.cr is the number of intercept (cross) terms.
#' 
#'
#' @param mppData An object of class \code{mppData}.
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function \code{\link{QTL_select}} or
#' vector of \code{character} marker positions names.
#' Default = NULL.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. For more details see
#' \code{\link{mpp_SIM}}. Default = "cr".
#' 
#' @param glb.only \code{Logical} value. If glb.only = TRUE, only the global and 
#' global adjusted R squared will be returned. Default = FALSE.
#'
#' @return Return:
#'
#' object of class \code{QR2Res} containing the following objects:
#'
#' \item{glb.R2 }{ Global R squared of all QTL terms.}
#'
#' \item{glb.adj.R2 }{ Global adjusted R squared of all QTL terms.}
#'
#' \item{part.R2.diff }{ Vector of partial R squared doing
#' the difference between the full model and a model minus the ith QTL.}
#' 
#' \item{part.adj.R2.diff }{ Vector of partial adjusted R squared doing
#' the difference between the full model and a model minus the ith QTL.}
#' 
#' \item{part.R2.sg }{ Vector of partial R squared using only the ith QTL.}
#' 
#' \item{part.adj.R2.sg }{ Vector of partial adjusted R squared using only the
#' ith QTL.}
#' 
#' 
#' @author Vincent Garin
#'
#' @references
#' 
#' Utz, H. F., Melchinger, A. E., & Schon, C. C. (2000). Bias and sampling error
#' of the estimated proportion of genotypic variance explained by quantitative
#' trait loci determined from experimental data in maize using cross validation
#' and validation with independent samples. Genetics, 154(4), 1839-1849.
#'
#' @seealso \code{\link{QTL_select}}, \code{\link{summary.QR2Res}}
#'
#' @examples
#'
#' data(mppData)
#' 
#' SIM <- mpp_SIM(mppData)
#' QTL <- QTL_select(Qprof = SIM, threshold = 3, window = 20)
#' Q_R2 <- QTL_R2(mppData = mppData, QTL = QTL, Q.eff = "cr")
#' summary(Q_R2)
#' 
#' 
#' @export
#'


QTL_R2 <- function(mppData, trait = 1, QTL = NULL, Q.eff = "cr",
                   glb.only = FALSE){
  
  # 1. check the data format
  ##########################
  
  check.model.comp(mppData = mppData, trait = trait, Q.eff = Q.eff, QTL = QTL,
                   fct = "R2")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  ### 2.2 Formation of the list of QTLs
  
  if(is.character(QTL)){
    
    Q.pos <- which(mppData$map[, 1] %in% QTL)
    
  } else {
    
    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])
    
  }
  
  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                   Q.eff = Q.eff, order.MAF = TRUE)
  
  n.QTL <- length(Q.list)
  
  # 3. Compute the R squared
  ##########################
  
  ### 3.1 Global adjusted and unadjusted linear R squared
  
  R2.all <- R2_lin(mppData = mppData, trait = t_val,
                   QTL = do.call(cbind, Q.list))
  
  R2 <- R2.all[[1]]
  R2.adj <- R2.all[[2]]
  
  
  if(glb.only) {
    
    QR2Res <- list(glb.R2 = R2, glb.adj.R2 = R2.adj)
    
    class(QR2Res) <- c("QR2Res", "list")
    
    return(QR2Res)
    
  } else {
    
    if(n.QTL > 1){
      
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
      
      QR2Res <- list(glb.R2 = R2,
                  glb.adj.R2 = R2.adj,
                  part.R2.diff = R2_i.dif,
                  part.adj.R2.diff = R2_i.dif.adj,
                  part.R2.sg = R2_i.sg,
                  part.adj.R2.sg = R2_i.sg.adj)
      
      class(QR2Res) <- c("QR2Res", "list")
      
      return(QR2Res)
      
    } else {
      
      names(R2) <- names(R2.adj) <- "Q1"
      
      QR2Res <- list(glb.R2 = R2,
                    glb.adj.R2 = R2.adj,
                    part.R2.diff = R2,
                    part.adj.R2.diff = R2.adj,
                    part.R2.sg = R2,
                    part.adj.R2.sg = R2.adj)
      
      class(QR2Res) <- c("QR2Res", "list")
      
      return(QR2Res)
    }
    
  }
  
}