###############
# QTL_pred_R2 #
###############

#' Predicted QTL global and partial R squared
#' 
#' Compute predicted R squared in a validation set using QTLs detected in a
#' training set. These values are corrected by the heritability \code{her}.
#'
#' Compute QTLs predicted R squared in a validation set  (\code{mppData.vs}).
#' These QTLs have been previously detected in a training set
#' (\code{mppData.ts}). The global R squared (R2 = cor(y.ts,y.pred.ts)^2) is
#' obtained using the Pearson squared correlation between the observed trait
#' values in the validation set (y.vs) and predicted values using estimated QTL
#' effects in the training set (y.pred.vs = X.vs * B.ts).
#' 
#' After that the values are corrected by the general or within cross
#' heritability \code{her}. By default \code{her = 1} which means that the
#' R squared represent the proportion of explained phenotypic variance. The
#' values are returned per cross (\code{R2.cr}) or averaged at the population
#' level (\code{glb.R2}).
#' 
#' Partial R squared statistics are also calculated for each individual position.
#' The partial R squared are computed by making the difference between the
#' global R squared and the R squared computed without the ith position.
#'
#' @param mppData.ts An object of class \code{mppData} for the training set.
#' 
#' @param mppData.vs An object of class \code{mppData} for the validation set.
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. For more details see
#' \code{\link{mpp_SIM}}. Default = "cr".
#' 
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function \code{\link{QTL_select}} or
#' vector of \code{character} marker positions names.
#' Default = NULL.
#'
#' @param her \code{Numeric} value between 0 and 1 representing the heritability
#' of the trait. \code{her} can be a single value or a vector specifying each
#' within cross heritability. Default = 1.
#' 
#' @return Return:
#' 
#' \code{List} containing the following objects:
#'
#' \item{glb.R2 }{Global predicted R squared corrected for the heritability
#' of all QTL terms. Doing the average of the within cross predicted R squared
#' (R2.cr)}
#' 
#' \item{R2.cr}{Within cross predicted R squared corrected for the heritability}
#'
#' \item{part.R2.diff }{ Vector of predicted partial R squared corrected
#' for the heritability doing the difference between the full model and a model
#' minus the ith QTL.}
#' 
#' 
#' @author Vincent Garin
#' 
#' 
#' @seealso \code{\link{QTL_R2}}, \code{\link{QTL_select}}
#' 
#' @examples
#' 
#' data(mppData)
#' 
#' folds <- CV_partition(cross.ind = mppData$cross.ind, k = 5)
#' 
#' mppData.ts <- subset(x = mppData, gen.list = folds[[1]]$train.set)
#' 
#' mppData.vs <- subset(x = mppData, gen.list = folds[[1]]$val.set)
#' 
#' SIM <- mpp_SIM(mppData = mppData)
#' QTL <- QTL_select(SIM)
#'
#' QTL_pred_R2(mppData.ts = mppData.ts, mppData.vs = mppData.vs, QTL = QTL)
#' 
#' @export
#' 


QTL_pred_R2 <- function(mppData.ts, mppData.vs, trait = 1, Q.eff = "cr",
                        QTL = NULL, her = 1) {
  
  # 1. test data format
  #####################
  
  check.model.comp(Q.eff = Q.eff, trait = trait, VCOV = "h.err", QTL = QTL,
                   mppData.ts = mppData.ts, fct = "R2_pred")

  
  if(is.character(QTL)){ n.QTL <- length(QTL) } else { n.QTL <- dim(QTL)[1] }
  
  # trait values
  
  t_val <- sel_trait(mppData = mppData.vs, trait = trait)
  
  
  # 2. obtain the genetic effects (Betas)
  #######################################
  
  if((Q.eff == "biall") || (Q.eff == "cr")){
    if(!is.null(mppData.ts$geno.par)) {mppData.ts$geno.par <- NULL}
    if(!is.null(mppData.ts$geno.par.clu)) {mppData.ts$geno.par.clu <- NULL}
  }
  
  effects <- QTL_gen_effects(mppData = mppData.ts, trait = trait, QTL = QTL,
                            Q.eff = Q.eff)[[1]]
  
  # need to re-order the row of the effects according to the parent list
  
  if((Q.eff == "par") || (Q.eff == "anc")){
    
    effects <- lapply(X = effects, function(x, ind) x[ind, ],
                      ind = mppData.vs$parents)
    
  }
  
  
  B.ts <- lapply(X = seq_along(effects), FUN = function(x, Qeff) Qeff[[x]][, 1],
                 Qeff = effects)
  
  # 3. obtain the QTL incidence matrices of the positions (X.vs)
  ##############################################################
  
  if(is.character(QTL)){
    
    Q.pos <- which(mppData.vs$map[, 1] %in% QTL)
    
  } else {
    
    Q.pos <- which(mppData.vs$map[, 1] %in% QTL[, 1])
    
  }
  
  # switch Qeff for parental QTL incidence matrix if ancestral model
  
  if(Q.eff == "anc") Q.eff.part <- "par" else Q.eff.part <- Q.eff
  
  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData.vs,
                   Q.eff = Q.eff.part)
  
  # 4. Predicted R squared computation cor(y.vs, X.vs * B.ts)^2
  ##############################################################
  
  # global R squared
  
  R2 <- R2_pred(mppData.vs = mppData.vs, y.vs = t_val, B.ts = B.ts,
                Q.list = Q.list, her = her)
  
  # partial R2
  
  if (n.QTL > 1) {
    
    part.R2.diff <- function(x, mppData.vs, y.vs, B.ts, Q.list, her) {
      R2_pred(mppData.vs = mppData.vs, y.vs = t_val, B.ts = B.ts[-x],
              Q.list =  Q.list[-x], her = her)[[1]]
    }
    
    
    R2_i.dif <- lapply(X = 1:n.QTL, FUN = part.R2.diff, mppData.vs = mppData.vs,
                       y.vs = t_val, B.ts = B.ts, Q.list = Q.list, her = her)
    
    R2_i.dif <- R2[[1]] - unlist(R2_i.dif) # difference full model and model minus i
    
    
    names(R2_i.dif) <- paste0("Q", 1:n.QTL)
    
    
    return(list(glb.R2 = R2[[1]], R2.cr = R2[[2]],
                part.R2.diff = R2_i.dif))
    
  } else {
    
    R2.Q1 <- R2[[1]]
    names(R2.Q1) <- "Q1"
    
    return(list(glb.R2 = R2[[1]], R2.cr = R2[[2]],
                part.R2.diff = R2.Q1))
    
  }
  
}