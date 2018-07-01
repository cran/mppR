##################
# QTL_gen_effects #
##################

#' QTL genetic effects
#' 
#' Computes a multi-QTL model with a list of QTL candidates (\code{QTL}) and
#' return the decomposed QTL effects per cross or per parents.
#' 
#' This function computes for each QTL position the genetic effects of the
#' cross, parental, ancestral or SNP allele components. For the cross-specific
#' model (\code{Q.eff = "cr"}), the genetics effects represent the substitution
#' effect of an single allele from the parent 2 (or B) with respect to an allele
#' coming from the parent 1 or A. All effects are given in absolute value with
#' the parent that carries the positive allele.
#' 
#' For the parental and the ancestral model (\code{Q.eff = "par" or "anc"}), it
#' is possible to estimate maximum n-1 parental or ancestral alleles per
#' interconnected part of the design. For these two models, one
#' parental (ancestral) allele is set as reference per interconnected part of the
#' design. Effects of the other alleles are estimated as deviation with respect
#' to the reference. Connected parts of the design can be determined using Weeks
#' and Williams (1964) method (\code{\link{design_connectivity}}). By default,
#' the reference allele is the most frequent one. The user can also specify a
#' parental allele that will be used as reference using the argument
#' \code{ref.par}. This option is only available if the MPP design is composed
#' of a unique connected part.
#' 
#' For the parental and ancestral model it is also possible to estimate the QTL
#' effects using a sum to zero constraint \code{sum_zero = TRUE}. In that case,
#' the effects of the different parental (ancestral) allele will represent the
#' deviation with respect to the average trait value.
#' 
#' For the bi-allelic model (\code{Q.eff = "biall"}), the genetic effects
#' represent the effects of a single allele copy of the least frequent allele.
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
#' @param ref.par Optional \code{Character} expression defining the parental
#' allele that will be used as reference for the parental model. For the
#' ancestral model, the ancestral class containing the reference parent will be
#' set as reference. \strong{This option can only be used if the MPP design is
#' composed of a unique connected part}. Default = NULL.
#' 
#' @param sum_zero Optional \code{Logical} value specifying if the QTL effect of
#' a parental or an ancestral model should be calculated using the sum to zero
#' constraint. Default = FALSE.
#' 
#' 
#' @return Return:
#' 
#' Object of class \code{QeffRes} containing the following elements:
#'
#' \item{Qeff}{\code{List} of \code{data.frame} (one per QTL) containing the
#' following information:
#' 
#' \enumerate{
#' 
#' \item{QTL genetic effects per cross or parent.}
#' \item{Standard error of the QTL effects.}
#' \item{Test statistics of the effects (t-test or Wald statistic).}
#' \item{P-value of the test statistics.}
#' \item{Significance of the QTL effects.}
#' \item{For cross-specific model, parent with the positive additive effects.}
#' \item{For parental and ancestral model, indicator of connected part of the
#' design and reference.}
#' \item{Allele scores of the parents if \code{geno.par} is non NULL
#' in the \code{mppData} object.}
#' 
#' }
#' 
#' }
#' 
#' \item{tab.Qeff}{\code{data.frame} with one column per QTL giving the
#' QTL genetic effects per cross or per parent with its significance. The
#' first two rows indicate the chromosome and the position in cM of each
#' QTL.}
#' 
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{QTL_select}}, \code{\link{summary.QeffRes}}
#' 
#' @references 
#' 
#' Weeks, D. L., & Williams, D. R. (1964). A note on the determination of
#' connectedness in an N-way cross classification. Technometrics, 6(3), 319-324.
#' 
#' @examples
#' 
#' data(mppData)
#' 
#' # QTL candidates
#' 
#' SIM <- mpp_SIM(mppData)
#' QTL <- QTL_select(SIM)
#' 
#' # Cross-specific model
#' 
#' QTL.effects <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "cr")
#' summary(QTL.effects)
#' 
#' # Parental model
#' 
#' QTL.effects <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "par")
#' summary(QTL.effects)
#' 
#' # Ancestral model
#' 
#' QTL.effects <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "anc")
#' summary(QTL.effects)
#' 
#' # Bi-allelic model
#' 
#' QTL.effects <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "biall")
#' summary(QTL.effects)
#' 
#' @export 
#'


QTL_gen_effects <- function(mppData, trait = 1,QTL = NULL, Q.eff = "cr",
                           ref.par = NULL, sum_zero = FALSE) {
  
  # 1. Check data format
  ######################
  
  check.model.comp(mppData = mppData, trait = trait, Q.eff = Q.eff,
                   VCOV = 'h.err', QTL = QTL, ref.par = ref.par,
                   sum_zero = sum_zero, fct = "QTLeffects")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  ### 2.3 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.4 Formation of the list of QTL
  
  if(is.character(QTL)){
    
    Q.pos <- which(mppData$map[, 1] %in% QTL)
    pos.info <- t(data.frame(mppData$map[Q.pos, c(1, 2, 4)]))
    
  } else {
    
    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])
    pos.info <- t(data.frame(QTL[, 1], QTL[, 2], QTL[, 4]))
    
  }
  
  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                  Q.eff = Q.eff)
  
  names(Q.list) <- paste0("Q", 1:length(Q.list))
  
  if(!sum_zero){
    
    ### 2.7 For the parental and ancestral model organise the matrix to get the
    # desired constraint.
    
    Q.eff_temp <- rep(Q.eff, length(Q.list))
    
    order.Qmat <- mapply(FUN = IncMat_QTL_MAF, QTL = Q.list,
                         Q.eff_i = Q.eff_temp, Q.pos_i = Q.pos,
                         MoreArgs = list(mppData = mppData, ref.par = ref.par),
                         SIMPLIFY = FALSE)
    
    
    Q.list <- lapply(X = order.Qmat, FUN = function(x) x$QTL)
    allele_order <- lapply(X = order.Qmat, FUN = function(x) x$allele_order)
    con.ind <- lapply(X = order.Qmat, FUN = function(x) x$con.ind)
    
  } else {
    
    # re-organise the QTL incidence matrices adding the constraint for sum to 0
    
    IncMat_ext <- IncMat_sum0_const(mppData = mppData, Q.eff = Q.eff,
                                    Q.list = Q.list, Q.pos = Q.pos,
                                    cross.mat = cross.mat, trait = t_val)
    
    Q.list <- IncMat_ext$Q.list
    names(Q.list) <- paste0("Q", 1:length(Q.list))
    cross.mat <- IncMat_ext$cross.mat
    t_val <- IncMat_ext$trait
    con.ind <- IncMat_ext$con.ind
    allele_order <- lapply(X = Q.list, FUN = colnames)
    
  }
  
  
  # 3. model computation
  ######################
  
  model <- QTLModelQeff(mppData = mppData, trait = t_val, cross.mat = cross.mat,
                        Q.list = Q.list, VCOV = 'h.err')
  
  
  # 4. data processing
  ####################
  
  
  results <- Qeff_res_processing(model = model, mppData = mppData,
                                 cross.mat =  cross.mat, Q.list = Q.list,
                                 QTL = QTL, Q.eff = Q.eff, VCOV = 'h.err',
                                 allele_order = allele_order, con.ind = con.ind)
  
  
  names(results) <- paste0("Q", 1:length(results))
  
  # table with all QTL effect per cross or per parents.
  
  if (Q.eff == "cr"){
    
    Qeff_sign <- lapply(results, `[`, c(1, 5))
    Qeff_sign <- lapply(Qeff_sign, function(x) paste(round(x[, 1], 3), x[, 2]))
    table.QTL <- data.frame(Qeff_sign)
    colnames(pos.info) <- paste0("Q", 1:length(results))
    table.QTL <- rbind.data.frame(pos.info, table.QTL)
    rownames(table.QTL) <- c("mk.names", "chr", "pos.cM", unique(mppData$cross.ind))
    
  } else if (Q.eff == "biall"){
    
    Qeff_sign <- lapply(results, `[`, c(1, 5))
    Qeff_sign <- lapply(Qeff_sign, function(x) paste(round(x[, 1], 3), x[, 2]))
    table.QTL <- data.frame(Qeff_sign)
    colnames(pos.info) <- paste0("Q", 1:length(results))
    table.QTL <- rbind.data.frame(pos.info, table.QTL)
    
    if (is.null(mppData$geno.par)){
      
      rownames(table.QTL) <- c("mk.names", "chr", "pos.cM", "Q.eff")
      
    } else {
      
      rownames(table.QTL) <- c("mk.names", "chr", "pos.cM", mppData$parents)
      
    }
    
  } else { # parental or ancestral
    
    Qeff_sign <- lapply(results, `[`, c(1, 5))
    
    # order according to parents
    
    Qeff_sign <- lapply(X = Qeff_sign, function(x, ind) x[ind, ],
                        ind = mppData$parents)
    
    Qeff_sign <- lapply(Qeff_sign, function(x) paste(round(x[, 1], 3), x[, 2]))
    table.QTL <- data.frame(Qeff_sign)
    colnames(pos.info) <- paste0("Q", 1:length(results))
    table.QTL <- rbind.data.frame(pos.info, table.QTL)
    rownames(table.QTL) <- c("mk.names", "chr", "pos.cM", mppData$parents)
    
  }
  
  QeffRes <- list(Qeff = results, tab.Qeff = table.QTL)
  
  class(QeffRes) <- c("QeffRes", "list")
  
  return(QeffRes)
  
}