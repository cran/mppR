##################
# MQE_gen_effects #
##################

#' QTL genetic effects multi-QTL effect model
#' 
#' Compute a multi-QTL model with a list of QTL candidates (\code{QTL}) and return
#' the decomposed QTL genetic effects per cross or per parents. The list of QTL
#' can be of different types (cross-specific, parental, ancestral or bi-allelic).
#' The type of QTL effects are specified in the vector \code{Q.eff}.
#' 
#' This function computes for each QTL position the genetic effects of the
#' cross, parental, ancestral or SNP allele components. For the cross-specific
#' model (\code{Q.eff = "cr"}), the genetics effects represent the substitution
#' effect of an single allele from the parent 2 (or B) with respect to an allele
#' coming from the parent 1 or A. All effects are given in absolute value with
#' the parent that cary the positive allele.
#' 
#' For the parental and the ancestral model (\code{Q.eff = "par" or "anc"}), the
#' reference allele is defined per interconneted part. The most frequent
#' parental (ancestral) allele is set as reference. Effects of the other alleles
#' are estimated as deviation with respect to the reference. For more details
#' about reference definition see \code{\link{QTL_gen_effects}} and
#' \code{\link{design_connectivity}}.
#' 
#' For the bi-allelic model (\code{Q.eff = "biall"}), the genetic effects
#' represent the effects of a single allele copy of the least frequent allele.
#'
#' @param mppData An object of class \code{mppData}
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#' 
#' @param QTL Vector of \code{character} markers positions
#' names. Default = NULL.
#' 
#' @param Q.eff \code{Character} vector indicating for each QTL position the
#' type of QTL effect among: "cr", "par", "anc" and "biall". For details look at
#' \code{\link{mpp_SIM}}.
#' 
#' @param ref.par Optional \code{Character} expression defining the parental
#' allele that will be used as reference for the parental model. For the
#' ancestral model, the ancestral class containing the reference parent will be
#' set as reference. \strong{This option can only be used if the MPP design is
#' composed of a unique connected part}. Default = NULL.
#'
#' @return Return:
#'
#' \item{results}{\code{List} of \code{data.frame} (one per QTL) containing the
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
#'   }
#' 
#' }
#' 
#'   
#' @author Vincent Garin
#' 
#' @seealso \code{\link{MQE_proc}}, \code{\link{mpp_SIM}},
#' \code{\link{QTL_gen_effects}}, \code{\link{design_connectivity}}
#'
#' @examples
#'
#' data(mppData)
#' 
#' SIM <- mpp_SIM(mppData = mppData)
#' QTL <- QTL_select(SIM)
#' 
#' QTL.eff <- MQE_gen_effects(mppData = mppData, QTL = QTL[, 1],
#'                           Q.eff = c("anc", "par", "biall"))
#'
#' @export
#'


MQE_gen_effects <- function(mppData = NULL, trait = 1, QTL = NULL, Q.eff,
                           ref.par = NULL){
  
  # 1. check the data format
  ##########################
  
  check.MQE(mppData = mppData, trait = trait, Q.eff = Q.eff,
            VCOV = "h.err", QTL = QTL, fct = "QTLeffects")
  
  if(!is.null(ref.par)){
    
   n_con_part <- length(design_connectivity(mppData$par.per.cross, plot_des = FALSE))
   
   if(n_con_part > 1){
     
     stop('ref.par option can only be used with design contaning a single interconnected part. See design_connectivity()')
   }
    
    if(!(ref.par %in% mppData$parents)){
      
      mes <- paste('ref.par must be one of:', paste(mppData$parents, collapse = ', '))
      
      stop(mes)
      
    }
    
  }
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  ### 2.3 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.4 Formation of the list of QTL
  
  # order list of QTL positions
  
  Q.pos <- vapply(X = QTL,
                  FUN = function(x, mppData) which(mppData$map[, 1] == x),
                  FUN.VALUE = numeric(1), mppData = mppData)
  
  Q.ord <- data.frame(QTL, Q.eff, Q.pos, stringsAsFactors = FALSE)
  
  Q.ord <- Q.ord[order(Q.pos), ]
  
  QTL <- Q.ord[, 1]; Q.eff <- Q.ord[, 2]
  
  Q.pos <- which(mppData$map[, 1] %in% QTL)
  
  # form a list of QTL incidence matrices with different type of QTL effect.
  
  Q.list <- mapply(FUN = inc_mat_QTL, x = Q.pos, Q.eff = Q.eff,
                   MoreArgs = list(mppData = mppData, order.MAF = TRUE),
                   SIMPLIFY = FALSE)
  
  # define the vector for reference parent
  
  n.QTL <- length(Q.list)
  
  if(!is.null(ref.par)){
    
    ref_par_l <- vector(mode = 'list', length = n.QTL)
    ref_par_l[Q.eff %in% c('par', 'anc')] <- ref.par
    
  } else {ref_par_l <- vector(mode = 'list', length = n.QTL)}
  
  # order the QTL incidence matrices
  
  order.Qmat <- mapply(FUN = IncMat_QTL_MAF, QTL = Q.list,
                       Q.eff_i = Q.eff, Q.pos_i = Q.pos, ref.par = ref_par_l,
                       MoreArgs = list(mppData = mppData),
                       SIMPLIFY = FALSE)
  
  Q.list <- lapply(X = order.Qmat, FUN = function(x) x$QTL)
  n.QTL <- length(Q.list)
  names(Q.list) <- paste0("Q", 1:n.QTL)
  allele_order <- lapply(X = order.Qmat, FUN = function(x) x$allele_order)
  con.ind <- lapply(X = order.Qmat, FUN = function(x) x$con.ind)
  
  n.allele <- lapply(X = Q.list, function(x) dim(x)[2])
  Q.ind <- rep(paste0("Q", 1:n.QTL), n.allele)
  
  
  # 3. Compute the model
  ######################
  
  model <- QTLModelQeff(mppData = mppData, trait = t_val, cross.mat = cross.mat,
                        Q.list = Q.list, VCOV = 'h.err')
  
  
  # 4. Results processing
  #######################
  
  results <- summary(model)$coefficients
  index <- (substr(rownames(results), 1, 1) == "Q")
  results <- subset(x = results, subset = index, drop = FALSE)
  
  ref.names <- names(coef(model))[-c(1:mppData$n.cr)]
  
  ref.mat <- matrix(rep(c(0, 0, 0, 1), length(ref.names)),
                    nrow = length(ref.names), byrow = TRUE)
  
  index <- match(rownames(results), ref.names)
  ref.mat[index, ] <- results
  rownames(ref.mat) <- ref.names
  
  # add sign stars
  
  Sign <- sapply(ref.mat[, 4], FUN = sign.star)
  results <- data.frame(ref.mat, Sign, stringsAsFactors = FALSE)
  
  # split the results per QTLs
  
  Q.res <- split(x = results,
                 f = factor(Q.ind, levels = paste0("Q", 1:n.QTL)))
  
  
  # 4.2 results processing
  
  results <- mapply(Qeff_res_processing_MQE, Q.res = Q.res, Q.eff = Q.eff,
                    Q.pos = Q.pos, con.ind = con.ind, allele_order = allele_order,
                    Q.nb = 1:n.QTL,
                    MoreArgs = list(mppData = mppData, VCOV = 'h.err'),
                    SIMPLIFY = FALSE)
  
  names(results) <- paste0("Q", 1:length(results))
  
  return(results)
  
}