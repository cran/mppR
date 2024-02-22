#############
# mppGE_CIM #
#############

#' MPP GxE Composite Interval Mapping
#' 
#' Computes multi-QTL models with cofactors along the genome using an approximate
#' mixed model computation. An initial variance covariance (VCOV) structure is
#' calculated using function from the \code{nlme} package. Then, this information
#' is used to estimate the QTL global and within parental effect significance using a
#' Wald test.
#' 
#' @details
#' 
#' The estimated model is the following:
#' 
#' \eqn{\underline{y}_{icj} = E_{j} + C_{cj} + \sum_{qc=1}^{n_{cof}} x_{i_{qc}p} + \beta_{pj} + x_{i_{q}p} * \beta_{pj} + \underline{GE}_{icj} + \underline{e}_{icj}}
#'
#' For further details see the vignette.
#' 
#' It is possible to calculate one initial VCOV using a null model with all
#' the cofactors (\code{VCOV_data = "unique"}) or one VCOV per combination of
#' cofactors (\code{VCOV_data = "minus_cof"}). In the later case, the cofactor
#' that fall witin a distance of \code{window} on the left and right of a QTL
#' position is removed for the calculation of the initial VCOV. Therefore,
#' N_cof + 1 VCOV are calculated.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits (environments) should be used.
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. 'CS' for compound symmetry assuming a unique
#' genetic covariance between environments. 'CSE' for cross-specific within
#' environment error term. 'CS_CSE' for both compound symmetry plus
#' cross-specific within environment error term. 'UN' for unstructured
#' environmental variance covariance structure allowing a specific genotypic
#' covariance for each pair of environments. Default = 'UN'
#' 
#' @param VCOV_data \code{Character} specifying if the reference VCOV should
#' be formed  taking all cofactors into consideration ("unique") or if different
#' VCOVs should be formed by removing the cofactor information that is too close
#' of a tested QTL position ("minus_cof"). Default = "unique"
#' 
#' @param cofactors Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function \code{QTL_select()} or
#' a vector of \code{character} marker positions names.
#' Default = NULL.
#' 
#' @param cof_red \code{Logical} value specifying if the cofactor matrix should
#' be reduced by only keeping the significant allele by environment interaction.
#' Default = FALSE
#' 
#' @param cof_pval_sign \code{Numeric} value specifying the p-value significance
#' of an allele by environment term to be kept in the model. Default = 0.1 
#'
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 20.
#' 
#' @param ref_par Optional \code{Character} expression defining the parental
#' allele that will be used as reference for the parental model. Default = NULL
#'
#' @param n.cores \code{Numeric}. Specify here the number of cores you like to
#' use. Default = 1.
#' 
#' @param maxIter maximum number of iterations for the lme optimization algorithm.
#' Default = 100.
#' 
#' @param msMaxIter maximum number of iterations for the optimization step inside
#' the lme optimization. Default = 100.
#'
#' @return Return:
#'
#' \item{CIM }{\code{Data.frame} of class \code{QTLprof}. with five columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) integer position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-val) of the global QTL effect
#' across environments 6) p-values of the within environment QTL effects
#' (one column per environment); and p-values of the within environment parental
#' QTL allelic effects (one column per parent environment combination).}
#'
#' @author Vincent Garin
#' 
#' @references
#' 
#' Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2021). nlme: Linear
#' and Nonlinear Mixed Effects Models_. R package version 3.1-152,
#' <URL: https://CRAN.R-project.org/package=nlme>.
#'
#' \code{\link{mppGE_SIM}},
#' \code{\link{mppGE_proc}}
#'
#' @examples
#'
#' data(mppData_GE)
#' 
#' cofactors <- mppData_GE$map$mk.names[c(35, 61)]
#'
#' CIM <- mppGE_CIM(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'),
#'                      cofactors = cofactors, window = 20)
#'                      
#' Qpos <- QTL_select(CIM)
#'                       
#' plot(CIM)
#'
#' plot_allele_eff_GE(mppData = mppData_GE, nEnv = 2, EnvNames = c('CIAM', 'TUM'),
#'                    Qprof = CIM, Q.eff = 'par', QTL = Qpos, text.size = 14)
#'
#' @export
#'


mppGE_CIM <- function(mppData, trait, VCOV = 'UN', VCOV_data = "unique",
                      cofactors = NULL, cof_red = FALSE, cof_pval_sign = 0.1,
                      window = 20, ref_par = NULL, n.cores = 1, maxIter = 100,
                      msMaxIter = 100)
{
  
  #### 1. Check data format and arguments ####
  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = "par", VCOV = VCOV,
                  QTL_ch = FALSE, fast = TRUE, CIM = TRUE,
                  cofactors = cofactors, ref_par = ref_par)
  
  #### 2. Form required elements for the analysis ####
  nEnv <- length(trait)
  TraitEnv <- c(mppData$pheno[, trait])
  NA_id <- is.na(TraitEnv)
  vect.pos <- 1:dim(mppData$map)[1]
  n_cof <- nrow(cofactors)
  
  #### 3. Cofactor matrices ####
  
  # function for cofactor partition 
  test.cof <- function(x, map, window) {
    
    t1 <- map$chr == as.numeric(x[1])
    t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
    !(t1 & t2)
    
  }
  
  if(is.character(cofactors)){
    
    cof.pos <- which(mppData$map$mk.names %in% cofactors)
    cof_inf <- mppData$map[mppData$map$mk.names %in% cofactors, ]
    cof.part <- apply(X = cof_inf[, c(2, 4)], MARGIN = 1, FUN = test.cof,
                      map = mppData$map, window = window)
    
  } else {
    
    cof.pos <- which(mppData$map$mk.names %in% cofactors$mk.names)
    cof.part <- apply(X = cofactors[, c(2, 4)], MARGIN = 1, FUN = test.cof,
                      map = mppData$map, window = window)
    
  }
  
  cof_list <- mapply(FUN = inc_mat_QTL, x = cof.pos,
                     MoreArgs = list(Q.eff = "par", mppData = mppData),
                     SIMPLIFY = FALSE)
  cof_list <- lapply(cof_list, function(x) x[, -1])
  
  ind.ref <- paste0('c', 1:dim(cof.part)[2])
  
  cof.comb <- apply(X = cof.part, MARGIN = 1,
                    FUN = function(x, ind.ref) paste(ind.ref[x], collapse = ""),
                    ind.ref = ind.ref)
  
  cof.comb[cof.comb == ""] <- "c999"
  
  # cofactor combination
  unique.cof.comb <- unique(cof.comb)
  n_cof_comb <- length(unique.cof.comb)
  cof_id_list <- strsplit(x = unique.cof.comb, split = 'c')
  cof_id_list <- lapply(X = cof_id_list, function(x) x[-1])
  
  # combined cofactor matrices
  cof_mat_list <- vector(mode = 'list', length = n_cof_comb)
  
  for(i in 1:n_cof_comb){
    
    cof_i <- as.numeric(cof_id_list[[i]])
    if(!(999 %in% cof_i)){
      cof_mat_list[[i]] <- do.call(cbind, cof_list[cof_i])
    } 
    
  }
  
  names(cof_mat_list) <- unique.cof.comb
  
  #### 4. Calculate the VCOV matrices ####
  
  # Calculate one VCOV for each combination of cofactors ...
  if(VCOV_data == 'minus_cof'){
    
    VCOV_list <- vector(mode = 'list', length = n_cof_comb)
    
    for(i in 1:n_cof_comb){
      
      m <- MM_comp(mppData = mppData, nEnv = nEnv, y = TraitEnv,
                   cof_mat = cof_mat_list[[i]], VCOV = VCOV,
                   maxIter = maxIter, msMaxIter = msMaxIter)
      
      VCOV_list[[i]] <- getVCOV(mppData = mppData, model = m$model, VCOV = VCOV,
                                data = m$data, nEnv = nEnv, inv = TRUE)
      
    }
    
    names(VCOV_list) <- unique.cof.comb
    
    # ... or calculate a unique VCOV with all cofactors.  
  } else {
    
    s_cof_comb <- unique.cof.comb[which.max(nchar(unique.cof.comb))]
    
    m <- MM_comp(mppData = mppData, nEnv = nEnv, y = TraitEnv,
                 cof_mat = cof_mat_list[[s_cof_comb]], VCOV = VCOV,
                 maxIter = maxIter, msMaxIter = msMaxIter)
    
    VCOV_m <- getVCOV(mppData = mppData, model = m$model, VCOV = VCOV,
                      data = m$data, nEnv = nEnv, inv = TRUE)
    
  }
  
  ##### 5. Calculate the QTL effects #####
  
  if(n.cores > 1){ parallel <- TRUE; cluster <- makeCluster(n.cores)
  } else { parallel <- FALSE; cluster <- NULL }
  
  nGeno <- dim(mppData$pheno)[1]
  env <- rep(paste0('E', 1:nEnv), each = nGeno)
  cross <- rep(mppData$cross.ind, nEnv)
  geno <- rep(rownames(mppData$pheno), nEnv)
  cross_env <- factor(paste0(cross, '_', env))
  
  cross_mat <- model.matrix(~ cross_env -1)
  cross_mat <- cross_mat[!NA_id, ]
  
  log_pval_tot <- c()
  
  for(c in 1:n_cof_comb){
    
    sel_cof_comb <- unique.cof.comb[c]
    vect.pos <- which(cof.comb %in% sel_cof_comb)
    
    if(VCOV_data == "minus_cof"){ # selection of VCOV calculated minus cof pos ...
      
      Vi <- VCOV_list[[sel_cof_comb]]
      
    } else { # ... or global unique VCOV
      
      Vi <- VCOV_m
      
    }
    
    # form the cofactor matrices
    cof_mat <-  cof_mat_list[[sel_cof_comb]]
    
    if(!is.null(cof_mat)){
      cof_mat <- diag(nEnv) %x% cof_mat
      cof_mat <- cof_mat[!NA_id, ]
      colnames(cof_mat) <- paste0('cof', 1:ncol(cof_mat))
      
      if(cof_red){
        
        cof_mat <- cof_mat_reduce(y = TraitEnv[!NA_id], cross_mat = cross_mat,
                                    cof_mat = cof_mat, cof_pval_sign = cof_pval_sign)
        
      }
      
    }
    
    #### possibility of PCA reduction
    
    if (parallel) {
      
      log.pval <- parLapply(cl = cluster, X = vect.pos, fun = W_QTL,
                            y = TraitEnv[!NA_id], Vi = Vi,
                            mppData = mppData,  nEnv = nEnv, 
                            Q.eff = "par", cross_mat = cross_mat,
                            cof_mat = cof_mat, NA_id = NA_id, ref_par = ref_par)
      
    } else {
      
      log.pval <- lapply(X = vect.pos, FUN = W_QTL,
                         y = TraitEnv[!NA_id], Vi = Vi,
                         mppData = mppData,  nEnv = nEnv, 
                         Q.eff = "par", cross_mat = cross_mat,
                         cof_mat = cof_mat, NA_id = NA_id, ref_par = ref_par)
      
    }
    
    log.pval <- t(data.frame(log.pval))
    log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
    log.pval[is.na(log.pval[, 1]), 1] <- 0
    log.pval <- cbind(log.pval, vect.pos)
    
    log_pval_tot <- rbind(log_pval_tot, log.pval)
    
  }
  
  log.pval <- log_pval_tot[order(log_pval_tot[, ncol(log_pval_tot)]), ]
  log.pval <- log.pval[, -ncol(log.pval)]
  
  if(n.cores > 1){stopCluster(cluster)}
  
  #### 6. format the results and return ####
  
  CIM <- Qprof_process(mppData = mppData, Q.eff = "par", log.pval = log.pval,
                       nEnv = nEnv)
  
  return(CIM)
  
}