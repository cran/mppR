#######################
# QTL_effect_main_QEI #
#######################

#' Main and QTL by environment interaction model
#'
#' The function estimate a QTL model where each parental QTL allelic effect is
#' decomposed into a main effect and a QTL by environment effect (QEI). It allows
#' the user to determine which parental allelic effects have a significant
#' interaction with the environment. 
#' 
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits (environments) should be used.
#' 
#' @param env_id \code{Character} vector specifying the environment names.
#' By default, E1, ... En
#' 
#' @param ref_env Optional \code{Character} expression defining the environment
#' that will be used as reference for the parental model. By default, the last
#' environment is set as reference.
#' 
#' @param ref_par Optional \code{Character} expression defining the parental
#' allele that will be used as reference for the parental model. Default = NULL
#'
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. 'CS' for compound symmetry assuming a unique
#' genetic covariance between environments. 'CSE' for cross-specific within
#' environment error term. 'CS_CSE' for both compound symmetry plus
#' cross-specific within environment error term. 'UN' for unstructured
#' environmental variance covariance structure allowing a specific genotypic
#' covariance for each pair of environments. Default = 'UN'
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function QTL_select() or
#' a vector of \code{character} marker positions names. Default = NULL.
#'
#' @param maxIter maximum number of iterations for the lme optimization algorithm.
#' Default = 100.
#' 
#' @param msMaxIter maximum number of iterations for the optimization step inside
#' the lme optimization. Default = 100.
#' 
#' @details
#' 
#' The function estimate the following model
#' 
#' \eqn{y_{icj} = E_j + C_{c_j} + \sum_{q=1}^{n_{QTL}}{x_{i_{q}p}*(\alpha_{p} + \beta_{pj})} + GE_{ijc} + e_{ijc}}
#' 
#' where the QTL effect is decomposed into \eqn{\alpha_{p}} that represent the
#' main parental allelic effect across environments and \eqn{\beta_{pj}} which is
#' the QEI effect. allelic effects must be interpreted as deviation with respect
#' to the reference parent ('ref_par') in the reference environment ('ref_env').
#' By default the reference parent is the one with the highest allelic frequency
#' (e.g. central parent in a NAM population).
#' 
#' The estimation is performed using an exact mixed model with function from R
#' package \code{nlme}. The significance of the allele effect is assessed using a 
#' Wald test.
#'
#' @return Return:
#'
#' \code{List} with one \code{data.frame} per QTL that contains the following
#' elements:
#' 
#' \enumerate{
#' 
#' \item{To be filled}
#' 
#' \item{To be filled}
#' 
#' \item{Significance of the parent main effect expressed as the -log10(p-val)}
#' 
#' \item{Significance of the parent QTLxE effect expressed as the -log10(p-val)}
#' 
#' }
#'
#' @author Vincent Garin
#' 
#' @references
#' 
#' Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2021). nlme: Linear
#' and Nonlinear Mixed Effects Models_. R package version 3.1-152,
#' <URL: https://CRAN.R-project.org/package=nlme>.
#'
#' @examples
#' 
#' \dontrun{
#' 
#' data(mppData_GE)
#'
#' Qpos <- c("PZE.105068880", "PZE.106098900")
#'
#' Qeff <- QTL_effect_main_QEI(mppData = mppData_GE,
#'                             trait = c('DMY_CIAM', 'DMY_TUM', 'DMY_INRA_P', 'DMY_KWS'),
#'                             env_id = c('CIAM', 'TUM', 'INRA', 'KWS'),
#'                             QTL = Qpos)
#'
#' Qeff
#' 
#' }
#' 
#' @export
#'

QTL_effect_main_QEI <- function(mppData, trait, env_id = NULL, ref_env = NULL,
                                ref_par = NULL, VCOV = "UN", QTL = NULL,
                                maxIter = 100, msMaxIter = 100){
  
  #### 1. Check data format and arguments ----
  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = "par", VCOV = VCOV,
                  QTL_ch = TRUE, QTL = QTL, fast = TRUE, CIM = FALSE)
  
  if(!is.null(ref_env) & is.null(env_id)){
    stop("To identify the reference environment, you must provide the vector of environment name 'env_id'.")
  }
  
  #### 2. trait vector ----
  nEnv <- length(trait)
  if(is.null(env_id)){ env_id <- paste0("E", 1:nEnv)}
  TraitEnv <- c(mppData$pheno[, trait])
  NA_id <- is.na(TraitEnv)
  
  #### 3. QTL matrices ----
  if(is.character(QTL)){
    
    QTL.pos <- which(mppData$map$mk.names %in% QTL)
    
  } else {
    
    QTL.pos <- which(mppData$map$mk.names %in% QTL$mk.names)
    
  }
  
  nQTL <- length(QTL.pos)
  
  QTL_list <- mapply(FUN = inc_mat_QTL, x = QTL.pos,
                     MoreArgs = list(Q.eff = "par", mppData = mppData,
                                     order.MAF = TRUE, ref_par = ref_par),
                     SIMPLIFY = FALSE)
  
  QTL_list <- lapply(QTL_list, function(x) x[, -ncol(x)])
  
  nAllele <- sapply(QTL_list, function(x) ncol(x))
  
  # modify the names and define ref parents
  ref_pars <- rep(NA, nQTL)
  p_nm <- mdf_par_name(mppData$parents)
  
  for(i in 1:nQTL){
    
    colnames(QTL_list[[i]]) <- mdf_par_name(nm = colnames(QTL_list[[i]]))
    ref_pars[i] <- p_nm[-which(p_nm %in% colnames(QTL_list[[i]]))]
    
  }
  
  # combined QTL matrices
  QTL_mat <- do.call(cbind, QTL_list)
  
  # QTL main effect matrix
  QTL_mat_main <- matrix(1, nEnv, 1) %x% QTL_mat 
  
  Q_nm <- colnames(QTL_mat)
  Q_id <- paste0('QTL', rep(1:nQTL, nAllele))
  QTL_nm <- paste0(Q_id, '_main_', Q_nm)
  colnames(QTL_mat_main) <- QTL_nm
  
  # QEI matrix
  QTL_mat_QEI <- vector(mode = 'list', length = nQTL)
  QTL_nm <- c()
  
  for(i in 1:nQTL){
    
    QTL_i_env <- diag(nEnv) %x% QTL_list[[i]]
    n_allele_i <- ncol(QTL_list[[i]])
    QTL_nm_i <- paste0('QTL', i, '_', rep(env_id, each = n_allele_i))
    QTL_nm_i <- paste0(QTL_nm_i, '_', rep(colnames(QTL_list[[i]]), nEnv))
    
    if(!is.null(ref_env)){
      env_nb <- which(env_id == ref_env)
      st_pt <- seq(1, length(QTL_nm_i), by = n_allele_i)[env_nb]
      rem_pos <- st_pt:(st_pt + n_allele_i - 1)
      QTL_nm_i <- QTL_nm_i[-rem_pos]
      QTL_i_env <- QTL_i_env[, -rem_pos]
    }
    
    QTL_nm <- c(QTL_nm, QTL_nm_i)
    QTL_mat_QEI[[i]] <- QTL_i_env
    
  }
  
  QTL_mat_QEI <- do.call(cbind, QTL_mat_QEI)
  colnames(QTL_mat_QEI) <- QTL_nm
  
  #### 4. General element to form the model ----
  
  nGeno <- dim(mppData$pheno)[1]
  env <- rep(paste0('E', 1:nEnv), each = nGeno)
  cross <- rep(mppData$cross.ind, nEnv)
  geno <- rep(rownames(mppData$pheno), nEnv)
  cross_env <- paste0(cross, '_', env)
  
  d <- data.frame(trait = TraitEnv, env = env, cross_env = cross_env, geno = geno)
  d[, 2:4] <- lapply(d[, 2:4], as.factor)
  
  #### 5. Computation of the mixed model ----
  d_m <- data.frame(d, QTL_mat_main, QTL_mat_QEI)
  d_m <- remove_singularities(d_m)
  
  Q_id <- colnames(d_m)[5:ncol(d_m)]
  fix_form <- paste0('trait~-1 + cross_env+', paste(Q_id, collapse = '+'))
  
  m <- lme_comp(fix_form = fix_form, VCOV = VCOV, data = d_m,
                maxIter = maxIter, msMaxIter = msMaxIter)
  
  #### 6. estimation of the Wald statistic for the main and QEI component for each parent ----
  
  Q_sign <- W_test_Qmain_QEI(m = m, nQTL = nQTL, n_par = length(mppData$parents),
                             ref_par = ref_pars)
  
  names(Q_sign) <- paste0('QTL', 1:nQTL)
  
  #### 7. Storage of the estimated QTL effects with significance ----
  Beta <- m$coefficients$fixed
  Q_ind <- !grepl(pattern = 'cross_env', x = names(Beta))
  B_QTL <- Beta[Q_ind]
  B_SD <- sqrt(diag(m$varFix)[Q_ind])
  
  # Wald test
  W_Qa <- (B_QTL/B_SD)^2
  W_pval <- pchisq(W_Qa, 1, lower.tail = FALSE)
  W_sign <- sapply(W_pval, sign.star) 
  
  # process the results
  res_tab <- data.frame(Effect = round(B_QTL, 3), Std.dev = round(B_SD, 3),
                        Wald = round(W_Qa, 2), df = 1, p.val = W_pval,
                        sign = W_sign)
  
  p_nm <- mdf_par_name(mppData$parents)
  
  # nEnv + 1: all env (one ref) + main effect 
  ref_QTL <- rep(paste0('QTL', 1:nQTL), each = (mppData$n.par * (nEnv + 1)))
  ref_QTL1 <- paste0(ref_QTL, '_', rep(p_nm, (nEnv+1)))
  ref_QTL1 <- paste0(ref_QTL1, '_', rep(c('main', env_id), each = mppData$n.par))
  
  ref_QTL2 <- paste0(ref_QTL, '_', rep(c('main', env_id), each = mppData$n.par))
  ref_QTL2 <- paste0(ref_QTL2, '_', rep(p_nm, (nEnv + 1)))
  
  # merge data to fill with NA non-estimated effects
  ref_QTL <- data.frame(ref_id = ref_QTL2)
  res_tab <- data.frame(ref_id = rownames(res_tab), res_tab)
  res_tab <- merge(x = ref_QTL, y = res_tab, by = 'ref_id', all.x = TRUE)
  
  rownames(res_tab) <- res_tab$ref_id
  res_tab <- res_tab[ref_QTL$ref_id, ]
  res_tab <- res_tab[, -1]
  
  Q_f <- strsplit(x = rownames(res_tab), split = '_')
  Q_f <- unlist(lapply(Q_f, `[[`, 1))
  Q_f <- factor(Q_f, levels = unique(Q_f))
  Q_eff <- split(x = res_tab, f = Q_f)
  
  return(list(Q_sign = Q_sign, Q_eff = Q_eff))
  
}