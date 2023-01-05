#######################
# QTL_effect_main_QxE #
#######################

#' Estimation of QTL main effect and QTLxE effect
#'
#'
#' Decomposition of the QTL effect into main component across environments and
#' QTLxE component.
#' 
#' The function estimate the QTL parent allele main effect across environments
#' as well the QTLxE effect. The significance of the QTL parental main effect
#' as well as the QTLxE effect are also estimated and returned as -log10(p-value).
#' 
#' The function use two models, one where the QTL parent allele effect are
#' considered to be different in each environments (QTLxE model) and a model
#' where the QTL parental effect are assumed to be constant across environment
#' (QTL main model). Concerning the model to estimate the QTL main effect, there
#' are two option, the first (default) option (QmainQi = TRUE), estimate a model
#' where only the ith QTL is defined with a main effect and the other position
#' are assumed to have parental effect that vary in each environment (same as
#' the QTLxE model). In that case, the function estimate as many QTL main
#' model as there are QTL positions to get the main effect estimate of each
#' QTL position. The alternative option (QmainQi = FALSE), calculate a single
#' model where all QTL are defined with a main effect term. The estimated
#' main effect obtained with the two options are generally very similar. The
#' second option is less time consumming.
#' 
#' The QTL main allelic effect is the deviation of the parental allelic effect
#' with respect to the reference parent (e.g. the central or recurrent parent
#' in a NAM population)
#' 
#' The estimation is performed using an exact mixed model with function from R
#' package \code{nlme}. The significance of the allele effect is assessed using a 
#' Wald test.
#' 
#' @param mppData An object of class \code{mppData}.
#'
#' @param trait \code{Character vector} specifying which traits (environments) should be used.
#' 
#' @param env_id \code{Character} vector specifying the environment names.
#' By default, E1, ... En
#' 
#' @param VCOV VCOV \code{Character} expression defining the type of variance
#' covariance structure used. 'CS' for compound symmetry assuming a unique
#' genetic covariance between environments. 'CSE' for cross-specific within
#' environment error term. 'CS_CSE' for both compound symmetry plus
#' cross-specific within environment error term. 'UN' for unstructured
#' environmental variance covariance structure allowing a specific genotypic
#' covariance for each pair of environments. Default = 'UN'
#' 
#' @param ref_par Optional \code{Character} expression defining the parental
#' allele that will be used as reference for the parental model. Default = NULL
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected marker positions obtained with the function QTL_select() or
#' a vector of \code{character} marker positions names. Default = NULL.
#' 
#' @param QmainQi \code{logical} value specifying how the QTL parental allele
#' main effects are estimated. For further explanation see the details section.
#' Default = TRUE
#'
#' @param maxIter maximum number of iterations for the lme optimization algorithm.
#' Default = 100.
#' 
#' @param msMaxIter maximum number of iterations for the optimization step inside
#' the lme optimization. Default = 100.
#'
#' @return Return:
#'
#' \code{List} with one \code{data.frame} per QTL that contains the following
#' elements:
#' 
#' \enumerate{
#' 
#' \item{QTL parent allele main effect expressed as deviation with respect to
#' the reference parent}
#' \item{QTL parent allele effect in environment j expressed as deviation with
#' respect to the reference parent}
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
#' Qeff <- QTL_effect_main_QxE(mppData = mppData_GE,
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

QTL_effect_main_QxE <- function(mppData, trait, env_id = NULL, VCOV = "UN",
                                ref_par = NULL, QTL = NULL, QmainQi = TRUE,
                                maxIter = 100, msMaxIter = 100){
  
  #### 1. Check data format and arguments ####
  check_mod_mppGE(mppData = mppData, trait = trait, Q.eff = "par", VCOV = VCOV,
                  QTL_ch = TRUE, QTL = QTL, fast = TRUE, CIM = FALSE,
                  ref_par = ref_par)
  
  if(!is.null(env_id)){
    
    if(!is.character(env_id)){
      
      stop('env_id is not a character vector.')
      
    }
    
    if(length(env_id) != length(trait)){
      
      stop('the length of env_id is not equal to the number of environments specified in trait.')
      
    }
    
  }
  
  #### 2. Form required elements for the analysis ####
  nEnv <- length(trait)
  TraitEnv <- c(mppData$pheno[, trait])
  NA_id <- is.na(TraitEnv)
  
  #### 3. QTL matrices ####
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
  
  # modify the names
  for(i in 1:length(QTL_list)){
    
    colnames(QTL_list[[i]]) <- mdf_par_name(nm = colnames(QTL_list[[i]]))
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
    QTL_nm_i <- paste0('QTL', i, '_E', rep(1:nEnv, each = n_allele_i))
    QTL_nm_i <- paste0(QTL_nm_i, '_', rep(colnames(QTL_list[[i]]), nEnv))
    QTL_nm <- c(QTL_nm, QTL_nm_i)
    QTL_mat_QEI[[i]] <- QTL_i_env
    
  }
  
  QTL_mat_QEI <- do.call(cbind, QTL_mat_QEI)
  colnames(QTL_mat_QEI) <- QTL_nm
  
  #### 4. General element to form the model ####
  
  nGeno <- dim(mppData$pheno)[1]
  env <- rep(paste0('E', 1:nEnv), each = nGeno)
  cross <- rep(mppData$cross.ind, nEnv)
  geno <- rep(rownames(mppData$pheno), nEnv)
  cross_env <- paste0(cross, '_', env)
  
  d <- data.frame(trait = TraitEnv, env = env, cross_env = cross_env, geno = geno)
  d[, 2:4] <- lapply(d[, 2:4], as.factor)
  
  Qterm_main <-  colnames(QTL_mat_main)
  Qterm_GxE <- colnames(QTL_mat_QEI)
  
  #### 5. Computation of the mixed model (main) ####
  
  if(QmainQi){ # Only QTLi fitted as main the rest as QxE
    
    d_m <- data.frame(d, QTL_mat_main, QTL_mat_QEI)
    
    Qeff_main <- vector(mode = 'list', length = nQTL)
    
    for(i in 1:nQTL){
      
      # remove the QxE effect for the ith QTL
      d_m_i <- d_m
      QTL_id_m <- grepl(pattern = paste0('QTL', i, '_main'), x = colnames(d_m_i))
      QTL_id <- grepl(pattern = paste0('QTL', i), x = colnames(d_m_i))
      d_m_i <- d_m_i[, !(QTL_id & !QTL_id_m)]
      
      d_m_i <- remove_singularities(d_m_i)
      
      Q_id <- colnames(d_m_i)[5:ncol(d_m_i)]
      fix_form <- paste0('trait~-1 + cross_env+', paste(Q_id, collapse = '+'))
      m <- lme_comp(fix_form = fix_form, VCOV = VCOV, data = d_m_i,
                    maxIter = maxIter, msMaxIter = msMaxIter)
      
      Qeff_main[[i]] <- W_test_Qpar_main(m = m, nQTL = nQTL)[[i]]
      
    }
    
  } else { # all QTL terms fitted as main effect
    
    d_m <- data.frame(d, QTL_mat_main)
    d_m <- remove_singularities(d_m)
    
    Q_id <- colnames(d_m)[5:ncol(d_m)]
    fix_form <- paste0('trait~-1 + cross_env+', paste(Q_id, collapse = '+'))
    m <- lme_comp(fix_form = fix_form, VCOV = VCOV, data = d_m,
                  maxIter = maxIter, msMaxIter = msMaxIter)
    
    Qeff_main <- W_test_Qpar_main(m = m, nQTL = nQTL)
    
  }
  
  #### 6. Computation of the mixed model (QTLxE) ####
  
  d_m <- data.frame(d, QTL_mat_QEI)
  d_m <- remove_singularities(d_m)
  
  Q_id <- colnames(d_m)[5:ncol(d_m)]
  fix_form <- paste0('trait~-1 + cross_env+', paste(Q_id, collapse = '+'))
  m <- lme_comp(fix_form = fix_form, VCOV = VCOV, data = d_m,
                maxIter = maxIter, msMaxIter = msMaxIter)
  
  Qeff_QxE <- W_test_Qpar_GxE(m = m, nQTL = nQTL, nEnv = nEnv, env_id = env_id)
  
  
  #### 7. Results processing ####
  
  Q_res <- vector(mode = 'list', length = nQTL)
  
  par_ref <- mdf_par_name(nm = mppData$parents)
  
  d_QTL_ref <- data.frame(par = par_ref)
  
  for(i in 1:nQTL){
    
    Qmain_i <- Qeff_main[[i]]
    QTLxE_i <- Qeff_QxE[[i]]
    
    # put the two dataset in the same (parent) order
    rownames(QTLxE_i) <- QTLxE_i$par
    QTLxE_i <- QTLxE_i[Qmain_i$par, ]
    
    B_QxE <- QTLxE_i[, 2:(1+nEnv)]
    
    Q_i <- data.frame(par = Qmain_i$par, Effect_main = Qmain_i$Effect, B_QxE,
                      logP_main = Qmain_i$log10P,
                      logP_QxE = QTLxE_i$log10P)
    
    d_QTL_i <- merge(x = d_QTL_ref, y = Q_i, by = 'par', all.x = TRUE)
    rownames(d_QTL_i) <- d_QTL_i$par
    d_QTL_i <- d_QTL_i[par_ref, ]
    d_QTL_i <- d_QTL_i[, -1]
    
    Q_res[[i]] <- round(d_QTL_i, 3)
    
  }
  
  names(Q_res) <- paste0('QTL', 1:nQTL)
  
  return(Q_res)
  
}