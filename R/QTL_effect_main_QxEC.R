########################
# QTL_effect_main_QxEC #
########################

#' Estimation of a model with main and QTL by environmental sensitivity terms
#'
#' 
#' After estimating which parental allelic effects have a significant interaction
#' with the environment (QEI), the function extends the model for the allelic
#' effect with a significant QEI to characterize this interaction in terms of
#' sensitivity to (a) specific environmental covariate(s).
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
#' allele that will be used as reference for the parental model. By default,
#' the parent with the largest MAF is set as reference.
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
#' @param thre_QTL \code{Numerical} value specifying the -log10(p-val) threshold
#' for a parental QTL allele to be considered as significant. By default,
#' thre_QTL = 2, which correspond to a p-value of 0.01.
#' 
#' @param EC \code{Numeric} matrix with environments as row and environmental
#' covariates (EC) as column. The cell i, j of EC specify the value of the
#' jth EC in environment i.
#' 
#' @param Qmain_QEI results from \code{\link{QTL_effect_main_QEI}}
#'
#' @param maxIter maximum number of iterations for the lme optimization algorithm.
#' Default = 100.
#' 
#' @param msMaxIter maximum number of iterations for the optimization step inside
#' the lme optimization. Default = 100.
#'
#' @details
#' 
#' The function first estimate the parental QTL allele main and QTLxE effect
#' using the function \code{\link{QTL_effect_main_QEI}}. Optionally the output
#' of \code{\link{QTL_effect_main_QEI}} can be passed through the `Qmain_QEI`
#' argument. The function consider that a parental QTL allele significantly
#' interacts with the environment if its QTLxE term is significant at the
#' `thre_QTL` level. Thre_QTL is expressed in terms of -log10(p-val).
#' For example, for p-val = 0.01, thre_QTL = -log10(p-val) = 2. Given this
#' information, the effect of the parental QTL allele with a significant QEI
#' are extended like that \eqn{\beta_{pj} = EC_j*S_p+l_{p\epsilon}} where
#' \eqn{EC_j} represents the EC value in environment j associated with the
#' sensitivity term \eqn{S_p}. The \eqn{S_{p}} determines the rate of change of
#' the parental QTL allelic additive effect given an extra unit of EC. Finally,
#' \eqn{l_{p\epsilon}} is a residual effect. The fitted model becomes:
#' 
#' \eqn{\underline{y}_{icj} = E_{j} + C_{cj} + \sum_{q=1}^{n_{QTL}} x_{i_{q}p} (\alpha_p + \beta_{pj}) + x_{i_{q}pxE} (\alpha_p + EC_j*S_p+l_{p\epsilon}) + \underline{GE}_{icj} + \underline{e}_{icj}}
#' 
#' The estimation is performed using an exact mixed model with function from R
#' package \code{nlme}. The significance of \eqn{S_{p}} is assessed using a 
#' Wald test.
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
#' @seealso
#' 
#' \code{\link{QTL_effect_main_QEI}}
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
#' EC <- matrix(c(180, 310, 240, 280), 4, 1)
#' rownames(EC) <- c('CIAM', 'TUM', 'INRA', 'KWS')
#' colnames(EC) <- 'cum_rain'
#'
#' Qeff <- QTL_effect_main_QxEC(mppData = mppData_GE,
#'                          trait = c('DMY_CIAM', 'DMY_TUM', 'DMY_INRA_P', 'DMY_KWS'),
#'                          env_id = c('CIAM', 'TUM', 'INRA', 'KWS'),
#'                          QTL = Qpos, EC = EC)
#'
#' Qeff
#' 
#' }
#' 
#' @export
#'

QTL_effect_main_QxEC <- function(mppData, trait, env_id = NULL, ref_env = NULL,
                                 ref_par = NULL, VCOV = "UN", QTL = NULL,
                                 thre_QTL = 2, EC, Qmain_QEI = NULL,
                                 maxIter = 100, msMaxIter = 100){
  
  #### 1. Check some arguments ----
  if(!is.matrix(EC)){stop('EC is not a matrix')}
  if(!is.numeric(EC)){stop('EC is not (fully) numeric.')}
  
  if(nrow(EC) != length(trait)){
    stop('The number of EC row is not equal to the number of environments specified in trait.')
  }
  
  #### 2. Estimation of the QTL main and QxE effects ----
  
  if(!is.null(Qmain_QEI)){
    
    Qeff <- Qmain_QEI
    
  } else {
    
    Qeff <- QTL_effect_main_QEI(mppData = mppData, trait = trait,
                                env_id = env_id, ref_env = ref_env, VCOV = VCOV,
                                ref_par = ref_par, QTL = QTL, maxIter = maxIter,
                                msMaxIter = msMaxIter)
  }
  
  
  
  #### 3. Determine the significance of QTLxE ----
  
  Q_sign <- determine_QEI_sign(Q_sign = Qeff$Q_sign,
                               par_id = mdf_par_name(mppData$parents),
                               thre_QTL = thre_QTL)
  
  # Check if there is at least some QTLxE term significant
  QxE_sign_vect <- c(Q_sign$QTLxE)
  if(any(QxE_sign_vect[!is.na(QxE_sign_vect)])){
    
    #### 4. Form required elements for the analysis ----
    nEnv <- length(trait)
    TraitEnv <- c(mppData$pheno[, trait])
    NA_id <- is.na(TraitEnv)
    
    if(is.null(ref_env)){
      env_id <- paste0('E', 1:nEnv)
    }
    
    EC_id <- paste0('EC', 1:ncol(EC))
    EC_org_nm <- colnames(EC)
    n_EC <- length(EC_id)
    par_id <- mdf_par_name(mppData$parents)
    
    #### 5. QTL matrices ----
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
    
    #### 6. modify the (parent) names of QTL list and Q_sign info ----
    
    for(i in 1:length(QTL_list)){
      colnames(QTL_list[[i]]) <- mdf_par_name(colnames(QTL_list[[i]]))
    }
    
    #### 7. General element for MM ----
    
    nGeno <- dim(mppData$pheno)[1]
    env <- rep(paste0('E', 1:nEnv), each = nGeno)
    cross <- rep(mppData$cross.ind, nEnv)
    geno <- rep(rownames(mppData$pheno), nEnv)
    cross_env <- paste0(cross, '_', env)
    d <- data.frame(trait = TraitEnv, env = env, cross_env = cross_env, geno = geno)
    d[, 2:4] <- lapply(d[, 2:4], as.factor)  
    
    
    #### 8. QTL effect design matrices ----
    
    QTL_mat <- do.call(cbind, QTL_list)
    
    # QTL main effect matrix
    QTL_mat_main <- matrix(1, nEnv, 1) %x% QTL_mat 
    
    Q_nm <- colnames(QTL_mat)
    Q_id <- paste0('QTL', rep(1:nQTL, nAllele))
    QTL_nm <- paste0(Q_id, '_main_', Q_nm)
    colnames(QTL_mat_main) <- QTL_nm
    
    # QTLxE and QTLxEC matrix
    QTL_mat_EC <- form_QEI_QxEC_mat(QTL_list = QTL_list, EC = EC, Q_sign = Q_sign,
                                    env_id = env_id, ref_env = ref_env)
    
    
    #### 9. Computation of the mixed model ----
    d_m <- data.frame(d, QTL_mat_main, QTL_mat_EC)
    d_m <- remove_singularities(d_m)
    
    Q_id <- colnames(d_m)[5:ncol(d_m)]
    fix_form <- paste0('trait~-1 + cross_env+', paste(Q_id, collapse = '+'))
    
    m <- lme_comp(fix_form = fix_form, VCOV = VCOV, data = d_m,
                  maxIter = maxIter, msMaxIter = msMaxIter)
    
    #### 10. results processing ----
    Beta <- m$coefficients$fixed
    Q_ind <- grepl(pattern = 'QTL', x = names(Beta))
    B_QTL <- Beta[Q_ind]
    B_SD <- sqrt(diag(m$varFix)[Q_ind])
    
    # Wald test
    W_Qa <- (B_QTL/B_SD)^2
    W_pval <- pchisq(W_Qa, 1, lower.tail = FALSE)
    W_sign <- sapply(W_pval, sign.star) 
    
    res_tab <- data.frame(B_QTL, Std.dev = B_SD, p_val = W_pval,
                          logP = -log10(W_pval), sign = W_sign)
    
    name_dis <- strsplit(x = names(B_QTL), split = '_')
    res_tab$Q_id <- unlist(lapply(name_dis, `[[`, 1))
    res_tab$type_Eff <- unlist(lapply(name_dis, `[[`, 2))
    res_tab$par <- unlist(lapply(name_dis, `[[`, 3))
    
    # table of the QTL effect per parents
    Q_res <- vector(mode = 'list', length = nQTL)
    d_ref <- data.frame(par = par_id)
    
    for(q in 1:nQTL){
      
      res_q <- res_tab[res_tab$Q_id == paste0('QTL', q), ]
      res_q_m <- res_q[res_q$type_Eff == 'main', c('par', 'B_QTL', 'logP')]
      d_q <- merge(x = d_ref, y = res_q_m, by = 'par', all.x = TRUE)
      colnames(d_q) <- c('par', 'B_main', 'logP_main')
      
      # add the EC
      for(e in 1:n_EC){
        
        res_EC <- res_q[res_q$type_Eff == EC_id[e], c('par', 'B_QTL', 'logP')]
        colnames(res_EC) <- c('par', paste0(c('B_', 'logP_'), colnames(EC)[e])) 
        d_q <- merge(x = d_q, y = res_EC, by = 'par', all.x = TRUE)
        
      }
      
      rownames(d_q) <- d_q$par
      d_q <- d_q[par_id, ]
      d_q <- d_q[, -which(colnames(d_q) == 'par')]
      
      Q_res[[q]] <- d_q
      
    }
    
    names(Q_res) <- paste0('QTL', 1:nQTL)
    
    return(list(Qeff_main_QxE = Qeff, Qeff_EC = Q_res))
    
  } else {
    
    cat('No QTL show significant QTLxE interaction.')
    
    return(list(Qeff_main_QxE = Qeff, Qeff_EC = NULL))
    
  }
  
}