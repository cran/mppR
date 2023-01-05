###################
# QTL_effect_QxEC #
###################

#' Estimation of QTL effect sensitivity to environmental covariates
#'
#'
#' Determination of which parental QTL effect show a significant interaction
#' with the environment. Then, the function try to characterize the nature
#' of the QTLxE effect by estimating the sensitivity of the parental allelic
#' effects showing significant QTLxE interaction to environmental covariates
#' provided by the user.
#' 
#' The function first estimate the parental QTL allele main and QTLxE effect
#' using the function \code{\link{QTL_effect_main_QxE}}. Then it determines
#' which parental allele shows a significant QTLxE effect by looking if the
#' -log10(p-val) of the parental QTLxE effect is superior or equal to
#' \code{thre_QTL} and if the -log10(p-val) of QTLxE term is superior to one of
#' the main effect. Finally, given this information, the function replaces the
#' QTLxE term of the parental QTL allelic effect showing a significant QTLxE
#' effect with a main effect and QTLxEC term representing interaction between
#' the parental QTL allele and the environmental covariate (EC). The QTLxEC term can
#' be interpreted as a sensitivity of the QTL to the variation of the EC in the
#' different environments.
#' 
#' Two options are possible concerning the inclusion of the parental QTL allele
#' as main effect in the QTLxEC model. Either all parental allele are introduced
#' (\code{all_main = TRUE}, default), or only the parental allele showing a
#' singificant main effect are introduced (\code{all_main = FALSE}). 
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
#' @param thre_QTL \code{Numerical} value specifying the -log10(p-val) threshold
#' for a parental QTL allele to be considered as significant. By default,
#' thre_QTL = 2, which correspond to a p-value of 0.01.
#' 
#' @param all_main \code{Logical} value specifying if all the parental alleles
#' should be set as main effect in the QTLxEC model or if only the significant
#' parental allele should be introduced in the model as main effect and QTLxEC
#' effect if the QTLxE term is significant. Default = TRUE.
#' 
#' @param EC \code{Numeric} matrix with environments as row and environmental
#' covariates (EC) as column. The cell i, j of EC specify the value of the
#' jth EC in environment i.
#' 
#' @param Qmain_QxE results from \code{\link{QTL_effect_main_QxE}}
#' 
#' @param QTLxEC_plot \code{Logical} value specifying if the data to
#' plot sensitivity curve with the function plot_QTLxEC should be returned.
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
#' @seealso
#' 
#' \code{\link{QTL_effect_main_QxE}}
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
#' Qeff <- QTL_effect_QxEC(mppData = mppData_GE,
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

QTL_effect_QxEC <- function(mppData, trait, env_id = NULL, VCOV = "UN",
                            ref_par = NULL, QTL = NULL, QmainQi = TRUE,
                            thre_QTL = 2, all_main = TRUE, EC,
                            Qmain_QxE = NULL, QTLxEC_plot = TRUE,
                            maxIter = 100, msMaxIter = 100){
  
  #### 1. Check some arguments ####
  if(!is.matrix(EC)){stop('EC is not a matrix')}
  if(!is.numeric(EC)){stop('EC is not (fully) numeric.')}
  
  if(nrow(EC) != length(trait)){
    stop('The number of EC row is not equal to the number of environments specified in trait.')
  }
  
  #### 2. Estimation of the QTL main and QxE effects ####
  
  if(!is.null(Qmain_QxE)){
    
    Qeff <- Qmain_QxE
    
  } else {
    
    Qeff <- QTL_effect_main_QxE(mppData = mppData, trait = trait,
                                env_id = env_id, VCOV = VCOV, ref_par = ref_par,
                                QTL = QTL, QmainQi = QmainQi, maxIter = maxIter,
                                msMaxIter = msMaxIter)
  }
  
  
  
  #### 3. Determine the significance of QTLxE ####
  
  Q_dist <- determine_Qmain_QxE_sign(Qeff = Qeff, par_id = mppData$parents,
                                     thre_QTL = thre_QTL)
  
  # Check if there is at least some QTLxE term significant
  QxE_sign_vect <- c(Q_dist$QTLxE)
  if(any(QxE_sign_vect[!is.na(QxE_sign_vect)])){
    
    if(all_main){
      
      Q_dist_main <- Q_dist[[1]]
      Q_dist_main[Q_dist_main == FALSE] <- TRUE
      Q_dist[[1]] <- Q_dist_main
      
    }
    
    #### 4. Form required elements for the analysis ####
    nEnv <- length(trait)
    TraitEnv <- c(mppData$pheno[, trait])
    NA_id <- is.na(TraitEnv)
    
    env_id <- paste0('E', 1:nEnv)
    EC_id <- paste0('EC', 1:ncol(EC))
    EC_org_nm <- colnames(EC)
    n_EC <- length(EC_id)
    par_id <- mdf_par_name(mppData$parents)
    
    #### 5. QTL matrices ####
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
    
    #### 6. modify the (parent) names of QTL list and Q_dist info ####
    
    for(i in 1:length(QTL_list)){
      
      colnames(QTL_list[[i]]) <- mdf_par_name(colnames(QTL_list[[i]]))
    }
    
    colnames(Q_dist$QTL_main) <- mdf_par_name(colnames(Q_dist$QTL_main))
    colnames(Q_dist$QTLxE) <- mdf_par_name(colnames(Q_dist$QTLxE))
    
    #### 7. General element for MM ####
    
    nGeno <- dim(mppData$pheno)[1]
    env <- rep(paste0('E', 1:nEnv), each = nGeno)
    cross <- rep(mppData$cross.ind, nEnv)
    geno <- rep(rownames(mppData$pheno), nEnv)
    cross_env <- paste0(cross, '_', env)
    d <- data.frame(trait = TraitEnv, env = env, cross_env = cross_env, geno = geno)
    d[, 2:4] <- lapply(d[, 2:4], as.factor)  
    
    
    #### 8. QTL main effect matrix ####
    
    QTL_main <- vector(mode = 'list', length = nQTL)
    QTL_main_nm <- c()
    M_env <- matrix(rep(1, nEnv), ncol = 1)
    
    for(i in 1:nQTL){
      
      Qmat_i <- QTL_list[[i]]
      P_sign <- Q_dist[[1]][i, ]
      P_sign <- P_sign[colnames(Qmat_i)]
      
      QTL_main[[i]] <- M_env %x% Qmat_i[, P_sign, drop = FALSE]
      
      Q_mn_nm_i <- paste0(paste0('QTL', i), '_main_', names(P_sign)[P_sign])
      QTL_main_nm <- c(QTL_main_nm, Q_mn_nm_i)
      
    }
    
    QTL_mat_main <- do.call(cbind, QTL_main) 
    colnames(QTL_mat_main) <- QTL_main_nm
    
    #### 9. QTL main + EC mixed model computation ####
    
    QTL_mat_EC <- form_QTLxEC_mat(QTL_list = QTL_list, EC = EC, Q_dist = Q_dist)
    M_i <- MM_QTL(d = data.frame(d, QTL_mat_main, QTL_mat_EC),
                  VCOV = VCOV, maxIter = 100, msMaxIter = 100)
    m <- M_i$m
    rm(M_i)
    
    #### 10. results processing ####
    
    Beta <- m$coefficients$fixed
    
    if(QTLxEC_plot){
      
      # reference table
      tab_ref <- expand.grid(paste0('E', 1:nEnv) , mppData$par.per.cross[, 1],
                             stringsAsFactors = FALSE)
      tab_ref <- tab_ref[, 2:1]
      colnames(tab_ref) <- c('cross', 'env')
      tab_ref$cr_env_int <- NA
      
      cr_env_int <- Beta[grep(pattern = 'cross_env', x = names(Beta))]
      names(cr_env_int) <- substr(x = names(cr_env_int), 10, nchar(names(cr_env_int)))
      tab_ref$cr_env_int <- cr_env_int[paste0(tab_ref$cross, '_', tab_ref$env)]
      
      par_lk <- mdf_par_name(mppData$par.per.cross[, 3])
      names(par_lk) <- mppData$par.per.cross[, 1]
      
      d_int <- data.frame(cross = tab_ref$cross, par = par_lk[tab_ref$cross],
                          env = tab_ref$env, EC = NA, cr_env_int = tab_ref$cr_env_int)
      
      # add EC information
      EC_lk <- EC[, 1]
      names(EC_lk) <- paste0('E', 1:nEnv)
      d_int$EC <- EC_lk[d_int$env]
      
    }
    
    # QTL effect
    Q_ind <- grepl(pattern = 'QTL', x = names(Beta))
    B_QTL <- Beta[Q_ind]
    B_QTL_var <- diag(m$varFix)[Q_ind]
    
    W_Qa <- rep(NA, length(B_QTL))
    for(q in 1:length(W_Qa)) W_Qa[q] <- (B_QTL[q]^2)/B_QTL_var[q]
    W_Qa <- pchisq(W_Qa, 1, lower.tail = FALSE)
    sign_star <- sapply(X = W_Qa, sign.star)
    
    res_tab <- data.frame(B_QTL, p_val = W_Qa, logP = -log10(W_Qa),
                          sign = sign_star)
    
    name_dis <- strsplit(x = names(B_QTL), split = '_')
    res_tab$Q_id <- unlist(lapply(name_dis, `[[`, 1))
    res_tab$type_Eff <- unlist(lapply(name_dis, `[[`, 2))
    res_tab$par <- unlist(lapply(name_dis, `[[`, 3))
    
    # table of the QTL effect per parents
    Q_res <- vector(mode = 'list', length = nQTL)
    d_ref <- data.frame(par = par_id)
    
    if(QTLxEC_plot){
      Q_res_plot <- matrix(NA, nrow = nrow(d_int), ncol = nQTL)
      colnames(Q_res_plot) <- paste0('QTL', 1:nQTL)
    }
    
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
      
      if(QTLxEC_plot){
        
        # QTL sensitivity coeff
        d_QS <- d_q[, c(1, 3, 4)]
        d_QS <- d_QS[!is.na(d_QS[, 2]) & d_QS[, 3] > thre_QTL, ]
        
        if(nrow(d_QS) > 0){
          
          Qp_EC <- EC %*% t(d_QS[, 2]) +  t(matrix(d_QS[, 1])) %x% matrix(rep(1, nEnv))
          Q_res_plot[d_int$par %in% rownames(d_QS), q] <- c(Qp_EC)
          
        }
        
      }
      
    }
    
    names(Q_res) <- paste0('QTL', 1:nQTL)
    
    if(QTLxEC_plot){
      
      Q_res_plot <- data.frame(d_int, Q_res_plot)
      
      return(list(Qeff_main_QxE = Qeff, Qeff_EC = Q_res, Q_res_plot = Q_res_plot))
      
    } else {
      
      return(list(Qeff_main_QxE = Qeff, Qeff_EC = Q_res))
      
    }
    
    
    
    
  } else {
    
    return(list(Qeff_main_QxE = Qeff, Qeff_EC = NULL))
    
    cat('No QTL show significant QTLxE interaction.')
    return(NULL)
    
  }
  
}