#############
# EC_effect #
#############

#' Determine EC effects
#' 
#' Determine the effect of environmental covariates (EC) on the mean of a
#' trait across environments and the time range where this effect is the
#' strongest. The procedure was originally proposed by Li et al. (2018).
#' 
#' @param trait_env_mean vector of trait mean over environment.
#' 
#' @param crop_duration numerical value indicating the crop duration.
#' 
#' @param EC_list list EC parameter matrix. one per environment. The
#' order of the environment must be the same as the one of the trait
#' mean.
#' 
#' @param type character string vector indicating the type of statistic
#' that correspond to the EC. Either the cumulated sum ('sum') or the
#' average ('mean').
#' 
#' @param min_win Numerical value indicating the minimum size of the range
#' between start and end day when the EC values are measured. Default = 20
#' 
#' @param sel_criteria Character specifying the selection criteria.
#' Default = 'global'
#' 
#' @param plot Logical value indicating if a plot of the EC effects over time
#' should be returned. Default = FALSE,
#' 
#' @param plot_dir Directory where the plot should be returned. Default = NULL
#' 
#' @param p_title Title of the plot. Default = 'EC_plot'
#' 
#' @param env_nm Optional vector of environment name. Default = NULL.
#' 
#' @return Return:
#'
#' \code{data.frame} that contains the following elements for each EC (line).
#' The first line is the value of the EC in the different environments:
#' 
#' \enumerate{
#' 
#' \item{Start and end date of the optimal window}
#' \item{R2 of correlation between trait and EC}
#' \item{Direction of the correlation}
#' \item{Average R2 value over all tested windows}
#' \item{EC value in the different environments for the optimal time window}
#' 
#' }
#'
#' @author Vincent Garin
#' 
#' @references
#' 
#' Li, X., Guo, T., Mu, Q., Li, X., & Yu, J. (2018). Genomic and environmental
#' determinants and their interplay underlying phenotypic plasticity.
#' Proceedings of the National Academy of Sciences, 115(26), 6679-6684.
#' 
#' @export


EC_effect <- function(trait_env_mean, crop_duration, EC_list, type,
                      min_win = 20, sel_criteria = 'global', plot = TRUE,
                      plot_dir = NULL, p_title = 'EC_plot', env_nm = NULL){
  
  # sub function to calculate the EC score
  EC_val_fct <- function(x, col, start, end, type){
    
    EC <- x[start:end, col]
    
    if(type == 'sum'){
      m_EC <- mean(EC, na.rm = TRUE)
      EC[is.na(EC)] <- m_EC
      sum(EC)
    } else {
      mean(EC, na.rm = TRUE)
    }
    
  }
  
  EC <- colnames(EC_list[[1]])
  n_EC <- length(EC)
  n_env <- length(trait_env_mean)
  
  # determine all the possible EC windows
  win_mat <- EC_window_mat(crop_duration = crop_duration)
  
  # table to store the end, start, R2 and EC value
  EC_res <- matrix(NA, nrow = n_EC+1, ncol = 5 + n_env)
  EC_res[1, ] <- c(1, crop_duration, NA, NA, NA, trait_env_mean)
  
  if(plot){
    
    pdf(file = file.path(plot_dir, paste0(p_title, '.pdf')))
    
  }
  
  # iterate over all EC
  for(e in 1:n_EC){
    EC_res_r2 <- matrix(NA, nrow = nrow(win_mat), ncol = ncol(win_mat))
    EC_res_sign <- matrix(NA, nrow = nrow(win_mat), ncol = ncol(win_mat))
    # iterate over all the possible window
    for(i in 1:nrow(win_mat)){
      for(j in 1:ncol(win_mat)){
        
        if(!is.na(win_mat[i, j])){
          
          EC_val <- lapply(X = EC_list, FUN = EC_val_fct, col = e, start = j,
                           end = win_mat[i, j], type = type[e])
          EC_val <- unlist(EC_val)
          cor_EC_e <- cor(trait_env_mean, EC_val)
          EC_res_r2[i, j] <- (cor_EC_e^2)
          EC_res_sign[i, j] <- sign(cor_EC_e)
          
        }
      }
    }
    
    # data for optim result find and plot
    x <- rep(1:ncol(EC_res_r2), each = nrow(EC_res_r2))
    y <- rep(as.numeric(rownames(win_mat)), time = ncol(EC_res_r2))
    end_d <- c(win_mat)
    
    d_res <- data.frame(start_day = x, win_size = y, end_day = end_d,
                        R2 = c(EC_res_r2), sign = c(EC_res_sign))
    
    d_sort <- d_res[order(d_res$R2, decreasing = TRUE), ]
    start_opt <- d_sort[1, 1]
    end_opt <- d_sort[1, 3]
    R2 <- d_sort$R2[1]
    R2_sign <- d_sort$sign[1]
    R2_glb <- mean(d_sort$R2, na.rm = TRUE)
    
    # calculate the EC_val at the optimum time
    
    EC_val_opt <- lapply(X = EC_list, FUN = EC_val_fct, col = e,
                         start =  start_opt, end = end_opt, type = type[e])
    
    # fill the table
    EC_res[e+1, ] <-  c(start_opt, end_opt, R2, R2_sign, R2_glb, unlist(EC_val_opt))
    
    if(plot){
      
      pl <- ggplot(d_res, aes_string(x=.data$start_day, y=.data$end_day)) +
        geom_tile(aes(fill = R2 * sign)) +
        scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                             guide = guide_colorbar(order = 1)) +
        ggtitle(EC[e])
      print(pl)
      
    }
    
  }
  
  if(plot){dev.off()}
  
  rownames(EC_res) <- c('trait_val', EC)
  if(!is.null(env_nm)){
    colnames(EC_res) <- c('start', 'end', 'R2', 'R2_sign', 'R2_glb', env_nm)
  } else {
    colnames(EC_res) <- c('start', 'end', 'R2', 'R2_sign', 'R2_glb', paste0('EC_E', 1:n_env))
  }
  
  EC_res <- data.frame(EC_res)
  
  return(EC_res)
  
}