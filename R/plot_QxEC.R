#############
# plot_QxEC #
############

#' plot QTLxEC effect
#' 
#' Plot allowing the visualization of the QTL parental allelic effect variation
#' given an environmental covariate (EC). The function plot the sensitivity
#' curve of the parent allelic effects.
#' 
#' @param Qeff output from the function \code{\link{QTL_effect_main_QxEC}}.
#' 
#' @param EC \code{Numeric} matrix containing the EC values of a single covariate
#' with environments as row and EC as column.
#' 
#' @param env_id \code{Character} vector specifying the environment names.
#' By default, E1, ... En
#' 
#' @param QTL \code{Numeric value} indicating which QTL to plot
#' 
#' @param sign_thre \code{Numeric value} indicating the significance threshold for a
#' parent sensitivity slope to be ploted. Default = 0.05
#' 
#' @param EC_id \code{Character} string indicating the name of the environmental covariate.
#' Default = 'EC'.
#' 
#' @param trait_id \code{Character} string indicating the name of the trait.
#' Default = 'trait'.
#'  
#' @param main \code{Character} string title of the plot. Default = 'QTLxEC'
#' 
#' @param col_vec \code{Character} vector specifying colors for the parent sensitivity lines.
#' Default = NULL
#' 
#' @param text_size \code{Numerical} value specifying the size of the text in the plot.
#' Default = 14.
#' 
#' @return
#' 
#' QTLxEC sensitivity plot
#' 
#' @author Vincent Garin
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
#' pl <- plot_QxEC(Qeff, EC = EC, env_id = c('CIAM', 'TUM', 'INRA', 'KWS'), 
#'                 QTL = 2, EC_id = 'cum rain', trait_id = 'DMY')
#' 
#' }
#' 
#' @export
#' 

plot_QxEC <- function(Qeff, EC, env_id = NULL, QTL, sign_thre = 0.05, EC_id = 'EC',
                      trait_id = "trait", main = "QTLxEC", col_vec = NULL,
                      text_size = 14){
  
  d_B_slope <- Qeff$Qeff_EC[[QTL]]
  if(all(is.na(d_B_slope[, 3]))){stop("No sensitivity slope were estimated for the selected QTL. Probably because there were no significant parental QEI.")}
  d_EC <- data.frame(EC)
  if(is.null(env_id)){env_id <- paste0("E", 1:nrow(EC))}
  d_EC$env <- env_id
  par_nm <- rownames(d_B_slope)
  
  # calculate the four points for each parents
  EC <- sort(unique(d_EC[, 1]))
  nEnv <- length(EC) 
  EC %*% t(d_B_slope[, 3])
  Qp_EC <- EC %*% t(d_B_slope[, 3]) + t(d_B_slope[, 1]) %x% matrix(rep(1, nEnv))
  
  d <- data.frame(EC = EC, Qeff = c(Qp_EC), parents = rep(par_nm, each = nEnv))
  d <- d[complete.cases(d), ]
  
  # modification of parent names (organised per significance)
  d_par <- data.frame(par = rownames(d_B_slope), sign = d_B_slope[, 4])
  d_par <- d_par[order(d_par$sign, decreasing = TRUE), ]
  d_par <- d_par[complete.cases(d_par), ]
  e_val <- sapply(X = 10^(-d_par$sign), sign.star)
  e_val[e_val == ""] <- "ns"
  
  par_sign_lk <- paste(d_par$par, e_val)
  names(par_sign_lk) <- d_par$par
  d[, 3] <- par_sign_lk[d[, 3]]
  d[, 3] <- factor(x = d[, 3], levels = par_sign_lk)
  
  n_par <- nrow(d_par)
  n_sign <- sum(d_par$sign > -log10(sign_thre))
  n_n_sign <- (n_par - n_sign)
  
  # plot
  p <- ggplot(d, aes(x = EC, y = Qeff, group = .data$parents, col = .data$parents)) +
    geom_line(aes(linetype = .data$parents)) +
    
    # format of the slope (colour, shape)
    {if(!is.null(col_vec)) scale_color_manual(values = col_vec) } +
    scale_linetype_manual(values = c(rep(1, n_sign), rep(3, n_n_sign))) +
    scale_size_manual(values = c(rep(1, n_sign), rep(1, n_n_sign))) +
    
    labs(x = EC_id, y = trait_id, title = main) +
    
    theme(axis.title.x = element_text(size = text_size-1),
          axis.title.y = element_text(size = text_size-2),
          axis.text.x = element_text(size = text_size-2),
          axis.text.y = element_text(size = text_size-2),
          plot.title = element_text(size = (text_size + 1)),
          strip.text.x = element_text(size = text_size),
          legend.title = element_text(size = (text_size-2)),
          legend.text = element_text(size = (text_size-3)))
  
  return(p)
  
}