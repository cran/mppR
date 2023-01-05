###############
# plot_QTLxEC #
###############

#' plot QTLxEC effect
#' 
#' Plot allowing the visualisation of the QTL allelic effect given an
#' environmental covariate (EC). It represents parental QTL effects that
#' significantly interact with the EC. Those values are added to the (central)
#' reference parent effect (intercept) which allow a comparison of the
#' parental allele contribution with respect to the reference and in the scale
#' of the trait.
#' 
#' @param Qeff output from the function \code{\link{QTL_effect_QxEC}} obtained
#' with option QTLxEC_plot = TRUE.
#' 
#' @param Q_id \code{Numeric value} indicating which QTL to plot
#' 
#' @param RP \code{Character} string indicating the name of the reference (central) parent.
#' Default = 'RP'.
#' 
#' @param EC_id \code{Character} string indicating the name of the environmental covariate.
#' Default = 'EC'.
#' 
#' @param trait_id \code{Character} string indicating the name of the trait.
#' Default = 'trait'.
#'  
#' @param main \code{Character} string title of the plot. Default = 'QTLxEC'
#' 
#' @param keep_par \code{Character} vector or string specifying the only parents that
#' should be ploted. Default = NULL.
#'  
#' @param rem_par \code{Character} vector or string specifying the parents that
#' should not be ploted. Default = NULL.
#' 
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
#' Qeff <- QTL_effect_QxEC(mppData = mppData_GE,
#'                          trait = c('DMY_CIAM', 'DMY_TUM', 'DMY_INRA_P', 'DMY_KWS'),
#'                          env_id = c('CIAM', 'TUM', 'INRA', 'KWS'),
#'                          QTL = Qpos, EC = EC)
#' 
#' pl <- plot_QTLxEC(Qeff, Q_id = 2, RP = "UH007", EC_id = 'cum rain',
#' trait_id = 'DMY')
#' 
#' }
#' 
#' @export
#' 

plot_QTLxEC <- function(Qeff, Q_id, RP = "RP", EC_id = 'EC', trait_id = 'trait',
                        main = 'QTLxEC', keep_par = NULL, rem_par = NULL,
                        text_size = 14){
  
  # define variables initially (silence note CRAN)
  env <- EC <- cr_env_int <- int <- slope <- min_EC <- max_EC <- p_max <- NULL
  parents <- tr_val <- desc <- NULL
  
  d <- Qeff$Q_res_plot
  
  # select the QTL column
  d <- d[, c(1:5, which(colnames(d) == paste0('QTL', Q_id)))]
  colnames(d)[ncol(d)] <- 'QTL'
  
  if (all(is.na(d$QTL))) {
    
    stop('no significant term for this QTL.')
    
  }
  
  d <- d[!is.na(d$QTL), ]
  d$tr_val <- d$cr_env_int + d$QTL
  
  RP_tr <- d %>% group_by(env) %>% summarise(EC = mean(EC), tr_val = mean(cr_env_int))
  RP_tr$par <- RP
  
  d_p <- data.frame(EC = c(RP_tr$EC, d$EC), par = c(RP_tr$par, d$par),
                    tr_val = c(RP_tr$tr_val, d$tr_val))
  d_p$par <- factor(d_p$par, levels = unique(c(RP_tr$par, d$par)))
  
  # estimate intercept and slope
  int_f <- function(x, y){ m <- lm(y ~ x); coef(m)[1]}
  slope_f <- function(x, y){ m <- lm(y ~ x); coef(m)[2]}
  
  d_slope <- d_p %>% group_by(par) %>% summarise(int = int_f(x = EC, y = tr_val),
                                                 slope = slope_f(x = EC, y = tr_val),
                                                 min_EC = min(EC),
                                                 max_EC = max(EC)) %>%
    mutate(p_min = int + (slope * min_EC), p_max = int + (slope * max_EC)) %>%
    arrange(desc(p_max))
  
  # sort per parents
  d_slope$par <- factor(d_slope$par, levels = unique(d_slope$par))
  d_slope <- data.frame(parents = rep(d_slope$par, 2),
                        EC = c(d_slope$min_EC, d_slope$max_EC),
                        tr_val = c(d_slope$p_min, d_slope$p_max))
  
  # keep or remove parents
  
  if(!is.null(keep_par)){
    d_slope <- d_slope[(d_slope$parents %in% c(RP, keep_par)), ]
  }
  
  if(!is.null(rem_par)){
    d_slope <- d_slope[!(d_slope$parents %in% rem_par), ]
  }
  
  linetype_val <- rep(1, length(unique(d_slope$parents)))
  linetype_val[unique(d_slope$parents) == RP] <- 3
  
  pl <- ggplot(d_slope, aes(x = EC, y = tr_val, group = parents, col = parents)) +
    geom_line(size = 1.2, aes(linetype = parents)) +
    scale_linetype_manual(values = linetype_val) +
    labs(x = EC_id, y = trait_id, title = main) +
    theme(axis.title.x = element_text(size=text_size),
          axis.title.y = element_text(size=text_size),
          axis.text.x  = element_text(size = text_size),
          axis.text.y = element_text(size = text_size),
          plot.title = element_text(size=(text_size+4)),
          strip.text.x =  element_text(size=text_size),
          legend.title = element_text(size=(text_size)),
          legend.text = element_text(size=(text_size)))
  
  return(pl)
  
}