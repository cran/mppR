######################
# plot_allele_eff_GE #
######################

#' plot of genome wide QTL allelic effect significance
#'
#' Plot of the genome wide significance of the QTL allelic effects in multiple
#' environments.
#'
#' @param mppData An object of class \code{mppData}.
#'
#' @param nEnv \code{Numeric} number of environment or trait.
#'
#' @param EnvNames \code{character} expression indicating the environment or trait
#' labels.
#'
#' @param Qprof object obtained with function \code{\link{mppGE_SIM}},
#' \code{\link{mppGE_CIM}},
#'
#' @param Q.eff one of "cr", "par", "anc" or "biall". For the moment only "par"
#' is available.
#'
#' @param QTL Optional argument. Object of class \code{QTLlist} representing a
#' list of selected position obtained with the function \code{\link{QTL_select}}
#' or a vector of \code{character} marker or in between marker positions names.
#' These positions will be plotted on the graph. Default = NULL.
#' 
#' @param ref_par \code{Charater} specifying the reference parent. Default = NULL.
#'
#' @param main Title of the graph. Default = "QTL genetic effects plot".
#'
#' @param text.size \code{Numeric} value specifying the size of graph axis text
#' elements. Default = 18.
#'
#'
#' @author Vincent Garin
#'
#' @seealso
#'
#' \code{\link{mppGE_CIM}}, \code{\link{mppGE_SIM}}
#'
#' @examples
#'
#' data(mppData_GE)
#'
#' SIM <- mppGE_SIM(mppData = mppData_GE, trait = c('DMY_CIAM', 'DMY_TUM'))
#'
#' Qpos <- QTL_select(Qprof = SIM, threshold = 3, window = 50)
#'
#' plot_allele_eff_GE(mppData = mppData_GE, nEnv = 2, EnvNames = c('CIAM', 'TUM'),
#'                    Qprof = SIM, Q.eff = 'par', QTL = Qpos, text.size = 14)
#'
#' @export

plot_allele_eff_GE <- function(mppData, nEnv, EnvNames, Qprof, Q.eff = 'par', QTL = NULL,
                               ref_par = NULL, main = "QTL genetic effects plot",
                               text.size = 18)
{
  
  # 1. check data format
  ######################
  
  stopifnot(inherits(Qprof, "QTLprof"))
  
  if(!(Q.eff %in% c("cr", "par", "anc", "biall"))){
    
    stop("The Q.eff argument must take value: 'cr', 'par', 'anc' or 'biall'.")
    
  }
  
  n.eff <- dim(Qprof)[2] - (5 + nEnv)
  
  if(n.eff == 0) {
    
    stop("The Qprof object does not contain any QTL p-value information.
         It was probably not obtained using plot.gen.eff = TRUE")
    
  }
  
  # Check if the reference parent is in the parent list
  if(!is.null(ref_par)){
    
    if(!(ref_par %in% mppData$parents)){
      
      av_par <- paste(mppData$parents, collapse = ', ')
      err_mes <- 'The reference parent specified (ref_par) is not in the parent list. It should be one of:'
      
      stop(paste(err_mes, av_par))
      
    }
    
  }
  
  # change environment order
  
  EnvNames <- rev(EnvNames)
  
  EnvPval <- vector(mode = "list", length = nEnv)
  
  pval <- Qprof[, (6 + nEnv):dim(Qprof)[2]]
  
  # replace NA by 1
  pval <- as.matrix(pval)
  pval[is.na(pval)] <- 1
  pval <- data.frame(pval)
  
  index <- split(1:dim(pval)[2], factor(sort(rank(1:dim(pval)[2])%%nEnv)))
  
  rev_id <- nEnv:1
  
  for(i in 1:nEnv){
    
    EnvPval[[i]] <- pval[, index[[rev_id[i]]]]
    
  }
  
  pval <- do.call(what = cbind, EnvPval)
  Qprof <- cbind(Qprof[, 1:5], pval)
  
  # 2. order columns within connected parts
  #########################################
  
  if((Q.eff == "par") || (Q.eff == "anc")){
    
    # split the pval into environment list
    
    EnvPval <- vector(mode = "list", length = nEnv)
    
    pval <- Qprof[, 6:dim(Qprof)[2]]
    
    index <- split(1:dim(pval)[2], factor(sort(rank(1:dim(pval)[2])%%nEnv)))
    
    for(i in 1:nEnv){
      
      EnvPval[[i]] <- pval[, index[[i]]]
      
    }
    
    # determine the connected parts
    
    con.part <- design_connectivity(par_per_cross = mppData$par.per.cross,
                                    plot_des = FALSE)
    
    if(length(con.part) > 1){
      
      stop("The function apply only to MPP designs composed of a single connected part")
      
    }
    
    all.ref <- apply(X = EnvPval[[1]], MARGIN = 2, FUN = function(x) sum(x == 1))
    names(all.ref) <- mppData$parents
    
    ord.par <- names(sort(all.ref))
    
    if(!is.null(ref_par)){
      ord.par <- c(ord.par[-which(ord.par == ref_par)], ref_par)
    }
    
    # order each environmental parts
    
    for(i in 1:nEnv){
      
      pval_i <- EnvPval[[i]]
      colnames(pval_i) <- mppData$parents
      pval_i <- pval_i[, ord.par]
      
      EnvPval[[i]] <- pval_i
      
    }
    
    pval <- do.call(what = cbind, EnvPval)
    
    Env_name <- rep(paste0(" ",EnvNames), each = length(ord.par))
    y.names <- paste0(rep(ord.par, nEnv), Env_name)
    
    Qprof <- cbind(Qprof[, 1:5], pval)
    
  }
  
  # 2. elements for the plot
  ##########################
  
  ### 2.1 list of QTL
  
  if(!is.null(QTL)){
    
    stopifnot(inherits(QTL, "QTLlist"))
    pos.Q <- QTL[, c(2, 4)]
    
  }
  
  
  ### 2.2 colour code from -5 red to 5 blue
  
  z <-  c(apply(X = Qprof[, 6:dim(Qprof)[2]], MARGIN = c(1, 2),
                FUN = color.code))
  
  ### 2.3 cross or parent indicator and chromosome
  
  y <- factor(rep(1:n.eff, each = dim(Qprof)[1]))
  
  chr <- rep(Qprof$chr, n.eff)
  
  ### 2.4 genetic map positions in cM with width between two positions.
  
  x <- rep(Qprof$pos.cM, n.eff)
  
  w <- tapply(X = Qprof$pos.cM, INDEX = Qprof$chr,
              FUN = function(x) c(diff(x), 1))
  w <- unlist(w)
  w <- rep(w, n.eff)
  
  pos.cM <- Qprof$pos.cM
  
  y_lab <- "parents"
  
  ### 2.5 legend for the y-axis (cross or parents)
  
  if(Q.eff == "cr") {
    
    cross.names <- unique(mppData$cross.ind)
    par.cross.names <- paste0("(", mppData$par.per.cross[, 2],
                              "x", mppData$par.per.cross[, 3], ")")
    
    cr.names <- paste(cross.names, par.cross.names, sep = "\n")
    
    Env_name <- rep(paste0(" ",EnvNames), each = length(cr.names))
    y.names <- paste0(rep(cr.names, nEnv), Env_name)
    
    y_lab <- "crosses"
    
  } else if (Q.eff == "biall"){
    
    y.names <- EnvNames
    
    y_lab <- "minor SNP"
    
  }
  
  ### 2.6 gather data for the plot
  
  data <- data.frame(x, y, z, chr, w)
  
  # coordinate for env separation
  
  y_env <- max(as.numeric(y))/nEnv
  y_env <- y_env * 1:nEnv
  y_env <- y_env[1:(nEnv-1)] + 0.5
  
  
  # 3. plot
  #########
  
  if(is.null(QTL)){ # no QTL position given
    
    pl <- ggplot(data, aes(x, y, z = z))
    pl + geom_tile(aes(fill = z, width = w)) +
      facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
      scale_fill_gradient2(limits = c(-6, 6), low = "red", mid = "white",
                           high = "blue") +
      geom_hline(yintercept= y_env) +
      theme_bw() + xlab("position [cM]") + ylab(y_lab) +
      scale_y_discrete(labels = y.names) + ggtitle(main) +
      theme(axis.title.x = element_text(size=text.size),
            axis.title.y = element_text(size=text.size),
            axis.text.x  = element_text(size=text.size),
            axis.text.y = element_text(size = text.size),
            plot.title = element_text(size=(text.size+4)),
            strip.text.x =  element_text(size=text.size),
            legend.title = element_text(size=(text.size-2)),
            legend.text = element_text(size=(text.size-2)))
    
  } else { # QTL position given
    
    pl <- ggplot(data, aes(x, y, z = z))
    pl + geom_tile(aes(fill = z, width = w)) +
      facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
      scale_fill_gradient2(limits = c(-6, 6), low = "red", mid = "white",
                           high = "blue") +
      geom_hline(yintercept= y_env) +
      geom_vline(aes(xintercept = pos.cM), pos.Q, linetype = "longdash") +
      theme_bw() + xlab("position [cM]") + ylab(y_lab) +
      scale_y_discrete(labels = y.names) + ggtitle(main) +
      theme(axis.title.x = element_text(size=text.size),
            axis.title.y = element_text(size=text.size),
            axis.text.x  = element_text(size=text.size),
            axis.text.y = element_text(size = text.size),
            plot.title = element_text(size=(text.size+4)),
            strip.text.x =  element_text(size=text.size),
            legend.title = element_text(size=(text.size-2)),
            legend.text = element_text(size=(text.size-2)))
    
  }
  
}
