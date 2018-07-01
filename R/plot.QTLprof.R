################
# plot.QTLprof #
################

#' plot QTL profile
#' 
#' Plots the -log10(p-val) profile of a QTL analysis or a genome-wide
#' genetic effect plot using package ggplot2.
#' 
#' The user can plot regular QTL profiles (\code{gen.eff = FALSE}) with
#' -log10(p-val) plotted against genetic position or genome-wide genetic
#' effects plots (\code{gen.eff = TRUE}). To plot the genome-wide genetic
#' effects, the SIM and CIM QTL profile must have been computed with
#' \code{plot.gen.eff = TRUE}.
#' 
#' The genome-wide genetic effects plots is a visualization of the significance
#' of the QTL effect per cross or per parents along the genome. For a
#' cross-specific QTL profile (\code{Q.eff = "cr"}): Red color means
#' that the allele coming from parent A(1) increases the phenotypic value and
#' parent B(2) decreases it and blue that parent A(1) decreases the trait and
#' parent B(2) increases it.
#' 
#' For a parental (\code{Q.eff = "par"}) or an ancestral model
#' (\code{Q.eff = "anc"}), the results are given per parents. The significance
#' of the effect must be interpreted as a deviation with respect to the
#' reference of each connected part. The reference allele is always defined as
#' the most frequent one. Red (Blue) colour means a significant negative
#' (positive) effect with respect to the reference of the connected part.
#' 
#' The reference parental allele can change at each position according to the
#' segregation rate. The parent are plotted from the top to the bottom according
#' to the number of time their allele is set as reference.  Therefore
#' interpretation of the genetic effect plot should be done with caution.
#' In that case, the plot should be taken as a rough indication of the signal
#' distribution.
#' 
#' The colour intensity increase with the significance of the effect.
#' There are five colour intensities according to the p-value of the QTL
#' effect: 0.05<p-val<0.01; 0.01<p-val<0.001; 0.001<p-val<0.0001;
#' 0.0001<p-val<0.00001 and p-val< 0.00001.
#' 
#' For both type of plot, the user can pass a list of cofactors or QTL position
#' to the argument \code{QTL}. These positions will be drawn on the graph using
#' dotted lines.
#' 
#' @param x Object of class \code{QTLprof} returned by the function
#' \code{\link{mpp_SIM}} or \code{\link{mpp_CIM}}.
#' 
#' @param gen.eff \code{Logical}. Specify the type of plot.
#' If \code{gen.eff = FALSE}, standard QTL profile. If \code{gen.eff = TRUE},
#' genome-wide genetic effect plot. In that case, the \code{QTLprof} object in
#' \code{x} must have been calculated with argument \code{plot.gen.eff = TRUE}.
#' Default = FALSE.
#' 
#' @param mppData An object of class \code{mppData}. Only required if
#' \code{gen.eff = TRUE}.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental
#' effects; 3) "anc" for ancestral effects. Only required if
#' \code{gen.eff = TRUE}
#' 
#' @param QTL Optional argument. List of QTL positions. Object of class
#' \code{QTLlist} representing a list of selected position obtained with the
#' function \code{\link{QTL_select}} or two columns numeric matrix with the
#' chromosome and the position in cM. These positions will be drawn on the
#' graph. Default = NULL.
#' 
#' @param type \code{Character} expression indicating the type of plot should be
#' drawn: "l" for lines , "h" for vertical bar. Default = "l".
#' 
#' @param main Title of the graph. Default = "QTL profile".
#' 
#' @param threshold \code{Numeric} QTL significance threshold value draw on
#' the plot. Default = 3.
#' 
#' @param text.size \code{Numeric} value specifying the size of graph axis text
#' elements. Default = 18.
#' 
#' @param ... Ignored.
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{mpp_SIM}}, \code{\link{mpp_CIM}}, \code{\link{QTL_select}}
#' 
#' @examples
#' 
#' data(mppData)
#' 
#' SIM <- mpp_SIM(mppData = mppData)
#' QTL <- QTL_select(SIM)
#' plot(x = SIM, QTL = QTL)
#' 
#' SIM <- mpp_SIM(mppData = mppData, Q.eff = "cr", plot.gen.eff = TRUE)
#' QTL <- QTL_select(SIM)
#' plot(x = SIM, gen.eff = TRUE, mppData = mppData, Q.eff = "cr", QTL = QTL)
#' 
#' @export
#' 

plot.QTLprof <- function(x, gen.eff = FALSE, mppData, Q.eff, QTL = NULL,
                         type = "l", main = "QTL profile", threshold = 3,
                         text.size = 18, ...)
{
  
  if(!inherits(x, "QTLprof")){stop("'x' is not of class ", dQuote("QTLprof"))}
  
  if(!gen.eff){ # Regular QTL profile.
    
    chr <- x$chr
    pos.cM <- x$pos.cM
    log10pval <- x$log10pval
    Qprof <- data.frame(chr, pos.cM, log10pval)
    
    # QTL positions
    
    if(!is.null(QTL)){
      
      if(inherits(QTL, "QTLlist")){
        
        pos.Q <- QTL[, c(2, 4)]
        
      } else{
        
        if(!((is.matrix(QTL)) & (dim(QTL)[2] == 2) & (is.numeric(QTL)))){

          stop("'QTL' must be a two columns numeric matrix with ",
                     "chromosome, and marker position")

        }
        
        pos.Q <- data.frame(QTL)
        colnames(pos.Q) <- c("chr", "pos.cM")
        
      }
      
    }
    
    if(is.null(QTL)){ # no QTL info given
      
      if (type == "l") {
        
        ggplot(Qprof, aes(x = pos.cM, y = log10pval, group = chr)) + geom_line() + 
          facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
          geom_hline(yintercept = threshold, colour = "red") + theme_bw() +
          xlab("position [cM]") + ylab("-log10(p.val)") + 
          ggtitle(main) +
          theme(axis.title.x = element_text(size=text.size),
                axis.title.y = element_text(size=text.size),
                axis.text.x  = element_text(size=text.size),
                axis.text.y = element_text(size = text.size),
                plot.title = element_text(size=(text.size + 4)),
                strip.text.x =  element_text(size=text.size))
        
        
      } else if (type == "h") {
        
        ggplot(Qprof, aes(x = pos.cM, xend = pos.cM, y = 0, yend = log10pval,
                          group = chr)) + geom_segment() +
          facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
          geom_hline(yintercept = threshold, colour = "red") + 
          theme_bw() + xlab("position [cM]") + ylab("-log10(p.val)") + 
          ggtitle(main) +
          theme(axis.title.x = element_text(size=text.size),
                axis.title.y = element_text(size=text.size),
                axis.text.x  = element_text(size=text.size),
                axis.text.y = element_text(size = text.size),
                plot.title = element_text(size=(text.size + 4)),
                strip.text.x =  element_text(size=text.size))
        
      }
      
    } else { # QTL info given
      
      if (type == "l") {
        
        ggplot(Qprof, aes(x = pos.cM, y = log10pval, group = chr)) + geom_line() + 
          facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
          geom_vline(aes(xintercept = pos.cM), pos.Q, linetype = "longdash",
                     colour = "black") +
          geom_hline(yintercept = threshold, colour = "red") + theme_bw() +
          xlab("position [cM]") + ylab("-log10(p.val)") + 
          ggtitle(main) + 
          theme(axis.title.x = element_text(size=text.size),
                axis.title.y = element_text(size=text.size),
                axis.text.x  = element_text(size=text.size),
                axis.text.y = element_text(size = text.size),
                plot.title = element_text(size=(text.size + 4)),
                strip.text.x =  element_text(size=text.size))
        
        
      } else if (type == "h") {
        
        ggplot(Qprof, aes(x = pos.cM, xend = pos.cM, y = 0, yend = log10pval,
                          group = chr)) + geom_segment() +
          facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
          geom_vline(aes(xintercept = pos.cM), pos.Q, linetype = "longdash",
                     colour = "black") +
          geom_hline(yintercept = threshold, colour = "red") + 
          theme_bw() + xlab("position [cM]") + ylab("-log10(p.val)") + 
          ggtitle(main) +
          theme(axis.title.x = element_text(size=text.size),
                axis.title.y = element_text(size=text.size),
                axis.text.x  = element_text(size=text.size),
                axis.text.y = element_text(size = text.size),
                plot.title = element_text(size=(text.size + 4)),
                strip.text.x =  element_text(size=text.size))
        
      }
      
    }
    
  } else { # Genome-wide genetic effect plot.
    
    if(!(Q.eff %in% c("cr", "par", "anc"))){
      
      stop("'Q.eff' must be ", dQuote("cr"), ', ', dQuote("par"), ', ',
           dQuote("anc"), ' or ', dQuote("biall"))
      
    }
    
    n.eff <- dim(x)[2] - 5
    
    if(n.eff == 0) {
      
      stop("'x' does not contain any QTL p-value information. ",
           "Use plot.gen.eff = TRUE, when you compute the QTL profile")
      
    }
    
    # 2. order columns within connected parts
    #########################################
    
    if((Q.eff == "par") || (Q.eff == "anc")){
      
      # determine the connected parts
      
      con.part <- design_connectivity(par_per_cross = mppData$par.per.cross,
                                      plot_des = FALSE)
      
      len.con <- unlist(lapply(X = con.part, FUN = function(x) length(x)))
      con.part <- con.part[names(sort(len.con, decreasing = FALSE))]
      
      
      allele_order <- c()
      pval <- data.frame(row.names = 1:dim(x)[1])
      ref.ind <- length(con.part) + 1 
      
      for(i in seq_along(con.part)){
        
        con.part_i <- con.part[[i]]
        
        # subset results of the connected part
        
        pval_i <- x[, con.part_i]
        
        # order according to number of reference values
        
        all.ref <- apply(X = pval_i, MARGIN = 2,
                         FUN = function(x) sum(x == 1))
        
        pval <- cbind.data.frame(pval, pval_i[, names(sort(all.ref))])
        
        allele_ord_i <- names(sort(all.ref))
        
        allele_ord_i <- c(paste(allele_ord_i, paste0("(c", (ref.ind-i),")"),
                                sep = "\n"))
        
        allele_order <- c(allele_order, allele_ord_i)
        
      }
      
      # Rename x
      
      x <- cbind(x[, 1:5], pval)
      
      y.names <- allele_order
      
    }
    
    # 2. elements for the plot
    ##########################
    
    ### 2.1 list of QTL
    
    if(!is.null(QTL)){
      
      stopifnot(inherits(QTL, "QTLlist"))
      
      pos.Q <- QTL[, c(2, 4)]
      
    }
    
    
    ### 2.2 colour code from -5 red to 5 blue
    
    z <-  c(apply(X = x[, 6:dim(x)[2]], MARGIN = c(1, 2),
                  FUN = color.code))
    
    ### 2.3 cross or parent indicator and chromosome
    
    y <- factor(rep(1:n.eff, each = dim(x)[1]))
    
    chr <- rep(x$chr, n.eff)
    
    ### 2.4 genetic map positions in cM with width between two positions.
    
    x.pos <- rep(x$pos.cM, n.eff)
    
    w <- tapply(X = x$pos.cM, INDEX = x$chr,
                FUN = function(x) c(diff(x), 1))
    w <- unlist(w)
    w <- rep(w, n.eff)
    
    pos.cM <- x$pos.cM
    
    y_lab <- "parents"
    
    ### 2.5 legend for the y-axis (cross or parents)
    
    if(Q.eff == "cr") {
      
      cross.names <- unique(mppData$cross.ind)
      par.cross.names <- paste0("(", mppData$par.per.cross[, 2], 
                                "x", mppData$par.per.cross[, 3], ")")
      
      y.names <- paste(cross.names, par.cross.names, sep = "\n")
      
      y_lab <- "crosses"
      
    } 
    
    ### 2.6 gather data for the plot
    
    data <- data.frame(x.pos, y, z, chr, w)
    
    
    # 3. plot
    #########
    
    if(is.null(QTL)){ # no QTL position given
      
      pl <- ggplot(data, aes(x.pos, y, z = z))
      pl + geom_tile(aes(fill = z, width = w)) +
        facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
        scale_fill_gradient2(limits = c(-5, 5), low = "red", mid = "white",
                             high = "blue") +
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
      
      pl <- ggplot(data, aes(x.pos, y, z = z))
      pl + geom_tile(aes(fill = z, width = w)) +
        facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
        scale_fill_gradient2(limits = c(-5, 5), low = "red", mid = "white",
                             high = "blue") +
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
  
}