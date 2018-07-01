############
# MQE_plot #
############

#
# MQE model QTL profile
# 
# Plot the QTL profile of a multi(-QTL effect (MQE) model determined by forward
# regression (\code{\link{MQE_forward}}).
# 
# This function can be used to plot the final results of a MQE model. Once the
# list of QTL and their type of effect have been determined using
# \code{\link{MQE_forward}}, the user can introduce these values in the function
# \code{\link{MQE_CIM}} and use the argument \code{chg.Qeff = TRUE}. This means
# that \code{\link{MQE_CIM}} will compute a CIM profile with the QTL list used
# as cofactor and change the type of QTL effect of the tested position
# for the one of the detected QTL when it enter their region. Then the profile
# obtained and the QTL list can be introduced in \code{MQE_plot}. This function
# will draw the QTL profile and colour the QTL regions according to the type of
# QTL effect at the QTL position (cross-specific: black, parental: red,
# ancestral: green and bi-allelic: blue).
#
# @param mppData An object of class \code{mppData}
# 
# @param Qprof Object of class \code{QTLprof} returned by the function
# \code{\link{MQE_CIM}} with argument \code{chg.Qeff = TRUE}.
# 
# @param QTL list of QTL positions with the corresponding QTL incidence
# matrix returned by the function \code{\link{MQE_forward}}.
# 
# @param window \code{Numeric} distance (cM) on the left and the right of a
# QTL position that will be coloured according to the type of QTL effect.
# Default = 20.
# 
# @param threshold \code{Numeric} QTL significance threshold value draw on
# the plots. Default = 3.
# 
# @param main Title of the graph. Default = "MQE QTL profile".
# 
# @author Vincent Garin
# 
# @seealso \code{\link{MQE_CIM}}, \code{\link{MQE_forward}}
# 
# @examples
# 
# data(mppData)
# 
# QTL <- MQE_forward(mppData = mppData, Q.eff = c("par", "anc", "biall"))
# 
# CIM <- MQE_CIM(mppData = mppData,  Q.eff = "cr", cofactors = QTL[, 1],
#                cof.Qeff = QTL[, 5], chg.Qeff = TRUE)
# 
# MQE_plot(mppData = mppData, Qprof = CIM, QTL = QTL)
#
# @export
# 


MQE_plot <- function(mppData, Qprof, QTL, window = 30, threshold = 3, 
                                main = "MQE QTL profile") {
  
  # form the partition for the colour
  
  cofactors2 <- mppData$map[mppData$map[, 1] %in% QTL[, 1], c(2,4)]
  
  test.cof <- function(x, map, window) {
    
    t1 <- map$chr == as.numeric(x[1])
    t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
    !(t1 & t2)
    
  }
  
  cof.part <- apply(X = cofactors2, MARGIN = 1, FUN = test.cof,
                    map = mppData$map, window = window)
  
  cof.part2 <- (!cof.part)*1
  
  Qeff.partition <- function(x, cof.Qeff){
    if(sum(x) == 0){"pos" } else { cof.Qeff[max(which(x == 1))] }
  }
  
  Qeff.part <- apply(X = cof.part2, MARGIN = 1, FUN = Qeff.partition,
                     cof.Qeff = QTL[, 5])
  
  org.part <- matrix(0, dim(cof.part)[1], 4)
  
  org.part[Qeff.part == "cr", 1] <- Qprof[Qeff.part == "cr", 5]
  org.part[Qeff.part == "par", 2] <- Qprof[Qeff.part == "par", 5]
  org.part[Qeff.part == "anc", 3] <- Qprof[Qeff.part == "anc", 5]
  org.part[Qeff.part == "biall", 4] <- Qprof[Qeff.part == "biall", 5]
  
  org.part[org.part == 0] <- NA
  colnames(org.part) <- c("cr", "par", "anc", "biall")
  
  data <- data.frame(Qprof, org.part, stringsAsFactors = FALSE)
  
  # redefine data within the function to suppress R CMD check notation
  
  chr <- Qprof$chr
  pos.cM <- Qprof$pos.cM
  log10pval <- Qprof$log10pval
  cr <- org.part[, 1]
  par <- org.part[, 2]
  anc <- org.part[, 3]
  biall <- org.part[, 4]
  
  data <- data.frame(chr, pos.cM, log10pval, cr, par, anc, biall)
  
  ggplot(data, aes(x = pos.cM, y = log10pval, group = chr)) +
    geom_line(linetype = 3) + 
    geom_path(mapping = aes(y = cr), na.rm = TRUE, lwd = 0.5, colour = 1) + 
    geom_path(mapping = aes(y = par), na.rm = TRUE, lwd = 0.5, colour = 2) + 
    geom_path(mapping = aes(y = anc), na.rm = TRUE, lwd = 0.5, colour = 3) + 
    geom_path(mapping = aes(y = biall), na.rm = TRUE, lwd = 0.5, colour = 4) + 
    facet_wrap(nrow = 1, ~chr, scales = "free_x") +
    geom_hline(yintercept = threshold, colour = "red") +
    labs(title = main) + theme_bw() + xlab("position [cM]") + 
    ylab("-log10(p.val)") +
    theme(axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.text.x  = element_text(size=18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size=22),
          strip.text.x =  element_text(size=18))
  
}