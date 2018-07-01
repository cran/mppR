###########
# plot_CV #
###########

# Plot cross-validation QTL profiles
# 
# Plot cross-validation (CV) QTL profiles using transparency. The
# black straight lines are proportional to the number of time a position has
# been detected during the entire CV process. Transparency plot are based on
# an idea by Pieter Jongsma (see reference).
# 
# @param CV.res Object returned by function \code{\link{mpp_CV}}.
# 
# @param main Title of the graph. Default = "CV QTL profiles".
# 
# @author Vincent Garin
# 
# 
# @references
#
# https://github.com/pieterjongsma/RProfilePlot
# 
# @seealso \code{\link{mpp_CV}}
# 
# @examples
# 
# \dontrun{
# 
# data(mppData)
# 
# my.loc <- "C:/..."
# 
# CV <- mpp_CV(mppData = mppData, her = .5, output.loc = my.loc, Rep = 3)
# 
# plot_CV(CV.res = CV)
# 
# }
# 
# 
# @export
#     


plot_CV <- function(CV.res, main = "CV QTL profiles") {
  
  # arrange data
  
  map <- CV.res$QTL.profiles[, 1:4]
  profiles <- CV.res$QTL.profiles[, 6:dim(CV.res$QTL.profiles)[2]]
  n.prof <- dim(profiles)[2]
  
  # form a map with cumulated length
  
  # get the chromosome lengths
  
  chr.l <- tapply(X = map[, 4], INDEX = factor(map[, 2]),
                  FUN = function(x) max(x))
  
  l.chr <- cumsum(chr.l[-length(chr.l)])
  
  cum.len <- c(0, l.chr)
  
  # addition the cumulated length to the individual positions indicators
  
  
  ind.chr.pos <- split(x = map[, 4], f = factor(map[, 2]))
  
  vect.cum.l <- mapply(FUN = function(pos, cum.len) pos + cum.len,
                       pos = ind.chr.pos, cum.len = cum.len)
  
  vect.cum.l <- unlist(vect.cum.l)
  
  
  # get the chromosome border
  
  chr.border <- c(cum.len, max(vect.cum.l))
  
  # get y max
  
  max.y <- max(unlist(profiles))
  
  # make the vector of number of time that a QTL appeared in the CV process
  # proportional to the max lod values.
  
  N_QTL.vect <- CV.res$QTL.profiles[, 5] * (max.y/n.prof)
  
  # points data to superimpose the number of QTLs on the profiles plots
  
  points.data <- cbind(vect.cum.l, N_QTL.vect)
  
  # plot the frame according to the most extreme positions values and LOD
  
  chr.border.dif <- diff(chr.border)
  
  pos.label <- c()
  pos.cum <- c()
  
  for (i in 1:length(chr.border.dif)) {
    
    seq.i <- seq(from = 0, to = chr.border.dif[i], by = 50)
    pos.label.i <- seq.i[ -1]
    pos.label <- c(pos.label, pos.label.i)
    
    pos.cum <- c(pos.cum, pos.label.i + chr.border[i])
    
  }
  
  # determine an extra-space above to add the chromosome number on the plot
  
  max.y2 <- 1.05*max.y
  y.pos.chr.lab <- max.y + ((max.y2 - max.y)/1.2)
  x.pos.chr.lab <- diff(chr.border)/2 + chr.border[-length(chr.border)]
  
  plot(0, xaxt = "n", xlim = c(min(vect.cum.l), max(vect.cum.l)),
       ylim = c(0, max.y2), type = "n", xlab = "position [cM]",
       ylab = "-log10(p-value)", main = main, cex.main = 1.8,
       cex.axis = 1.5, cex.lab = 1.5)
  
  # abline(h = max.y)
  segments(x0=0, y0=max.y, x1=max(vect.cum.l), y1=max.y)
  
  
  text(x = x.pos.chr.lab, y = y.pos.chr.lab, labels = 1:length(chr.l),
       cex = 1.5)
  
  axis(side = 1, at = pos.cum, labels = as.character(pos.label), cex.axis = 1.5)
  
  fill <- rgb(1, 0, 0, 0.15)
  
  # x vector stays the same
  
  x <- vect.cum.l
  
  for (i in 1:n.prof) {
    
    # changes sucessively the vector of -log10(pval)
    
    ys <- profiles[, i]
    
    x.adj <- c(min(x), x, max(x))
    y.adj <- c(0, ys, 0)
    
    polygon(x.adj, y.adj, col = fill, border = rgb(0, 0, 0, 0))
    
  }
  
  # plot the delimitation of the crhomosomes
  
  # abline(v = chr.border, col = 1, lwd = 2)
  segments(x0=chr.border, y0=0, x1=chr.border, y1=max.y, lwd = 2)
  
  
  # add the number of QTL detected infos
  
  points(x = points.data[, 1], y = points.data[, 2], type = "h")
  
  
}