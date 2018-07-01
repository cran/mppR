########################
# design_connectivity #
########################

#'
#' Connected parts of a MPP design
#' 
#' Determine the connected parts of a MPP design using the method of Weeks and
#' Williams (1964) and the package igraph.
#' 
#' @param par_per_cross Three columns \code{character matrix} specifying :
#' 1) the cross indicators ; 2) the parents 1 identifiers of the crosses;
#' 3) the parents 2 identifiers of the crosses.
#' 
#' @param plot_des \code{Logical} value indicating if connected part of the
#' design should be plotted. Default = TRUE.
#' 
#' @param output_loc Path where the plot of the design will be saved if the
#' argument is given. Default = NULL.
#' 
#' @return
#' 
#' Return a list with each element representing one connected part of the design
#' and the list of parents contained in this part.
#' 
#' If \code{plot_des = TRUE} and \code{output_loc} has been specified. A plot
#' of the graph (con_plot.pdf) will be saved at the specified location. 
#' 
#' @author Vincent Garin
#' 
#' @references 
#' 
#' Weeks, D. L., & Williams, D. R. (1964). A note on the determination of
#' connectedness in an N-way cross classification. Technometrics, 6(3), 319-324.
#' 
#' @examples
#' 
#' data(mppData)
#' 
#' par_per_cross <- mppData$par.per.cross
#' 
#' con.part <- design_connectivity(par_per_cross)
#' 
#' @export
#' 

design_connectivity <- function(par_per_cross, plot_des = TRUE,
                                output_loc = NULL){
  
  # 1. Test format of the arguments
  #################################
  
  if(!is.matrix(par_per_cross)){
    
    stop("'par_per_cross' is not a matrix")
  }
  
  if(!is.character(par_per_cross)){
    
    stop("'par_per_cross' is not a character matrix")
    
  }
  
  if(!is.null(output_loc)){
    
    if(!file.exists(output_loc)){
      
      stop("'output_loc' is not a valid path")
      
    }
    
  }
  
  if((!is.null(output_loc) && !plot_des)){ plot_des <- TRUE }
  
  # 2. Determine design conectedness and plot
  ###########################################
  
  vertices <- apply(X = par_per_cross[, 2:3], MARGIN = 1, FUN = function(x) x)
  
  vertices <- c(vertices)
  
  g <- graph(vertices)
  
  if(plot_des) {
    
    plot(g)
    
    if(!is.null(output_loc)){
      
      pdf(file.path(output_loc, "con_plot.pdf"), height = 10, width = 16)
      
      print(plot(g))
      
      dev.off()
      
    }
    
  }
  
  
  
  clu <- components(g)
  
  grp <- groups(clu)
  
  return(grp)
  
  
}