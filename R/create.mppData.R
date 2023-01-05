##################
# create.mppData #
##################

#' Create a multi-parent population data object
#' 
#' This function combines all raw data sources in a single data object of class
#' \code{mppData}.
#' 
#' @param geno.off \code{Character} marker score \code{matrix} of the offspring
#' with genotypes as row and markers as column.
#' \strong{Rows names must be the offspring genotypes identifiers similar to
#' the one used in \code{pheno}. The columns names must be the marker names
#' similar to the one used in \code{map}. Marker scores must be coded using one
#' letter per allele. For example, AA, CC, GG, etc. Missing values must be coded
#' \code{NA}.}
#' 
#' @param geno.par \code{Character} marker score \code{matrix} of the parents
#' with genotypes as row and markers as column.
#' \strong{Rows names must be the parents genotypes identifiers similar to
#' the one used in \code{par.per.cross}. The columns names must be the marker
#' names similar to the one used in \code{map}. Marker scores must be coded
#' using one letter per allele. For example, AA, CC, GG, etc. Missing values
#' must be coded \code{NA}.}
#' 
#' @param map Three columns \code{data.frame} with: 1) \code{character} marker
#' identifiers; 2) \code{numeric} chromosome; 3) \code{numeric} positions in
#' centi-Morgan.\strong{The marker identifiers must be identical to the column
#' names of the maker matrices (\code{geno.off} and \code{geno.par}).
#' The chromosome identifiers must start by 1 and increase by 1 unit,
#' e.g. 1, 2, 3, ...}
#' 
#' @param pheno \code{Numeric matrix} with one column per trait and rownames
#' as genotypes identifiers. \strong{The genotypes identifiers must be identical
#' to the rownames of the offspring marker matrix (\code{geno.off}).}
#' 
#' @param cross.ind \code{Character} vector indicating to which cross each
#' genotype belongs.
#' 
#' @param par.per.cross Three columns \code{Character matrix} specifying :
#' 1) the cross indicators; 2) the parents 1 identifiers
#' of the crosses; 3) the parents 2 identifiers of the crosses.
#' \strong{The list of crosses must contain the same cross indicators as in
#' \code{cross.ind} and they must appear in the same order.
#' The list of parent identifiers must be the same to the rownames of
#' the argument \code{geno.par}}.
#' 
#' @return 
#' 
#' a \code{list} of class \code{mppData} which contains the following elements
#' 
#' \item{geno.off}{\code{Matrix} with offspring marker scores.}
#' 
#' \item{geno.par}{\code{Matrix} with parents marker scores.}
#' 
#' \item{pheno}{\code{Matrix} with phenotypic trait values.}
#' 
#' \item{map}{\code{Data.frame} with genetic marker information.}
#' 
#' \item{cross.ind}{Cross indicator.}
#' 
#' \item{par.per.cross}{\code{Character matrix} information about cross and
#' the parents of the crosses.}
#' 
#' The \code{list} also contain other arguments that will be filled later in
#' the data processing.
#' 
#' @author Vincent Garin
#' 
#' @examples
#' 
#' data(USNAM_geno)
#' data(USNAM_map)
#' data(USNAM_pheno)
#' 
#' geno.off <- USNAM_geno[7:506, ]
#' geno.par <- USNAM_geno[1:6, ]
#' map <- USNAM_map
#' pheno <-  USNAM_pheno
#' cross.ind <- substr(rownames(pheno), 1, 4)
#' par.per.cross <- cbind(unique(cross.ind), rep("B73", 5),
#'                        rownames(geno.par)[2:6])
#' 
#' mppData <- create.mppData(geno.off = geno.off, geno.par = geno.par,
#'                           map = map, pheno = pheno, cross.ind = cross.ind,
#'                           par.per.cross = par.per.cross)
#'                           
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @import igraph
#' @import nlme
#' @import parallel
#' @import qtl
#' @import methods
#' @import utils
#' @importFrom dplyr arrange group_by mutate summarise
#' @importFrom Matrix Matrix t
#' @importFrom stats anova as.formula coef coefficients complete.cases cor df.residual lm
#' @importFrom stats model.matrix na.omit pchisq pt quantile runif vcov                           
#' 
#' @export
#' 

create.mppData <- function(geno.off = NULL, geno.par = NULL, map = NULL,
                           pheno = NULL, cross.ind = NULL,
                           par.per.cross = NULL){
  
  # 1. Check the correct format of the raw data
  #############################################
  
  # Test if one object is missing
  
  if(is.null(geno.off)){stop("'geno.off' is not provided")}
  
  if(is.null(geno.par)){stop("'geno.par' is not provided")}
  
  if(is.null(map)){stop("'map' is not provided")}
  
  if(is.null(pheno)){stop("'pheno' is not provided")}
  
  if(is.null(cross.ind)){stop("'cross.ind' is not provided")}
  
  if(is.null(par.per.cross)){stop("'par.per.cross' is not provided")}
  
  # test format of the marker matrices (geno.off, geno.par)
  
  if(!is.matrix(geno.off)){ stop("'geno.off' is not a matrix") }
  
  if(!is.character(geno.off)){
    
    stop("'geno.off' is not character")
    
  }
  
  if(!is.matrix(geno.par)){ stop("'geno.par' is not a matrix") }
  
  if(!is.character(geno.par)){
    
    stop("'geno.par' is not character")
    
  }
  
  # test if the marker identifiers are the same in the map and in the marker
  # matrices
  
  if(!identical(colnames(geno.off), map[, 1])){
    
    stop("the marker identifiers in 'geno.off' and 'map' are not the same")
    
  }
  
  if(!identical(colnames(geno.par), map[, 1])){
    
    stop("the marker identifiers in 'geno.par' and 'map' are not the same")
    
  }
  
  # test map format
  
  if(!is.character(map[, 1])){
    
    stop("the marker identifier in 'map' must be character")
    
  }
  
  if(!is.numeric(map[, 2])){
    
    stop("the chromosome in 'map' must be numeric")
    
  }
  
  # test that chr. id goes from 1 to n_cr increasing by one
  
  chr_id <- unique(map[, 2])
  chr_ref <- 1:length(chr_id)
  
  if(sum(chr_id - chr_ref) != 0){
    
    stop("the chromosome identifier in 'map' must start by 1 and increase by 1 unit, e.g. 1, 2, 3, ...")
    
  }
  
  if(!is.numeric(map[, 3])){
    
    stop("the genetic positions in 'map' must be numeric")
    
  }
  
  # test phenotype values
  
  if(!is.matrix(pheno)){stop(" 'pheno' is not a matrix")}
  
  if(!is.numeric(pheno)){stop(" 'pheno' is not numeric")}
  
  # test if the list of genotype is the same between the offspring marker matrix
  # and the phenotypic values.
  
  if(!identical(rownames(geno.off), rownames(pheno))){
    
    stop("the genotype identifiers of 'geno.off' and 'pheno' are not identical")
               
    
  }
  
  # Add trait title is not provided
  
  if(is.null(colnames(pheno))){
    
    trait_names <- paste0('trait', 1:dim(pheno)[2])
    colnames(pheno) <- trait_names
    
  }
  
  # length cross.ind same as list of genotypes
  
  if(length(cross.ind) != dim(geno.off)[1]){
    
    stop("'cross.ind' length is not equal to the number of genotype in 'geno.off'")
  }
  
  # test par.per.cross
  
  if(!is.matrix(par.per.cross)){
    
    stop("'par.per.cross' is not a matrix")
    
  }
  
  if(!is.character(par.per.cross)){
    
    stop("'par.per.cross' is not character")
    
  }
  
  # remove the eventual rownames of par.per.cross
  
  if(!is.null(rownames(par.per.cross))){
    
    rownames(par.per.cross) <- NULL
    
  }
  
  if (!identical(unique(cross.ind), par.per.cross[, 1])){
    
    stop("the cross identifiers in 'cross.ind' and 'par.per.cross' are not identical")
    
  }
  
  # test the similarity of parents list between par.per.cross and
  # rownames(geno.par)
  
  parents <- union(par.per.cross[, 2], par.per.cross[, 3])
  
  if(sum(!(parents %in% rownames(geno.par)))>0){
    
    list.par <- parents[!(parents %in% rownames(geno.par))]
    pbpar <- paste(list.par, collapse = ", ")
    
    message <- sprintf(ngettext(length(list.par),
                                "parent %s is used in 'par.per.cross' but not in 'geno.par'",
                                "parents %s are used in 'par.per.cross' but not in 'geno.par'"),
                       pbpar)
    
    stop(message)
    
  }
  
  # check if the marker are grouped by chromosome and are in increasing position
  
  map <- map[order(map[, 2]), ]
  
  map_temp <- c()
  
  chr_id <- unique(map[, 2])
  
  for(i in 1:length(chr_id)){
    
    map_i <- map[map[, 2] == chr_id[i], ]
    map_i_ord <- map_i[order(map_i[, 3]), ]
    
    map_temp <- rbind(map_temp, map_i_ord)
    
  }
  
  map <- map_temp
  
  # Re-order the marker matrix in case of
  
  geno.par <- geno.par[, map[, 1]]
  geno.off <- geno.off[, map[, 1]]
  
  # Check the number of connected parts
  
  con_part <- design_connectivity(par.per.cross, plot_des = FALSE)
  n_con <- length(con_part)
  
  n_geno <- dim(geno.off)[1]
  n_par <- dim(geno.par)[1]
  n_cr <- dim(par.per.cross)[1]
  n_pheno <- dim(pheno)[2]
  
  # Possibility to check here the necessity to have at least 2 crosses and 3
  # different parents.
  
  ####### end check format
  
  mppData <- list(geno.off = geno.off, geno.IBS = NULL, geno.IBD = NULL,
                  geno.id = NULL, ped.mat = NULL, allele.ref = NULL,
                  geno.par = geno.par, geno.par.clu = NULL, par.clu = NULL,
                  pheno = pheno, map = map, haplo.map = NULL,
                  cross.ind = cross.ind, par.per.cross = par.per.cross,
                  type = NULL, parents = NULL, n.cr = NULL, n.par = NULL,
                  n.zigo = NULL, rem.mk = NULL, rem.gen = NULL,
                  status = 'init')
  
  class(mppData) <- c("mppData", "list")
  
  # Return a short message that everything was successful
  
  cat("\n")
  cat("mppData object created!\n")
  cat("\n")
  cat(paste(n_geno, 'genotypes\n'))
  cat(paste(n_cr, 'crosses\n'))
  cat(paste(n_par, 'parents\n'))
  cat(paste(n_pheno, 'phenotype(s)\n'))
  
  if(n_con == 1){
    
    cat('1 connected part\n')
    
  } else {message(paste(n_con, 'independent connected parts\n'))}
  
  
  return(mppData)
  
}