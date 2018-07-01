##################
# QC_matchMarker #
##################

# Match markers in the genotype matrix and the map
# 
# Determine which markers are common between the genotypes marker matrix and
# makers present in the map. Return the marker matrix and the map that
# correspond to these markers.
# 
# 
# @param mk.mat Marker score \code{matrix} or \code{data.frame} with genotypes
# as row and markers as column. \strong{Marker names must be column names of the
# matrix.}
# 
# @param map \code{Data.frame} \strong{with first column containing the marker
# names}.
# 
# @return Return:
# 
# \code{List} containing the following objects:
# 
# \item{new.mk.mat}{Marker \code{matrix} with only the common list of
# marker between the marker matrix and the map in the same order as in the map.}
# 
# \item{new.map}{Map with only the common markers list.}
# 
# \item{mk.mat.only}{Markers present in the marker matrix but not in the map.}
# 
# \item{map.only}{Markers present in the map but not in the marker matrix.}
# 
# @author Vincent Garin
# 
# @seealso \code{\link{QC_matchGeno}}
# 
# @examples
# 
# data(USNAM_map)
# data(USNAM_geno)
# 
# # Remove one marker from the map
# map.red <- USNAM_map[-1,]
# 
# match <- QC_matchMarker(mk.mat = USNAM_geno, map = map.red)
# 
# # Marker only in the map
# match$map.only
# 
# # Marker only in the marker matrix
# match$mk.mat.only
# 
# # Get the new map and marker matrix
# new.mk.mat <- match$new.mk.mat
# new.map <- match$new.map
# 
# rm(match)
# 
# @export
# 


QC_matchMarker <- function(mk.mat, map) {
  
  # identify common marker in marker matrix and in map (intersection A n B)
  
  inter.mk.map <- intersect(map[, 1], colnames(mk.mat))
  
  # marker only in the genotype matrix (A \B)
  
  mk.mat.only <- setdiff(colnames(mk.mat), map[, 1])
  
  # marker only in the map (B \A)
  
  map.only <- setdiff(map[, 1], colnames(mk.mat))
  
  # subset the map and the matrix
  
  new.mk.mat <- subset(x = mk.mat, select = inter.mk.map, drop = FALSE)
  
  new.map <- subset(x = map, subset = map[, 1] %in% inter.mk.map, drop = FALSE)
  
  # make sure that the marker matrix is in the same order as the map
  
  new.mk.mat <- new.mk.mat[, as.character(new.map[, 1])]
  
  
  return(list(new.mk.mat = new.mk.mat, new.map = new.map,
              mk.mat.only = mk.mat.only, map.only = map.only))
  
}