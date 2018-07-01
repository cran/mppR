##############
# cross_ABH #
##############


# ABH assignment per cross
# 
# Transform offspring genotype scores into A, B, H or NA (missing) according
# to the scores of parents 1 and 2 of each cross.
# 
# The function transforms offspring genotype data of each cross. The function
# takes successively the parents of the different cross as
# reference and assign the following scores: "A" if the offspring score is
# equivalent to parent 1; "B" if it is equivalent to parent 2; "H" if it is
# heterozygous. The function attributes NA for missing when: 1) the offspring
# score is missing; 2) the two parents have the same score; or
# 3) when at least one parental score is missing.
# 
# @param par.sc \code{Character} marker scores \code{matrix} of the
# parents of the different cross. \strong{The rownames
# must be the parents identifiers and correspond to the one used in argument
# \code{par.par.cross}. Parent scores must be homozygous.
# Missing value must be coded NA.}
# 
# @param off.sc \code{Character} marker scores \code{matrix} of the
# offspring. \strong{The possible values must be: homozygous like parent 1 or 2,
# heterozygous or missing. Missing values must be coded NA.}
# 
# @param cross.ind \code{Character} vector with the same length as the number
# of offspring genotypes which specifies to which cross each offspring genotype
# belongs.
# 
# @param par.per.cross Three columns \code{Character matrix} specifying :
# 1) the cross indicators (\strong{The cross indicators must be  similar to
# the one used in \code{cross.ind} and appear in the same order}); 2) the
# parents 1 identifiers of the crosses; 3) the parents 2 identifiers of the
# crosses. \strong{The list of parent identifiers must be similar to the
# rownames of the argument \code{par.sc}}.
# 
# @return Return:
# 
# \item{geno.ABH}{\code{Character matrix} with ABH coding for the different
# genotypes of the population}
# 
# @author Vincent Garin
# 
# @examples
# 
# # Genotypes
# data("USNAM_geno")
# geno <- USNAM_geno
# 
# # Remove markers for which parents are monomorphic
# par.mono <- QC_MAF(mk.mat = geno[1:6, ])
# geno <- geno[, -which(par.mono == 0)]
# 
# # Remove markers for which at least one parent is heterozygous
# par.het <- QC_hetero(mk.mat = geno[1:6, ])
# geno <- geno[, -which(par.het != 0)]
# 
# # Parents scores
# par.sc <- geno[1:6, ]
# 
# # Offspring scores
# off.sc <- geno[7:506, ]
# 
# # Cross indicator
# cross.ind <- substr(rownames(geno)[7:506], 1, 4)
# 
# # Parent of the crosses matrix
# par.per.cross <- cbind(unique(cross.ind), rep("B73",5), rownames(par.sc)[-1])
# 
# geno.ABH <- cross_ABH(par.sc = par.sc, off.sc = off.sc, cross.ind = cross.ind,
#                       par.per.cross = par.per.cross)
#
# 
# @export
# 


cross_ABH <- function(par.sc, off.sc, cross.ind, par.per.cross) {
  
  # 1. check data format
  ######################
  
  check.cr.ABH(par.sc = par.sc, off.sc = off.sc, cross.ind = cross.ind, 
               par.per.cross = par.per.cross)
  
  
  # 2. function to test if parent are different (segregate) 
  #########################################################
  
  seg.test <- function(x){
    
    if(NA %in% x){ res <- FALSE } else {
      
      if(x[1] != x[2]) { res <-  TRUE } else { res <- FALSE  }
    }
    return(res)
  }
  
  # 3. function to test if the offspring scores are monomorphic
  #############################################################
  
  test.mono <- function(x){ length(unique(x[!is.na(x)])) == 1 }
  
  # 4. Iteration along the different crosses 
  ##########################################
  
  geno.ABH <- c() # empty matrix to store results
  
  crosses <- unique(cross.ind)
  
  for (k in seq_along(crosses)) {
    
    
    # select the parents and the offsprings that correspond the considered
    # cross
    
    parents.line <- par.per.cross[k, 2:3]
    
    parent1.gen <- par.sc[which(rownames(par.sc) == parents.line[1]), ]
    parent2.gen <- par.sc[which(rownames(par.sc) == parents.line[2]), ]
    parents.gen <- rbind(parent1.gen, parent2.gen)
    
    # test if the parent have a different marker score (segregate)
    
    mk.inf <- apply(X = parents.gen, MARGIN =  2, FUN = seg.test)
    
    off.gen <- subset(x = off.sc, subset = cross.ind == crosses[k], drop = FALSE)
    
    off.mono <- apply(X = off.gen, MARGIN = 2, FUN = test.mono)
    
    empty.mat <- matrix(0, dim(off.gen)[1], dim(off.gen)[2])  # empty matrix
    
    # we scan per column (per marker position)
    
    for (i in 1:dim(off.gen)[2]) {
      
      # we test if the marker is informative (segregate within the parents)
      
      if ((!mk.inf[i]) || off.mono[i]) {
        
        # the marker is not informative (parent monomorphic or missing) or the
        # offspring scores are monomorphic so all genotypes are put as missing
        
        empty.mat[, i] <- NA
        
      } else { 
        
        # if the marker segregates, then we look at each offspring position to
        # see from which parent it comes from.
        
        # test if the position are missing
        
        empty.mat[is.na(off.gen[, i]), i] <- NA
        
        # test if it looks like parent A (or 1)
        
        empty.mat[off.gen[, i] %in% parents.gen[1, i], i] <- "A"
        
        # test if it looks like parent B (or 2)
        
        empty.mat[off.gen[, i] %in% parents.gen[2, i], i] <- "B"
        
        # put the rest heterozygous
        
        empty.mat[empty.mat[, i] %in% "0", i] <- "H"
        
      }
      
    }
    
    # append the different crosses sucessively
    
    geno.ABH <- rbind(geno.ABH, empty.mat)
    
  }
  
  colnames(geno.ABH) <- colnames(off.sc)
  rownames(geno.ABH) <- rownames(off.sc)
  
  return(geno.ABH)
  
}