#################
# cross_ABH_het #
#################


# ABH assignment per cross with heterozygous or missing parents
# 
# Transform offspring genotype scores into A, B, H or NA (missing) according
# to the scores of parents 1 and 2 of each cross when parents have markers with
# heterozygous scores.
# 
# For marker without heterozygous parent scrores, the assignment rules is the
# same as in \code{\link{cross_ABH}}. If a parent score is heterozygous or
# missing, the assignment rules are the following. If the two parents are
# heterozygous or one parent is heterozygous and the other missing, the
# offspring get NA since the parental origin can not be inferred with certainty.
# If one parent is heterozygous or missing and the second parent is
# homozygous, the function looks at offspring segregating pattern to try to
# infer which allele was transmitted by the heterozygous parent. If this is
# possible we consider the heteroxzygous parent as homozygous for the
# transmitted allele and use this allele score for ABH assignment.
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
# geno.ABH <- cross_ABH_het(par.sc = par.sc, off.sc = off.sc,
#                           cross.ind = cross.ind,
#                           par.per.cross = par.per.cross)
# 
# @export



cross_ABH_het <- function(par.sc, off.sc, cross.ind, par.per.cross) {
  
  # 1. check data format
  ######################
  
  check.cr.ABH(par.sc = par.sc, off.sc = off.sc, cross.ind = cross.ind, 
               par.per.cross = par.per.cross)
  
  
  # 2. Determine the reference genotype score
  ###########################################
  
  allele.sc <- function(x) {
    
    x <- x[!is.na(x)] # remove missing values
    
    if (length(x) == 0){ # Only missing values
      
      return(c(NA, NA, NA, NA))
      
    } else {
      
      # split all marker score into single allele scores
      
      alleles <- c(vapply(X = x, FUN = function(x) unlist(strsplit(x, split = "")),
                          FUN.VALUE = character(2)))
      
      alleles.sum <- table(alleles)
      
      if (length(alleles.sum) == 1) { # monomorphic case
        
        return(c(NA, NA, NA, NA))
        
      } else {
        
        all.sc <- unlist(attr(sort(alleles.sum, decreasing = TRUE), "dimnames"))
        all.sc <- c(paste0(all.sc[1], all.sc[1]), paste0(all.sc[2], all.sc[2]),
                    paste0(all.sc[1], all.sc[2]), paste0(all.sc[2], all.sc[1]))
        
        return(all.sc)
        
      }
      
    } 
    
  }
  
  allele.ref <- apply(X = rbind(par.sc, off.sc), MARGIN = 2,
                      FUN = allele.sc)[1:2, ]
  
  
  # 3. Iteration along the different crosses 
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
    
    ### 3.1 Segregation pattern of the parents
    
    par.segretation <- function(x){
      
      # sub-function to test homozygosity
      
      al.sc.homo <- function(x) {
        
        strsplit(x, split = "")[[1]][1] == strsplit(x, split = "")[[1]][2]
        
      }
      
      miss <- is.na(x)
      
      if(sum(miss) == 2){return(1)
        
      } else if (sum(miss) == 1) {
        
        # test if the non-missing value is hetero or homozygous
        
        if(al.sc.homo(x = x[!miss])){return(2)} else {return(1)}
        
      }
      
      all.sc <- unique(x)
      
      if(length(all.sc) == 1) { return(1)
        
      } else {
        
        # test if homo/homo (AA/TT) or homo/hetero (AA/AT)
        
        test <- unlist(lapply(X = x, FUN = al.sc.homo))
        
        if (sum(!test) == 1) {return(2)} else{return(3)}
        
      }
      
      
    }
    
    par.seg <- apply(X = parents.gen, MARGIN = 2, FUN = par.segretation)
    
    ############# End parent segregation determination
    

    
    ### 3.2 Fill ABH score per cross (sub-cross)
    
    off.gen <- subset(x = off.sc, subset = cross.ind == crosses[k], drop = FALSE)
    
    empty.mat <- matrix(0, dim(off.gen)[1], dim(off.gen)[2])  # empty matrix
    
    # we scan per column (per marker position)
    
    for (i in 1:length(par.seg)) {
      
      ### 3.2.1: Cases where parent orignin can not identified with certainty:
      # 1) two allele missing; 2) two parents heterozygous; 3) two parents are
      # monomorphic; 4) or one parent is missing and the other is heterozygous
      
      if (par.seg[i] == 1){
        
        empty.mat[, i] <- NA
        
      ### 3.2.2: Situations where the allele of the heterogeneous or missing
      # parent can be inferred looking at the segregation pattern of the
      # offspring.
        
      } else if (par.seg[i] == 2){
        
        # test if fully missing
        
        if(sum(!is.na(off.gen[, i])) == 0){
          
          empty.mat[, i] <- NA
          
        } else {
          
          # test if there is segregation
          
          if(length(unique(off.gen[, i])) == 1){ # no segregation
            
            empty.mat[, i] <- NA
            
          } else { # segregation
            
            parents.gen[, i]
            
            # find the parent that is heterozygous or missing
            
            al.sc.homo <- function(x) {
              
              strsplit(x, split = "")[[1]][1] == strsplit(x, split = "")[[1]][2]
              
            }
            
            test.homo <- unlist(lapply(X = parents.gen[, i], FUN =  al.sc.homo))
            
            ind <- (is.na (parents.gen[, i]) | !test.homo)
            homo.sc <- parents.gen[!ind, i]
            opp.sc <- allele.ref[, i][allele.ref[, i] != homo.sc]
            parents.gen[ind, i] <- opp.sc
            
            # then assignement according to the new references
            
            empty.mat[is.na(off.gen[, i]), i] <- NA
            
            # test if it looks like parent A (or 1)
            
            empty.mat[off.gen[, i] %in% parents.gen[1, i], i] <- "A"
            
            # test if it looks like parent B (or 2)
            
            empty.mat[off.gen[, i] %in% parents.gen[2, i], i] <- "B"
            
            # put the rest heterozygous
            
            empty.mat[empty.mat[, i] %in% "0", i] <- "H"
            
          }
          
        }
        
        ### 3.2.3: Perfectly identifiable situation where the two parents are
        # homozygous with the opposite allele score.
        
        
      } else if (par.seg[i] == 3){
        
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
