##################
# subset.mppData #
##################

#' Subset \code{mppData} object
#' 
#' Pull out a specified set of markers and/or genotypes from a \code{mppData}
#' object.
#' 
#' @param x An object of class \code{mppData}.
#' 
#' @param mk.list Optional \code{character} vector, \code{numeric} position
#' vector or \code{logical} vector representing marker to keep. Default = NULL.
#' 
#' 
#' @param gen.list Optional \code{character} vector, \code{numeric} position
#' vector or \code{logical} vector representing genotypes to keep.
#' Default = NULL.
#' 
#' @param ... Ignored.
#' 
#' @return Return:
#' 
#' The mppData object but with only the specified subset of data.
#' 
#' 
#' @author Vincent Garin
#' 
#' @examples
#' 
#' ### Marker subset
#' 
#' data(mppData)
#' 
#' # Random selection of markers
#' mk.list <-  sample(mppData$map[, 1], 50)
#' mppData_sub <- subset(x = mppData, mk.list = mk.list)
#' 
#' # Selection of chromosome 1 marker
#' mk.list <-  (mppData$map[, 2] == 1)
#' mppData_sub <- subset(x = mppData, mk.list = mk.list)
#' 
#' ### Genotype subset
#' 
#' # Random selection of genotypes
#' gen.list <-  sample(mppData$geno.id, 200)
#' mppData_sub <- subset(x = mppData, gen.list = gen.list)
#' 
#' # Selection of genotype from cross 2 and 5
#' crosses <- unique(mppData$cross.ind)
#' gen.list <-  mppData$geno.id[mppData$cross.ind %in% crosses[c(2, 5)]]
#' mppData_sub <- subset(x = mppData, gen.list = gen.list)
#' 
#' ### Marker and genotype subset
#' 
#' mk.list <-  sample(mppData$map[, 1], 50)
#' gen.list <-  sample(mppData$geno.id, 200)
#' mppData_sub <- subset(x = mppData, mk.list = mk.list,
#' gen.list = gen.list)
#' 
#' @export
#' 



subset.mppData <- function(x, mk.list = NULL, gen.list = NULL, ...) {
  
  # 1. Check data.format and arguments
  ####################################
  
  check_mppData(mppData = x)
  
  # user must at least specify one argument (mk.list or gen.list)
  
  if(is.null(mk.list) && is.null(gen.list)){
    
    stop("You must specify either mk.list or gen.list.")
    
  }
  
  
  # test the format of the marker and genotype list given by the user
  
  if (!is.null(mk.list)){
    
    # test format marker list
    
    if (!(is.character(mk.list) || is.logical(mk.list) || is.numeric(mk.list))) {
      
      stop("mk.list must be a character, logical or numeric vector.")    
      
    }
    
    if(is.logical(mk.list)){
      
      if(length(mk.list) != dim(x$map)[1])
        
        stop("The mk.list does not have the same length as the map")
      
    }
    
    # Convert different type of marker list into a character vector.
    
    if (is.logical(mk.list) || is.numeric(mk.list)) {
      
      mk.list <- x$map[mk.list, 1]
      
    } 
    
    
  }
  
  if (!is.null(gen.list)){
    
    # test format marker list
    
    if (!(is.character(gen.list) || is.logical(gen.list) || is.numeric(gen.list))) {
      
      stop("gen.list must be a character, logical or numeric vector.")    
      
    }
    
    if(is.logical(gen.list)){
      
      if(length(gen.list) != length(x$geno.id))
        
        stop("gen.list does not have the same length as the genotypes list")
      
    }
    
    # Convert different type of marker list into a character vector. Then in a
    # logical list to keep the same order.
    
    if (is.logical(gen.list) || is.numeric(gen.list)) {
      
      gen.list <- x$geno.id[gen.list]
      
      geno.ind <- x$geno.id %in% gen.list 
      
    } else if(is.character(gen.list)) {
      
      geno.ind <- x$geno.id %in% gen.list
      
    }
    
  }
  
  # 2. subset by markers
  ######################
  
  if(!is.null(mk.list)){
    
    # geno.IBS
    
    x$geno.IBS <- x$geno.IBS[, x$map[, 1] %in% mk.list,
                                 drop = FALSE]
    
    # geno.IBD
    
    for (i in 1:length(x$geno.IBD$geno)) {
      
      # get the list of marker of the considered chromosome
      
      chr.mk.names <- attr(x$geno.IBD$geno[[i]]$prob, "dimnames")[[2]]
      
      # indicate which marker of chromosome i is in the reduced list
      
      indicator <- chr.mk.names %in% mk.list
      
      x$geno.IBD$geno[[i]]$prob <- x$geno.IBD$geno[[i]]$prob[, indicator, ]
      
      
    }
    
    # allele.ref
    
    x$allele.ref <- x$allele.ref[, x$map[, 1] %in% mk.list,
                                             drop = FALSE]
    
    # geno.par
    
    x$geno.par <- x$geno.par[(x$geno.par[, 1] %in% mk.list), ,
                                         drop = FALSE]
    
    # par.clu
    
    if(!is.null(x$par.clu)){
      
      x$par.clu <- x$par.clu[rownames(x$par.clu) %in% mk.list, ,
                             drop = FALSE]
      
    }
    
    # map
    
    x$map <- x$map[(x$map[, 1] %in% mk.list), , drop = FALSE]
    
    # recalculate the position indicators
    
    chr.ind <- factor(x = x$map[, 2], levels = unique(x$map[, 2]))
    
    x$map[, 3] <- sequence(table(chr.ind))
    
  }
  
  # 3. subset by genotypes
  ########################
  
  if(!is.null(gen.list)){
    
    # geno.IBS
    
    x$geno.IBS <- x$geno.IBS[geno.ind, ]
    
    # geno.IBD
    
    for (i in 1:length(x$geno.IBD$geno)) {
      
      x$geno.IBD$geno[[i]]$prob <- x$geno.IBD$geno[[i]]$prob[geno.ind, , ]
      
    }
    
    # pheno
    
    x$pheno <- x$pheno[geno.ind, , drop = FALSE]
    
    # geno.id
    
    x$geno.id <- x$geno.id[geno.ind]
    
    # ped.mat
    
    ped.temp <- as.matrix(x$ped.mat)
    
    ped.mat.found <- ped.temp[ped.temp[, 1] == "founder", , drop = FALSE]
    ped.mat.off <- ped.temp[ped.temp[, 1] == "offspring", , drop = FALSE]
    ped.mat.off <- ped.mat.off[geno.ind, , drop = FALSE]
    
    found.sub <- unique(c(ped.mat.off[, 3], ped.mat.off[, 4]))
    ped.mat.found <- ped.mat.found[ped.mat.found[, 2] %in% found.sub, ,
                                   drop = FALSE]
    
    pedigree.new <- rbind(ped.mat.found, ped.mat.off)
    
    x$ped.mat <- data.frame(pedigree.new, stringsAsFactors = FALSE)
    
    # cross.ind
    
    x$cross.ind <- x$cross.ind[geno.ind]
    
  }
  
  # 4. adapt the different elements
  #################################
  
  # review the par.per.cross, parents, n.cr and n.par objects
  
  list.cr <- unique(x$cross.ind)
  
  ppc <- x$par.per.cross
  ppc <- ppc[ppc[, 1] %in% list.cr, , drop = FALSE]
  x$parents <- union(ppc[, 2], ppc[, 3])
  x$n.cr <- length(list.cr)
  x$n.par <- length(x$parents)
  x$par.per.cross <- ppc
  
  # modify the geno.par argument. Remove the parent that are not used anymore
  
    new_par <- c("mk.names", "chr", "pos.ind", "pos.cM", x$parents)
    
    x$geno.par <- x$geno.par[, colnames(x$geno.par) %in% new_par ,
                                         drop = FALSE]
    
    # check the par.clu
    
    if(!is.null(x$par.clu)){
      
      par.clu <- x$par.clu
      par.clu <- par.clu[, x$parents]
      
      par.clu <- parent_clusterCheck(par.clu = par.clu)[[1]]
      
      x$par.clu <- par.clu
      
    }
  
    class(x) <- c("mppData", "list")
    
    return(x)

  
}