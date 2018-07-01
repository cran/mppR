################
# check.cr.ABH #
################

# function to check the data that are introduced in cross.ABH function

check.cr.ABH <- function(par.sc, off.sc, cross.ind, par.per.cross){
  
  # check data
  ############
  
  # check if the parent and offspring matrix are characters
  
  if(!is.character(par.sc)){
    
    stop("'par.sc' argument is not character")
    
  }
  
  if(!is.character(off.sc)){
    
    stop("'off.sc' is not character")
    
  }
  
  # check if the cross.ind is character
  
  if(!is.character(cross.ind)){
    
    stop("'cross.ind' is not character")
    
  }
  
  # check if the cross.indicator as the same length as the list of genotype
  
  if(length(cross.ind) != dim(off.sc)[1]){
    
    stop("'cross.ind' length is not equal to the number of genotype in 'off.sc'")
    
  }
  
  # check if the cross.ind is character
  
  if(!is.character(par.per.cross)){
    
    stop("'par.per.cross' is not character")
    
  }
  
  # remove the eventual rownames of par.per.cross
  
  if(!is.null(rownames(par.per.cross))){
    
    rownames(par.per.cross) <- NULL
    
  }
  
  
  
  if (!identical(unique(cross.ind), par.per.cross[, 1])){
    
    stop("the cross identifiers in 'cross.ind' and in 'par.per.cross' differ")
    
  }
  
  
  # check if the parent list is the same in the par.per.cross
  # and in the parents score rownames
  
  parents <- union(par.per.cross[,2], par.per.cross[,3])
  
  if(sum(!(parents %in% rownames(par.sc)))>0){
    
    list.par <- parents[!(parents %in% rownames(par.sc))]
    pbpar <- paste(list.par, collapse = ", ")
    
    message <- sprintf(ngettext(length(list.par),
                                "parent %s is used in 'par.per.cross' but not in 'par.sc'",
                                "parents %s are used in 'par.per.cross' but not in 'par.sc'"),
                       pbpar)
    
    stop(message)
    
  }
  
  if(sum(!(rownames(par.sc) %in% parents))>0){
    
    list.par <- rownames(par.sc)[!(rownames(par.sc) %in% parents)]
    pbpar <- paste(list.par, collapse = ", ")
    
    message <- sprintf(ngettext(length(list.par),
                                "parent %s is used in 'par.sc' but not in 'parents'",
                                "parents %s are used in 'par.sc' but not in 'parents'"),
                       pbpar)
    
    stop(message)
    
  }
  
}