#############
# check_IBD #
#############

# Function to check the argument provided to the function IBS.mppData which
# convert the genotype data into 0, 1, 2 format and do an optional
# genotype imputation.

check_IBD <- function(mppData, het.miss.par, subcross.ind, par.per.subcross,
                      type, F.gen, BC.gen, type.mating, map.function){
  
  # 1. check mppData format
  #########################
  
  if(!is_mppData(mppData)){
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  # test if correct step in the mppData processing
  
  if(mppData$status != 'IBS'){
    
    stop("you have to process 'mppData' in a strict order: ",
               "create.mppData, QC.mppData, IBS.mppData, IBD.mppData, ",
               "parent_cluster.mppData. You can only use IBD.mppData ",
               "after create.mppData, QC.mppData, and IBS.mppData")
    
  }
  
  geno.par <- mppData$geno.par
  
  # 2. check ABH coding arguments
  ###############################
  
  # Test if there are some heterozygous parents and the user want to make the
  # ABH assignement without specifying het_miss_par = TRUE
  
  if(!het.miss.par){
    
    # test if there are some heterozygous parents
    
    par.het <- QC_hetero(mk.mat = geno.par)
    
    het.mk <- which(par.het != 0)
    
    if(length(het.mk) > 0){
      
      pb_mk <- paste(colnames(geno.par)[het.mk], collapse = ", ")
      
      stop("The following markers: ", pb_mk,
           " are heterozygous for at least one parent. In order ",
          "to perform the ABH assignement either remove them ",
          "or use het.miss.par = TRUE")
      
    }
    
  }
  
  # test par.subcross
  
  if((!is.null(subcross.ind)) || (!is.null(par.per.subcross))){
    
    if((!is.null(subcross.ind)) && (is.null(par.per.subcross))){
      
      stop("'par.per.subcross' is not provided")
      
    }
    
    if((!is.null(par.per.subcross)) && (is.null(subcross.ind))){
      
      stop("'subcross.ind' is not provided")
      
    }
    
    
    if(!is.matrix(par.per.subcross)){
      
      stop("'par.per.subcross' is not a matrix")
      
    }
    
    if(!is.character(par.per.subcross)){
      
      stop("'par.per.subcross' is not a character matrix")
      
    }
    
    # remove the eventual rownames of par.per.subcross
    
    if(!is.null(rownames(par.per.subcross))){
      
      rownames(par.per.subcross) <- NULL
      
    }
    
    if (!identical(unique(subcross.ind), par.per.subcross[, 1])){
      
      stop("subcross identifiers used in 'subcross.ind' and 'par.per.subcross' ",
           "are not identical")
      
    }
    
    # test the similarity of parents list between par.per.subcross and
    # rownames(geno.par)
    
    parents <- union(par.per.subcross[, 2], par.per.subcross[, 3])
    
    if(sum(!(parents %in% rownames(geno.par))) > 0){
      
      list.par <- parents[!(parents %in% rownames(geno.par))]
      pbpar <- paste(list.par, collapse = ", ")
      
      message <- sprintf(ngettext(length(list.par),
                                  "parent %s is used in 'par.per.subcross' but not in 'geno.par'",
                                  "parents %s are used in 'par.per.subcross' but not in 'geno.par'"),
                         pbpar)
      
      stop(message)
      
    }
    
    
  }
  
  
  # 3. check argument for IBD computation
  #######################################
  
  
  if (!(type %in% c("F", "BC", "RIL", "DH", "BCsFt"))) {
    
    stop("'type' must be ", dQuote("F"), ', ', dQuote("BC"), ', ',
         dQuote("RIL"), ', ', dQuote("DH"), ' or ', dQuote("BCsFt"))
    
  }
  
  ### Check specification of generation number for BC F and BCsFt populations
  
  if (type == "F"){
    
    if(is.null(F.gen)){
      
      stop("'F.gen' is not specified") }
    
  }
  
  if (type == "BC"){
    
    if(is.null(BC.gen)){
      
      stop("'BC.gen' is not specified") }
    
  }
  
  if (type == "BCsFt"){
    
    if(is.null(F.gen)){ stop("'F.gen' is not specified") }
    if(is.null(BC.gen)){stop("'BC.gen' is not specified") }
    
  }
  
  
  if((type == "RIL") && is.null(type.mating)){
    
    stop("'type.mating' is missing")
    
  }
  
  # type.mating
  
  if(!is.null(type.mating)){
    
    if(!is.character(type.mating)){stop("'type.mating' must be character")}
    
    if(!(type.mating %in% c('selfing', 'sib.mat'))){
      
      stop("'type.mating' must be ", dQuote("selfing"), ' or ', dQuote("sib.mat"))
      
    }
    
  }
  
  # map.function
  
  if(!is.character(map.function)){
    
    stop("'map.function' must be character")
    
  }
  
  if(!(map.function %in% c('haldane', 'kosambi', 'c-f', 'morgan'))){
    
    stop("'map.function' must be ", dQuote("haldane"), ', ', dQuote("kosambi"), ', ',
         dQuote("c-f"), ' or ', dQuote("morgan"))
    
  }
  
}