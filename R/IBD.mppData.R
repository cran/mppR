###############
# IBD.mppData #
###############

#' IBD coding for \code{mppData} objects
#' 
#' The function first converts genotype data into ABH format. Then it calculates
#' within cross identical by descent (IBD) probabilities.
#' 
#' The function first transforms genotype data into within cross ABH format.
#' The function takes the parents of the different cross as
#' reference and assigns the following scores: "A" if the offspring score is
#' equivalent to parent 1; "B" if it is equivalent to parent 2; "H" if it is
#' heterozygous. The function attributes NA for missing when: 1) the offspring
#' score is missing; 2) the two parents have the same score; or
#' 3) when at least one parental score is missing.
#' 
#' If a parent score is heterozygous or missing (\code{het.miss.par = TRUE}),
#' the assignment rules are the following. If the two parents are
#' heterozygous or one parent is heterozygous and the other missing, the
#' offspring get NA since the parental origin can not be inferred with certainty.
#' If one parent is heterozygous or missing and the second parent is
#' homozygous, the function looks at offspring segregating pattern to
#' infer which allele was transmitted by the heterozygous parent. If this is
#' possible we consider the heteroxzygous parent as homozygous for the
#' transmitted allele and use this allele score for ABH assignment.
#' 
#' The ABH assignment can be performed using sub-cross structure providing
#' information about sub-cross in arguments \code{subcross.ind} and
#' \code{par.per.subcross}. 
#' 
#' Then the function calculates the IBD probabilities using \code{read.cross()}
#' and \code{calc.genoprob()} functions from the R/qtl package
#' (Broman et al. 2009).
#' 
#' The type of population must be specified in argument \cite{type}. Different
#' population types are possible: F-type ('F'), back-cross ('BC'), backcross
#' followed by selfing ('BCsFt'), double haploid ('DH'), and recombinant imbred
#' lines ('RIL'). The number of F and BC generations can be specified using
#' \code{F.gen} and \code{BC.gen}. The argument \code{type.mating} specifies if
#' F and RIL populations were obtained by selfing or by sibling mating.
#' 
#' DH and RIL populations are read as back-cross by R/qtl. For these two
#' population types, heterozygous scores will be treated as missing values.
#' 
#' 
#' @param mppData  An object of class \code{mppData}. the \code{mppData} must
#' have been processed using: \code{\link{create.mppData}},
#' \code{\link{QC.mppData}}, and \code{\link{IBS.mppData}}.
#' 
#' @param het.miss.par \code{Logical} value. if \code{het.miss.par = TRUE},
#' the function will use the offspring segregation to try to infer the allele
#' that was transmitted by the heterozygous or missing parent at a particular
#' locus to make the ABH conversion. Default = TRUE.
#' 
#' @param subcross.ind Optional \code{character} vector specifying to which
#' sub-cross each genotype belong. Default = NULL.
#' 
#' @param par.per.subcross Optional three columns \code{Character matrix}
#' specifying : 1) the sub-cross indicators; 2) the parents 1 identifiers
#' of the sub-crosses; 3) the parents 2 identifiers of the sub-crosses.
#' Default = NULL.
#' 
#' @param type \code{Character} indicator for the type of population analyzed:
#' type = "F" for Fn (F cross n generations); type = "BC" for BCn (backcross
#' n generations); type = "BCsFt" for backcross followed by selfing;
#' type = "DH" for double haploids; and type = "RIL"
#' for recombinant inbred lines. For RIL type specify if the population was
#' obtain using selfing or sibling mating using \code{type.mating}.
#' If type = "RIL" or "DH", the function does not assume heterozygous marker
#' scores for these populations and convert them into missing (NA).
#' 
#' @param F.gen \code{Numeric} integer representing the number of F generations.
#' For example F.gen = 2 for F2. Default = NULL.
#' 
#' @param BC.gen \code{Numeric} integer representing the number of
#' backcross generations. For example BC.gen = 1 for single backcross.
#' Default = NULL.
#' 
#' @param type.mating \code{Character} specifying for a RIL population if it was
#' obtained by selfing ("selfing") or by sibling mating ("sib.mat").
#' Default = NULL.
#' 
#' @param error.prob \code{Numeric} value for assumed genotyping error rate
#' used in the calculation of the penetrance Pr(observed genotype | true genotype).
#' Default = 0.0001.
#' 
#' @param map.function \code{Character} expression specifying the type of map
#' function used to infer the IBD probabilities. possibility to choose
#' between "haldane", "kosambi","c-f","morgan". Default = "haldane".
#' 
#' 
#' @return
#' 
#' an increased \code{mppData} object containing the the same elements
#' as the \code{mppData} object provided as argument and the
#' following new elements:
#' 
#' \item{geno.IBD}{A R/qtl \code{cross.object} containing the IBD probabilities.}
#' 
#' \item{n.zigo}{\code{Numeric} value Indicating the number of different
#' genotypes: 2 (AA/BB) or 3 (AA/AB/BB)}
#' 
#' \item{type}{\code{Character} expression indicating the type of population.}
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{create.mppData}}, \code{\link{QC.mppData}},
#' \code{\link{IBS.mppData}}
#' 
#' @references
#' 
#' Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping
#' in experimental crosses. Bioinformatics 19:889-890.
#' 
#' Broman, K. W., & Sen, S. (2009). A Guide to QTL Mapping with R/qtl (Vol. 46).
#' New York: Springer.
#' 
#' @examples
#' 
#' data(mppData_init)
#' 
#' mppData <- QC.mppData(mppData_init)
#' mppData <- IBS.mppData(mppData = mppData)
#' 
#' mppData <- IBD.mppData(mppData = mppData, het.miss.par = TRUE, type = 'RIL',
#'                        type.mating = 'selfing')
#' 
#' 
#' @export



IBD.mppData <- function(mppData, het.miss.par = TRUE, subcross.ind = NULL,
                        par.per.subcross = NULL, type, F.gen = NULL,
                        BC.gen = NULL, type.mating = NULL,
                        error.prob = 1e-04, map.function = "haldane"){
  
  # 1. check the format of the data
  #################################
  
  check_IBD(mppData = mppData, het.miss.par = het.miss.par,
            subcross.ind = subcross.ind, par.per.subcross = par.per.subcross,
            type = type, F.gen = F.gen, BC.gen = BC.gen,
            type.mating = type.mating, map.function = map.function)
  
  
  # 2. Restore the necessary objects from the mppData object
  ##########################################################
  
  geno.off <- mppData$geno.off
  geno.par <- mppData$geno.par
  pheno <- mppData$pheno
  map <- mppData$map
  cross.ind <- mppData$cross.ind
  par.per.cross <- mppData$par.per.cross
  parents <- mppData$parents
  
  # 3. General elements
  #####################
  
  ### type of population
  
  if (type == "F") {
    
    type.pop <- paste0("F", "(n = ", F.gen,")")
    
  } else if (type == "BC") {
    
    type.pop <- paste0("Back-cross ", "(n = ", BC.gen,")")
    
  } else if (type == "DH") {
    
    type.pop <- "Double haploid"
    
  } else if (type == "RIL") {
    
    type.pop <- "Recombinant inbred line"
    
    if(type.mating == "selfing") {
      
      type.pop <- paste (type.pop, "by selfing")
      
    } else {
      
      type.pop <- paste (type.pop, "by sibling mating")
      
    }
    
  } else if (type == 'BCsFt'){
    
    type.pop <- paste0('Back-cross followed by selfing ', '(',
                       paste0('BC', BC.gen, 'F', F.gen), ')')
    
  }
  
  ### number of allele class
  
  if ((type == "BC") | (type == "DH") | (type == "RIL")) {
    
    n.zigo <- 2
    
  } else if ((type == "F")| (type == "BCsFt")) {
    
    n.zigo <- 3
    
  }
  
  # 4. Convert data into ABH format
  #################################
  
  
  # check if a cross was completely removed. Then remove it from the
  # par.per.cross argument.
  
  if(!is.null(subcross.ind)){
    
    cr.list <- unique(subcross.ind)
    par.per.cross.ABH <- par.per.subcross[par.per.subcross[, 1] %in% cr.list, ]
    cross.ind.ABH <- subcross.ind
    
  } else {
    
    par.per.cross.ABH <- par.per.cross 
    cross.ind.ABH <- cross.ind
    
  }
  
  ### 11.1 case with heterozygous markers scores.
  
  if(het.miss.par){
    
    # ABH assignement with heterogeneous parents
    
    geno.off <- cross_ABH_het(par.sc = geno.par, off.sc = geno.off,
                              cross.ind = cross.ind.ABH,
                              par.per.cross = par.per.cross.ABH)
    
    
  } else {
    
    geno.off <- cross_ABH(par.sc = geno.par, off.sc = geno.off,
                          cross.ind = cross.ind.ABH,
                          par.per.cross = par.per.cross.ABH)
    
  }
  
  
  # 5. Form the cross object 
  ##########################
  
  # convert heterozygous into missing for RIL and DH population
  
  if(type %in% c('RIL', 'DH')){
    
    geno.off[geno.off == 'H'] <- NA
    
  }
  
  # format data to form a cross object
  
  chr.info <- t(map[, 2:3])
  colnames(chr.info) <- map[, 1]
  
  geno.aug <- rbind(chr.info, geno.off)
  
  n_pheno <- dim(pheno)[2]
  
  empty_mat <- matrix("", nrow = 2, ncol = n_pheno)
  
  trait.aug <- rbind(empty_mat, pheno)
  geno.aug <- cbind(trait.aug, geno.aug)
  
  # Export the data in a .csv file in a temporary file.
  
  tmp <- tempfile(fileext = "Cross_object.csv")
  
  write.csv(geno.aug, file = tmp, row.names = FALSE)
  
  # form a R/qtl cross object reading the data using the specifiec type of
  # population
  
  if (type == "F") {
    
    cross.object <- read.cross("csv", , tmp, F.gen = F.gen, 
                               crosstype = "bcsft")
    
    
  } else if (type == "BC") {
    
    cross.object <- read.cross("csv", , tmp, BC.gen = BC.gen, 
                               crosstype = "bcsft")
    
    
  } else if (type == "RIL") {
    
    # need to read the object as a backcross
    
    cross.object <- read.cross("csv", file = tmp, genotypes = c("A", "B"),
                               alleles = c("A", "B"))
    
    # then convert it following the type of mating
    
    if (type.mating == "selfing") {
      
      cross.object <- convert2riself(cross.object)
      
    }
    
    if (type.mating == "sib.mat") {
      
      cross.object <- convert2risib(cross.object)
      
    }
    
  } else if (type == "DH") {
    
    # need to read the object as a backcross
    
    cross.object <- read.cross("csv", file = tmp, genotypes = c("A", "B"),
                               alleles = c("A", "B"))
    
    class(cross.object)[1] <- "dh"
    
  } else if (type == "BCsFt"){
    
    cross.object <- read.cross("csv", , tmp, F.gen = F.gen,
                               BC.gen = BC.gen, crosstype = "bcsft")
    
  }
  
  # 6. Compute the IBD probabilities
  ##################################
  
  # make sure no need to provide stepwidth
  
  cross.object <- calc.genoprob(cross.object, step = 0,
                                error.prob = error.prob,
                                map.function = map.function)
  
  # delete the temporary directory
  
  unlink(tmp)
  
  # 7. transform the map and geno.par
  ###################################
  
  c(rep(1, 3), rep(2, 5),rep(10, 4))
  
  chr.ind <- factor(x = map[, 2], levels = unique(map[, 2]))
  
  new.map <- data.frame(map[, 1:2], sequence(table(chr.ind)),
                        map[, 3], stringsAsFactors = FALSE)
  
  colnames(new.map) <- c("mk.names", "chr", "pos.ind", "pos.cM")
  
  ### geno.par
  
  geno.par.new <- t(geno.par)
  
  geno.par.new <- data.frame(new.map, geno.par.new, stringsAsFactors = FALSE)
  colnames(geno.par.new)[5:dim(geno.par.new)[2]] <- parents
  
  # 8. fill the mppData object
  #############################
  
  mppData$geno.IBD <- cross.object
  
  mppData$type <- type.pop
  
  mppData$n.zigo <- n.zigo
  
  mppData$map <- new.map
  
  mppData$geno.par <- geno.par.new
  
  mppData$status <- 'IBD' 
  
  class(mppData) <- c("mppData", "list")
  
  return(mppData)
  
}