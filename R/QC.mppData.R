##############
# QC.mppData #
##############

#' Quality control for \code{mppData} objects
#' 
#' Perform different operations of quality control (QC) on the marker data of an
#' \code{mppData} object.
#' 
#' The different operations of the quality control are the following:
#' 
#' \enumerate{
#' 
#' \item{Remove markers with more than two alleles.}
#' 
#' \item{Remove markers that are monomorphic or fully missing in the parents.}
#' 
#' \item{Remove markers with a missing rate higher than \code{mk.miss}.}
#' 
#' \item{Remove genotypes with more missing markers than \code{gen.miss}.}
#' 
#' \item{Remove crosses with less than \code{n.lim} genotypes.}
#' 
#' \item{Keep only the most polymorphic marker when multiple markers map at the
#' same position.} 
#' 
#' \item{Check marker minor allele frequency (MAF). Different strategy can be
#' used to control marker MAF:
#' 
#' A) A first possibility is to filter marker based on MAF at the whole population
#' level using \code{MAF.pop.lim}, and/or on MAF within crosses using
#' \code{MAF.cr.lim}.
#' 
#' The user can give the its own vector of critical values for MAF within cross
#' using \code{MAF.cr.lim}. By default, the within cross MAF values are defined
#' by the following function of the cross-size n.ci: MAF(n.ci) = 0.5 if n.ci c
#' [0, 10] and MAF(n.ci) = (4.5/n.ci) + 0.05 if n.ci > 10. This means that up
#' to 10 genotypes, the critical within cross MAF is set to 50%. Then it
#' decreases when the number of genotype increases until 5% set as a lower bound. 
#' 
#' If the within cross MAF is below the limit in at least one cross, then marker
#' scores of the problematic cross are either put as missing
#' (\code{MAF.cr.miss = TRUE}) or the whole marker is discarded
#' (\code{MAF.cr.miss = FALSE}). By default, \code{MAF.cr.miss = TRUE} which
#' allows to include a larger number of markers and to cover a wider genetic
#' diversity.
#' 
#' B) An alternative is to select only markers that segregate in at least
#' on cross at the \code{MAF.cr.lim2} rate.
#' 
#'   }
#' 
#' }
#' 
#' @param mppData  An object of class \code{mppData} formed with
#' \code{\link{create.mppData}}.
#' 
#' @param mk.miss \code{Numeric} maximum marker missing rate at the whole
#' population level comprised between 0 and 1. Default = 0.1.
#' 
#' @param gen.miss \code{Numeric} maximum genotype missing rate at the whole
#' population level comprised between 0 and 1. Default = 0.25.
#' 
#' @param n.lim \code{Numeric} value specifying the minimum cross size.
#' Default = 15.
#' 
#' @param MAF.pop.lim \code{Numeric} minimum marker minor allele frequency at
#' the population level. Default = 0.05.
#' 
#' @param MAF.cr.lim \code{Numeric vector} specifying the critical within cross
#' MAF. Marker with a problematic segregation rate in at least
#' one cross is either set as missing within the problematic cross
#' (\code{MAF.cr.miss = TRUE}), or remove from the marker matrix
#' (\code{MAF.cr.miss = FALSE}). For default value see details.
#' 
#' @param MAF.cr.miss \code{Logical} value specifying if maker with a too low
#' segregation rate within cross (\code{MAF.cr.lim}) should be put as missing
#' (\code{MAF.cr.miss = TRUE}) or discarded (\code{MAF.cr.miss = FALSE}).
#' Default = TRUE.
#' 
#' @param MAF.cr.lim2 \code{Numeric}. Alternative option for marker MAF
#' filtering. Only markers segregating with a MAF larger than \code{MAF.cr.lim2}
#' in at least one cross will be kept for the analysis. Default = NULL.
#' 
#' @param verbose \code{Logical} value indicating if the steps of the QC should
#' be printed. Default = TRUE.
#' 
#' @param n.cores \code{Numeric}. Specify here the number of cores you like to
#' use. Default = 1.
#' 
#' @return
#' 
#' a filtered \code{mppData} object containing the the same elements
#' as \code{\link{create.mppData}} after filtering. It contains also the
#' following new elements:
#' 
#' \item{geno.id}{ \code{Character} vector of genotpes identifiers.}
#' 
#' \item{ped.mat}{Four columns \code{data.frame}: 1) the type of genotype:
#' "offspring" for the last genration and "founder" for the genotypes above
#' the offspring in the pedigree; 2) the genotype indicator; 3-4) the parent 1
#' (2) of each line.}
#' 
#' \item{geno.par.clu}{Parent marker matrix without monomorphic or completely
#' missing markers.}
#' 
#' \item{haplo.map}{Genetic map corresponding to the list of marker of the
#' \code{geno.par.clu} object.}
#' 
#' \item{parents}{List of parents.}
#' 
#' \item{n.cr}{Number of crosses.}
#' 
#' \item{n.par}{Number of parents.}
#' 
#' \item{rem.mk}{Vector of markers that have been removed.}
#' 
#' \item{rem.geno}{Vector of genotypes that have been removed.}
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{create.mppData}}
#' 
#' @examples
#' 
#' data(mppData_init)
#' 
#' mppData <- QC.mppData(mppData = mppData_init, n.lim = 15, MAF.pop.lim = 0.05,
#'                       MAF.cr.miss = TRUE, mk.miss = 0.1,
#'                       gen.miss = 0.25, verbose = TRUE)      
#' 
#' @export



QC.mppData <- function(mppData, mk.miss = 0.1, gen.miss = 0.25, n.lim = 15,
                       MAF.pop.lim = 0.05, MAF.cr.lim = NULL,
                       MAF.cr.miss = TRUE, MAF.cr.lim2 = NULL, verbose = TRUE,
                       n.cores = 1){
  
  
  # 1. check the format of the data
  #################################
  
  check_QC2(mppData = mppData, n.lim = n.lim, MAF.pop.lim = MAF.pop.lim,
            mk.miss = mk.miss, gen.miss = gen.miss, MAF.cr.lim = MAF.cr.lim,
            MAF.cr.lim2 = MAF.cr.lim2, n.cores)
  
  # 2. Restore the necessary objects from the mppData object
  ##########################################################
  
  geno.off <- mppData$geno.off
  geno.par <- mppData$geno.par
  map <- mppData$map
  cross.ind <- mppData$cross.ind
  pheno <- mppData$pheno
  par.per.cross <- mppData$par.per.cross
  
  
  # 3. Space for results and build cluster
  #########################################
  
  init.nb.mk <- dim(geno.off)[2]
  nb.mk.rem <- c()
  init.nb.gen <- dim(geno.off)[1]
  nb.gen.rem <- c()
  
  prob.mk.list <- c()
  prob.gen.list <- c()
  
  if(n.cores > 1){
    
    parallel <- TRUE
    cluster <- makeCluster(n.cores)
    
  } else {
    
    parallel <- FALSE
    cluster <- NULL
    
  }
  
  # 4. Remove markers with genotyping error
  #########################################
  
  prob.mk <- QC_GenotypingError(mk.mat = rbind(geno.par, geno.off),
                                parallel = parallel, cluster = cluster)
  
  if(is.null(prob.mk)) {rem.mk_i <- 0 } else {rem.mk_i <- length(prob.mk)}
  
  if(!is.null(prob.mk)){
    
    prob.mk.list <- c(prob.mk.list, prob.mk)
    
    ind.prob <- which(colnames(geno.off) %in% prob.mk)
    geno.par <- geno.par[, -ind.prob]
    geno.off <- geno.off[, -ind.prob]
    nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
    
  }
  
  if(verbose){
    
    cat("\n")
    cat(paste("Check genotyping error                        :",
              rem.mk_i, "markers removed", "\n"))
    
  }
  
  # 5. remove monomorphic markers in the parents
  ##############################################
  
  parent.MAF <- QC_MAF(mk.mat = geno.par, parallel = parallel,
                       cluster = cluster)
  
  mono <- which(parent.MAF == 0)
  miss <- which(is.na(parent.MAF))
  
  prob.mk.id <- c(mono, miss)
  
  if(is.null(prob.mk.id)) {rem.mk_i <- 0 } else {rem.mk_i <- length(prob.mk.id)}
  
  if(length(prob.mk.id) > 0){
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
    
    geno.par <- geno.par[, -prob.mk.id]
    geno.off <- geno.off[, -prob.mk.id]
    nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
    
  }
  
  if(verbose){
    
    cat(paste("Remove monomorphic/missing marker in parents  :",
              rem.mk_i, "markers removed", "\n"))
    
  }
  
  ### keep parent genotype and corresponding map to be used for
  # clustering. later
  
  geno.par.clu <- geno.par
  
  map.par.clu <- QC_matchMarker(mk.mat = geno.par.clu, map = map)[[2]]
  
  
  # 6. Remove markers with too high missing rate at the population level
  ######################################################################
  
  miss.ind.mk <- QC_missing(mk.mat = geno.off, threshold = mk.miss)
  
  if(dim(miss.ind.mk)[1] > 0) {rem.mk_i <- dim(miss.ind.mk)[1]
  } else {rem.mk_i <- 0}
  
  if(dim(miss.ind.mk)[1] > 0){
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[miss.ind.mk[, 2]])
    
    geno.par <- geno.par[, -miss.ind.mk[, 2]]
    geno.off <- geno.off[, -miss.ind.mk[, 2]]
    nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
    
  }
  
  if(verbose){
    
    cat(paste("Remove marker with missing rate >", mk.miss, "        :",
              rem.mk_i, "markers removed", "\n"))
    
  }
  
  # 7. Remove genotypes with too high missing rate at the population level
  ########################################################################
  
  geno.ref <- rownames(geno.off) # make a reference list of genotypes
  
  miss.ind.gen <- QC_missing(mk.mat = geno.off, MARGIN = 1,
                             threshold = gen.miss)
  
  if(dim(miss.ind.gen)[1] > 0) {rem.gen_i <- dim(miss.ind.gen)[1]
  } else {rem.gen_i <- 0}
  
  if(dim(miss.ind.gen)[1] > 0){
    
    geno.off <- geno.off[-miss.ind.gen[, 2], ]
    nb.gen.rem <- c(nb.gen.rem, rem.gen_i)
    prob.gen.list <- c(prob.gen.list, as.character(miss.ind.gen[, 1]))
    
    # adapt the other arguments which depend on the genotype list
    
    ind.geno <- geno.ref %in% rownames(geno.off)
    cross.ind <- cross.ind[ind.geno]
    pheno <- pheno[ind.geno, ]
    
  }
  
  if(verbose){
    
    cat(paste("Remove genotype with missing rate >", gen.miss, "     :",
              rem.gen_i, "genotypes removed", "\n"))
    
  }
  
  
  # 8. Remove cross with a too small size
  #######################################
  
  geno.ref <- rownames(geno.off) # make a reference list of genotypes
  
  geno.off <- QC_minCrSize(mk.mat = geno.off, cross.ind = cross.ind,
                           n.lim = n.lim)
  
  rem.gen_i <- length(geno.ref) - dim(geno.off)[1]
  
  if(rem.gen_i > 0){
    
    nb.gen.rem <- c(nb.gen.rem, rem.gen_i)
    gen.removed <- geno.ref[!(geno.ref %in% rownames(geno.off))]
    prob.gen.list <- c(prob.gen.list, gen.removed)
    
    # adapt the other arguments which depend on the genotype list
    
    ind.geno <- geno.ref %in% rownames(geno.off)  
    cross.ind <- cross.ind[ind.geno]
    pheno <- pheno[ind.geno, ]
    
  }
  
  if(verbose){
    
    cat(paste("Remove crosses with less than", n.lim, "observations :",
              rem.gen_i, "genotypes removed", "\n"))
    
  }
  
  
  # 9. Remove less polymorphic marker(s) if some markers are at the same position
  ###############################################################################
  
  # select maximum 1 marker per position
  
  map <- QC_matchMarker(mk.mat = geno.off, map = map)[[2]]
  
  difference <- diff(map[, 3])
  difference <- c(1, difference) # add 1 for the first position.
  
  rem.mk_j <- sum(difference == 0)
  
  if(rem.mk_j > 0){
    
    MAF.pop <- QC_MAF(mk.mat = geno.off, parallel = parallel, 
                      cluster = cluster) 
    
    # Identify blocks of marker that are at the same position.
    
    list.pos <- list()
    i <- 1
    max <- length(difference)
    list.pos.id <- 1
    
    while(i <= max){
      
      # test if the value is zero -> same position
      
      if(difference[i] == 0){
        
        vec <- c(i-1, i) # start the vector
        i <- i + 1
        
        if(i > max){
          list.pos[[list.pos.id]] <- vec
          break()
        }
        
        while(difference[i] == 0){ # continue as long as the values are zero
          
          vec <- c(vec, i)
          i <- i + 1
          
          if(i > max){
            list.pos[[list.pos.id]] <- vec
            break()
          }
          
        } # see if there are other positions
        
        list.pos[[list.pos.id]] <- vec # store the result
        list.pos.id <- list.pos.id + 1
        
      } else {
        
        i <- i + 1
        
      }
      
    }
    
    # identify the most polymorphic marker within the pre-selected blocks and
    # therefore the positions to remove.
    
    ref.MAF <- MAF.pop[colnames(geno.off)]
    
    prob.mk.id <- c()
    
    for(i in 1:length(list.pos)){
      
      rem.mk_i <- list.pos[[i]][-which.max(ref.MAF[list.pos[[i]]])]
      
      prob.mk.id <- c(prob.mk.id, rem.mk_i)
      
    }
    
    prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
    
    geno.par <- geno.par[, -prob.mk.id]
    geno.off <- geno.off[, -prob.mk.id]
    map <- map[-prob.mk.id, ]
    
    nb.mk.rem <- c(nb.mk.rem, rem.mk_j)
    
  }
  
  if(verbose){
    
    cat(paste("Remove markers at the same position           :",
              rem.mk_j, "markers removed", "\n"))
    
  }
  
  # 10. Marker MAF sorting
  ########################
  
  ### 10.1 : Select markers segregating in at least one cross
  
  if(!is.null(MAF.cr.lim2)){
    
    
    off.MAF <- QC_MAF(mk.mat = geno.off, cross.ind = cross.ind,
                      parallel = parallel, cluster = cluster)
    
    MAF.cr <- off.MAF[[2]]
    
    lim <- rep(MAF.cr.lim2, length(unique(cross.ind)))
    mk.sel <- rep(TRUE, dim(geno.off)[2])
    
    for(i in 1:dim(geno.off)[2]){
      
      test <- MAF.cr[, i] > lim
      test[is.na(test)] <- FALSE
      
      if(sum(test) > 0){ mk.sel[i] <- TRUE } else {mk.sel[i] <- FALSE }
      
    }
    
    prob.mk.id <- which(!mk.sel)
    
    if(is.null(prob.mk.id)) {rem.mk_i <- 0 } else {rem.mk_i <- length(prob.mk.id)}
    
    
    if(length(prob.mk.id) > 0){
      
      prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
      
      geno.par <- geno.par[, -prob.mk.id]
      geno.off <- geno.off[, -prob.mk.id]
      nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
      
    }
    
    ### message
    
    if(verbose){
      
      cat(paste("Remove markers with MAF <", MAF.cr.lim2, "in at least one cross :",
                rem.mk_i, "markers removed", "\n"))
      
    }
    
    
  } else {
    
    ### 10.2 : Select markers segregating whole population and within crosses
    
    ### 10.2.1: Marker at the whole population level
    
    off.MAF <- QC_MAF(mk.mat = geno.off, cross.ind = cross.ind,
                      parallel = parallel, cluster = cluster)
    
    
    MAF.pop <- off.MAF[[1]]
    MAF.cr <- off.MAF[[2]]
    
    prob.mk.id <- which(MAF.pop < MAF.pop.lim)
    
    if(is.null(prob.mk.id)) {rem.mk_i <- 0 } else {rem.mk_i <- length(prob.mk.id)}
    
    
    if(length(prob.mk.id) > 0){
      
      prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
      
      geno.par <- geno.par[, -prob.mk.id]
      geno.off <- geno.off[, -prob.mk.id]
      MAF.cr <- MAF.cr[, -prob.mk.id]
      MAF.pop <- MAF.pop[-prob.mk.id]
      nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
      
    }
    
    if(verbose){
      
      cat(paste("Remove markers with MAF <", MAF.pop.lim,"               :",
                rem.mk_i, "markers removed", "\n"))
      
    }
    
    ### 10.2.2: Marker at the within cross level
    
    
    MAF.pop <-  list(MAF.pop = MAF.pop, MAF.cr = MAF.cr)
    class(MAF.pop) <- c("list", "mafRes")
    
    # functions to determine the MAF limit within crosses
    
    if(is.null(MAF.cr.lim)){
      
      MAF.lim <- function(x, floor){
        
        if(x <= 10){ 0.5 } else { (4.5/x) + floor }
        
      }
      
      # determine the number of observation per cross. First transform into
      # factor with specified order
      
      n.cr <- table(factor(cross.ind, levels = unique(cross.ind)))
      
      lim <- unlist(lapply(X = n.cr, FUN = MAF.lim, floor = 0.05))
      
    } else {
      
      lim <- MAF.cr.lim
      
    }
    
    MAF.cr.ind <- QC_tagMAFCr(MAF = MAF.pop, MAF.lim = lim, tag.mono = FALSE,
                              parallel = parallel, cluster = cluster)
    
    # two options to manage the marker with problementic MAF within cross.
    # 1.: Put these markers as missing within the cross; 2.: remove the marker 
    
    prob.mk.id <- which(MAF.cr.ind)
    rem.mk_i <- length(prob.mk.id)
    
    if(MAF.cr.miss){ # put NA markers with prob. within cross MAF in at least 1 cross
      
      cr.id <- unique(cross.ind)
      
      for(i in 1:dim(geno.off)[2]){
        
        test <- (MAF.cr[, i] < lim) & (MAF.cr[, i] != 0)
        test[is.na(test)] <- FALSE
        
        if(sum(test) > 0){ # at least one cross has a problematic MAF
          
          geno.off[cross.ind %in% cr.id[test], i] <- NA
          
        }
        
      }
      
      
    } else { # remove markers with problematic within cross MAF in at least 1 cross
      
      
      if(rem.mk_i > 0){
        
        prob.mk.list <- c(prob.mk.list, colnames(geno.par)[prob.mk.id])
        
        geno.par <- geno.par[, -prob.mk.id]
        geno.off <- geno.off[, -prob.mk.id]
        nb.mk.rem <- c(nb.mk.rem, rem.mk_i)
        
      }
      
      if(verbose){
        
        cat(paste("Remove markers critical MAF witin cross            :",
                  rem.mk_i, "markers removed", "\n"))
        
      }
      
    }
    
  }
  
  
  
  
  # 11. equalize the list of markers and the list of genotypes
  ############################################################
  
  
  ### 11.1 markers
  
  match <- QC_matchMarker(mk.mat = geno.off, map = map)
  map <- match$new.map
  
  
  ### 11.2 genotypes (equalize also cross.ind and subcross.ind)
    
    pheno.aug <- data.frame(pheno, cross.ind, stringsAsFactors = FALSE)
    
    inter.geno.pheno <- intersect(rownames(geno.off), rownames(pheno.aug))
    
    pheno.aug <- pheno.aug[inter.geno.pheno, ]
    cross.ind <- pheno.aug[, dim(pheno.aug)[2]]
    pheno <- pheno.aug[, -dim(pheno.aug)[2], drop = FALSE]
    pheno <- as.matrix(pheno)
    colnames(pheno) <- colnames(mppData$pheno)
   
  
  ### 11.3 Warning message if user want to use less than 15 individual per
    # crosses
  
  freq <- table(cross.ind)
  
  if(sum(freq < 15)){
    
    warning(paste("It is still possible to perform a MPP QTL analysis with",
                  "crosses containing less that 15 genotypes. However, we",
                  "advice to use minimum 15 individuals per cross to have",
                  "enough information to estimate within crosses QTL effects."))
    
  }
  
  
  ### 11.4 modify the par.per.cross argument and geno.par
  
  cr.list <- unique(cross.ind)
  par.per.cross <- par.per.cross[par.per.cross[, 1] %in% cr.list, ]
  
  par.list <- union(par.per.cross[, 2], par.per.cross[, 3])
  geno.par <- geno.par[rownames(geno.par) %in% par.list, ]
  
  # Modify also the geno.par.clu by removing unused parents
  
  geno.par.clu <- geno.par.clu[rownames(geno.par.clu) %in% par.list,]
  
  # stop the eventual clusters
  
  if(n.cores > 1){stopCluster(cluster)}
  
  # 12. Form the pedigree information matrix (ped.mat)
  ####################################################
  
  geno.id <- rownames(geno.off)
  
  p1 <- par.per.cross[, 2]
  p2 <- par.per.cross[, 3]
  
  names(p2) <- names(p1) <- par.per.cross[, 1]
  
  ped.mat <- data.frame(rep("offspring", length(geno.id)), geno.id,
                        p1[cross.ind], p2[cross.ind], stringsAsFactors = FALSE)
  colnames(ped.mat) <- c("type" ,"genotypes", "parent1", "parent2")
  
  parents <- union(par.per.cross[, 2], par.per.cross[, 3])
  
  n.par <- length(parents)
  
  n.cr <- dim(par.per.cross)[1]
  
  if(n.cr < 2){
    
    stop('The QC reduced your MPP to less than 2 crosses. You must modify your QC parameters to have at least 2 crosses in your MPP.')
    
  }
  
  # 13. re-form the mppData object with the new elements
  ######################################################
  
  mppData <- list(geno.off = geno.off, geno.IBS = NULL, geno.IBD = NULL,
                  geno.id = geno.id, ped.mat = ped.mat, allele.ref = NULL,
                  geno.par = geno.par, geno.par.clu = geno.par.clu, par.clu = NULL,
                  pheno = pheno, map = map, haplo.map = map.par.clu,
                  cross.ind = cross.ind, par.per.cross = par.per.cross,
                  type = NULL, parents = parents, n.cr = n.cr, n.par = n.par,
                  n.zigo = NULL, rem.mk = prob.mk.list, rem.gen = prob.gen.list,
                  status = 'QC')
  
  class(mppData) <- c("mppData", "list")
  
  ##### final message
  
  if(verbose){
    
    tot.rem.mk <- init.nb.mk - sum(nb.mk.rem)
    tot.rem.gen <- init.nb.gen - sum(nb.gen.rem)
    
    cat("\n")
    cat("   End             :", tot.rem.mk ,
        "marker(s) remain after the check\n")
    cat("                    ",
        tot.rem.gen , "genotypes(s) remain after the check\n")
    
    
  }
  
  return(mppData)
  
}