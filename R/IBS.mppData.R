###############
# IBS.mppData #
###############

#' IBS coding for \code{mppData} objects
#' 
#' Transform the genotype marker matrix of a \code{mppData} object into
#' Identical by state (IBS) 0, 1, 2 format. The IBS score represent the number
#' of copies of the minor allele.
# The 0, 1, 2 coding can be preceded by a
# marker imputation (\code{impute = TRUE}) to fill the missing values. The
# imputation of the missing values is performed using the \code{codeGeno()}
# function from synbreed (Wimmer et al., 2012).
#' 
#' 
#' @param mppData  An object of class \code{mppData}. The \code{mppData} must
#' have been processed using: \code{\link{create.mppData}} and
#' \code{\link{QC.mppData}}.
#' 
# @param impute \code{Logical} value. if \code{impute = TRUE}, the function
# will impute missing values using the \code{codeGeno()} function from the
# synbreed package. Default = FALSE.
# 
# @param impute.type \code{character} with one out of \code{"fix"},
# \code{"random"}, \code{"family"}, \code{"beagle"}, \code{"beagleAfterFamily"},
# \code{"beagleAfterFamilyNoRand"}. For details see synbreed package
# documentation. \strong{To be able to use Beagle for imputation, Please load
# the synbreed package using \code{library(synbreed)}} Default = "random".
# 
# @param map_bp \code{data.frame} with three columns specifying for each marker
# position the marker identifier, the \code{numeric} or \code{character}
# chromosome the and physical bp position. This argument is necessary for
# imputation using Beagle. Default = NULL.
# 
# @param replace.value \code{numeric} scalar to replace missing value in case
# \code{impute.type = fix}. Only 0, 1, 2. Should be chosen. Default = NULL.
# 
# @param label.heter This is either a scalar or vector of characters to identify
# heterozygous genotypes or a function returning TRUE if an element of the
# marker matrix is the heterozygous genotype. Defining a function is useful,
# if number of unique heterozygous genotypes is large, i.e. if genotypes are
# coded by alleles. If the heterozygous genotype is coded like
# "A/T","G/C", ..., "AG", "CG", ..., "T:C", "G:A", ... or "G|T", "A|C", ...
# then label.heter="alleleCoding" can be used. Note that
# heterozygous values must be identified unambiguously by label.heter. Use
# label.heter=NULL if there are only homozygous genotypes, i.e. in DH lines,
# to speed up computation and restrict imputation to values 0 and 2.
# Default = "alleleCoding".
#' 
#' @return
#' 
#' an increased \code{mppData} object containing the the same elements
#' as the \code{mppData} object provided as argument and the
#' following new elements:
#' 
#' \item{geno.IBS}{Marker \code{matrix} with marker scores coded as 0, 1, 2
#' corresponding to the number of copies of the least frequent SNP allele.}
#' 
#' \item{allele.ref}{\code{matrix} with reference allele scores. The first row
#' represents the minor allele (lowest frequency), the second the one represent
#' the major allele (largest frequency) and the two others the heterozygous
#' scores.}
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{create.mppData}}, \code{\link{QC.mppData}}
#' 
# @references
# 
# Wimmer, V., Albrecht, T., Auinger, H. J., & Schon, C. C. (2012). synbreed: a
# framework for the analysis of genomic prediction data using R.
# Bioinformatics, 28(15), 2086-2087.
# 
# Browning, B. L., & Browning, S. R. (2013). Improving the accuracy and
# efficiency of identity-by-descent detection in population data. Genetics,
# 194(2), 459-471.
#' 
#' @examples
#' 
#' data(mppData_init)
#' 
#' mppData <- QC.mppData(mppData_init)
#' 
#' mppData <- IBS.mppData(mppData = mppData)
#'       
#' 
#' @export

# IBS.mppData <- function(mppData, impute = FALSE, impute.type = "random",
#                         map_bp = NULL, replace.value = NULL,
#                         label.heter = "alleleCoding"){

IBS.mppData <- function(mppData){
  
  # 1. check the format of the data
  #################################
  
  # check_IBS(mppData = mppData, impute = impute, impute.type = impute.type,
  #           map_bp = map_bp, replace.value = replace.value)
  
  check_IBS(mppData = mppData)
  
  
  # 2. Restore the necessary objects from the mppData object
  ##########################################################
  
  geno.off <- mppData$geno.off
  map <- mppData$map
  cross.ind <- mppData$cross.ind
  
  geno_names <- rownames(geno.off)
  mk_names <- colnames(geno.off)
  
  # 3. Space for results and build cluster
  #########################################
  
  # if(n.cores > 1){
  #   
  #   parallel <- TRUE
  #   cluster <- makeCluster(n.cores)
  #   
  # } else {
  #   
  #   parallel <- FALSE
  #   cluster <- NULL
  #   
  # }
  
  # 4. Convert to 012 format
  ##########################
  
  geno.trans <- geno_012(mk.mat = geno.off, parallel = FALSE,
                         cluster = NULL)
  
  geno.IBS <- geno.trans[[1]]
  
  allele.ref <- geno.trans[[2]]
  
  rm(geno.trans)
  
  # if(n.cores > 1){
  #   
  #   stopCluster(cluster)
  #   rm(cluster)
  #   
  # }
  
  
  # 4. Imputation
  ###############
  
  # if(impute){
  #   
  #   # separate imputation for fixed value
  #   
  #   if(impute.type == "fix"){
  #     
  #     geno.IBS[is.na(geno.IBS)] <- replace.value
  #     
  #   } else {
  #     
  #     family <- data.frame(cross.ind, stringsAsFactors = FALSE)
  #     rownames(family) <- rownames(geno.off)
  #     
  #     if(impute.type %in% c("random", "family", "fix")) {
  #       
  #       map_gp <- map[, 2:3]
  #       map.unit <- "cM"
  #       
  #     } else { # Cases with Beagle
  #       
  #       message("the imputation using Beagle can take several minutes")
  #       
  #       map_gp <- map_bp[map_bp[, 1] %in% map[, 1], 2:3]
  #       map.unit <- "bp"
  #       
  #     }
  #     
  #     colnames(map_gp) <- c("chr", "pos")
  #     rownames(map_gp) <- map[, 1]
  #     
  #     # Uniformize the heterozygous values in a single type: AT and TA -> AT
  #     
  #     geno_temp <- geno.off
  #     het_ref <- allele.ref[3:4, ]
  #     
  #     unif_het <- function(x, ref){
  #       
  #       if(length(unique(x[!is.na(x)])) > 3){
  #         
  #         x[x == ref[1]] <- ref[2]
  #         
  #         return(x)
  #         
  #       } else { return(x) }
  #       
  #     }
  #     
  #     for(i in 1:dim(geno_temp)[2]){
  #       
  #       geno_temp[, i] <- unif_het(x = geno_temp[, i], ref = het_ref[, i])
  #       
  #     }
  #     
  #     gp <- synbreed::create.gpData(geno = geno_temp, map = map_gp, family = family,
  #                         map.unit = map.unit)
  #     
  #     gp.imp <- tryCatch(expr = synbreed::codeGeno(gpData = gp, impute = TRUE,
  #                                        impute.type = impute.type,
  #                                        replace.value = replace.value, maf = 0,
  #                                        nmiss = 1, label.heter = label.heter,
  #                                        verbose = FALSE),
  #                        error = function(e) NULL)
  #     
  #     
  #     if(is.null(gp.imp)){
  #       
  #       if(map.unit == "cM"){ warning("the missing values imputation failed")
  #         
  #       } else{ warning("the missing values imputation using Beagle failed") }
  #       
  #     } else {
  #       
  #       geno.IBS <- gp.imp$geno
  #       
  #       # put the marker and genotype in the same order
  #       
  #       geno.IBS <- geno.IBS[geno_names, ]
  #       geno.IBS <- geno.IBS[, mk_names]
  #       
  #       # Need to check that the reference allele is correct
  #       
  #       ref_all <- cbind(gp.imp$map$alter, gp.imp$map$refer) 
  #       
  #       allele.ref[1:2, ] <- t(ref_all)
  #       
  #     }
  #     
  #   }
  #   
  # }
  
  
  # 13. fill the mppData object
  #############################
  
  mppData$geno.IBS <- geno.IBS
  
  mppData$allele.ref <- allele.ref
  
  mppData$status <- 'IBS' 
  
  class(mppData) <- c("mppData", "list")
  
  return(mppData)
  
}