###############
# mpp_forward #
###############

#' MPP QTL analysis using forward regression
#' 
#' Multi-parent population QTL analysis model using a forward regression.
#' 
#' The function run a full MPP QTL detection using models with different possible
#' assumptions concerning the number of alleles at the QTL position. For more details about
#' the different models, see documentation of the function \code{\link{mpp_SIM}}.
#' The procedure is the following:
#' 
#' \enumerate{
#' 
#' \item{Forward regression to determine the a multi-QTL model. The function
#' selects successively QTL positions with -log10(pval) above the threshold.
#' Those positions are added as cofactors for following detection run.
#' The procedure stop when no more position has a -log10(pval) above the
#' threshold (\code{\link{QTL_forward}}).}
#' 
#' \item{If \code{backward = TRUE}, backward elimination  on the final
#' list of detected QTLs.} 
#' 
#' \item{Estimation of the QTL genetic effects and R squared statistics
#' (\code{\link{QTL_gen_effects}} and \code{\link{QTL_R2}}).}
#' 
#' \item{If \code{CI = TRUE}, confidence interval calculation  
#' based on a CIM- (composite interval mapping removing all cofactors on the
#' scanned chromosome) of the last run of the forward regression.}
#' 
#' \item{If \code{plot.Qprof = TRUE}, plot  of the last run of the forward
#' regression using \code{\link{plot.QTLprof}}.}
#' 
#' \item{If \code{plot.gen.eff = TRUE}, plot of the genetic effect distribution
#' along the genome of the last run of the forward regression using
#' \code{\link{plot.QTLprof}}.}
#' 
#' }
#' 
#'
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP".
#' 
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#' 
#' @param mppData An object of class \code{mppData}.
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects; 4) "biall" for a bi-allelic
#' effects. For more details see \code{\link{mpp_SIM}}. Default = "cr".
#' 
#' @param ref.par Optional \code{Character} expression defining the parental
#' allele that will be used as reference to calculate the allelic effects of 
#' the parental model. For the ancestral model, the ancestral class containing
#' the reference parent will be set as reference. \strong{This option can only
#' be used if the MPP design is composed of a unique connected part}.
#' Default = NULL.
#' 
#' @param sum_zero Optional \code{Logical} value specifying if the QTL effect of
#' a parental or an ancestral model should be caculated using the sum to zero
#' constraint. Default = FALSE.
#' 
#' @param threshold \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be considered as significant.
#' Default = 4.
#' 
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 30.
#' 
#' @param backward \code{Logical} value. If \code{backward = TRUE},
#' the function performs
#' a backward elimination on the list of selected QTLs. Default = TRUE.
#' 
#' @param alpha.bk \code{Numeric} value indicating the significance level for
#' the backward elimination. Default = 0.05.
#' 
#' @param plot.Qprof \code{Logical} value. If \code{plot = TRUE},
#' the function will plot the QTL profile last run of the forward regression.
#' Default = FALSE.
#' 
#' @param plot.gen.eff \code{Logical} value. If \code{plot.gen.eff = TRUE},
#' the function will save the decomposed genetic effects per cross/parent.
#' These results can be ploted with the function \code{\link{plot.QTLprof}}
#' to visualize a genome-wide decomposition of the genetic effects. This
#' plot will be realized on the last run of the forward regression.
#' \strong{This functionality is ony available for the cross-specific,
#' parental and ancestral models.}
#' Default value = FALSE.
#' 
#' @param CI \code{Logical} value. If \code{CI = TRUE}, the function will
#' compute a -log10(pval) drop confidence interval for each QTL using
#' the QTL profile of the last iteration. Default = FALSE.
#' 
#' @param drop \code{Numeric} -log10(p-value) drop value at the limits of the
#' interval. Default = 1.5.
#' 
#' @param text.size \code{Numeric} value specifying the size of graph axis text
#' elements. Default = 18.
#' 
#' @param n.cores \code{Numeric}. Specify here the number of cores you like to
#' use. Default = 1.
#' 
#' @param verbose \code{Logical} value indicating if the steps of the procedure
#' should be printed. Default = TRUE.
#' 
#' @param output.loc Path where a folder will be created to save the results.
#' By default the function uses the current working directory.
#' 
#' 
#' @return Return:
#' 
#' List containing the following items:
#' 
#' \item{n.QTL}{Number of detected QTLs.}
#' 
#' \item{QTL}{\code{Data.frame} with QTL positions.}
#' 
#' \item{R2}{\code{List} containing R squared statistics of the QTL effects.
#' For details see \code{\link{QTL_R2}} output section.}
#' 
#' \item{QTL.effects}{\code{List} of QTLs genetic effects. For details see
#' \code{\link{QTL_gen_effects}} output section.}
#' 
#' \item{QTL.CI}{If \code{CI = TRUE}, confidence interval information of
#' the QTLs.}
#' 
#' Some output files are also saved at the specified location
#' (\code{output.loc}):
#' 
#' \enumerate{
#' 
#' \item{A QTL report (QTL_REPORT.txt) with: 1) the number of detected QTLs;
#' 2) the global R squared statistics; 3) for each QTL, position information
#' (plus confidence interval if \code{CI = TRUE}) and estimated QTL genetic
#' effects per cross or parents (for details see \code{\link{QTL_gen_effects}}).}
#' 
#' \item{The list of QTL (QTL.txt).}
#' 
#' \item{The QTL R squared statistics (QTL_R2.txt) (for details see
#' \code{\link{QTL_R2}}).}
#' 
#' \item{If \code{CI = TRUE}, the QTL confidence intervals (QTL_CI.txt).}
#' 
#' \item{General results of the QTL detection process: number of QTLs and
#' global adjusted and non-adjusted R squared statistics (QTL_genResults.txt).}
#' 
#' \item{If \code{plot.Qprof = TRUE}, the plot of the last regression run
#' (QTL_profile.pdf). If \code{plot.gen.eff = TRUE}, plot of the genetic
#' effects per cross or parents (gen_eff.pdf) with dashed lines representing
#' the QTL positions. For more details see \code{\link{plot.QTLprof}}}
#' 
#' }
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mpp_SIM}}, \code{\link{plot.QTLprof}},
#' \code{\link{QTL_gen_effects}}, \code{\link{QTL_forward}}, \code{\link{QTL_R2}}
#'
#' @examples
#' 
#' \dontrun{
#' 
#' data(mppData)
#' 
#' # Specify a location where your results will be saved
#' my.loc <- "C:/.../..."
#' 
#' # Cross-specific model
#' 
#' USNAM_cr <- mpp_forward(pop.name = "USNAM", trait.name = "ULA",
#'                         mppData = mppData, plot.gen.eff = TRUE,
#'                         plot.Qprof = TRUE, CI = TRUE, output.loc = my.loc)
#' 
#' 
#' }
#'
#' @export
#'

mpp_forward <- function(pop.name = "MPP", trait.name = "trait1", mppData,
                        trait = 1, Q.eff = "cr", ref.par = NULL, sum_zero = FALSE,
                        threshold = 4, window = 30,
                        backward = TRUE, alpha.bk = 0.05, plot.Qprof = FALSE,
                        plot.gen.eff = FALSE, CI = FALSE, drop = 1.5,
                        text.size = 18, n.cores = 1, verbose = TRUE, output.loc) {
  
  
  # 1. Check the validity of the parameters that have been introduced
  ###################################################################
  
  check.mpp.proc(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = 'h.err',
                 plot.gen.eff = plot.gen.eff, ref.par = ref.par,
                 sum_zero = sum_zero, n.cores = n.cores,
                 output.loc = output.loc)
  
  
  # 2. Create a directory to store the results
  ############################################
  
  # create a directory to store the results of the QTL analysis
  
  folder.loc <- file.path(output.loc, paste("ForwardQTL", pop.name, trait.name,
                                            Q.eff, sep = "_"))
  
  dir.create(folder.loc)
  
  # 3. Forward regression
  #######################
  
  if(verbose) {
    
    cat("\n")
    cat("Forward regression")
    cat("\n")
    cat("\n")
    
  }
  
  QTL <- QTL_forward(mppData = mppData, trait = trait, Q.eff = Q.eff,
                     threshold = threshold, window = window,
                     n.cores = n.cores, verbose = verbose)
  
  
  # 5. Backward elimination
  #########################
  
  if (backward){
    
    if(verbose){
      
      cat("\n")
      cat("Backward elimination")
      cat("\n")
      cat("\n")
      
    }
    
    QTL <- mpp_back_elim(mppData = mppData, trait = trait, QTL = QTL,
                         Q.eff = Q.eff, alpha = alpha.bk)
    
    if (is.null(QTL)) { # test if QTL have been selected
      
      stop("no QTL position stayed in the model after the backward elimination. ",
           "This is probably due to an error in the computation of the model ",
           "in asreml function")
      
    }
    
  }
  
  # save the final list of QTLs
  
  write.table(QTL, file = file.path(folder.loc, "QTL.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  
  # 6. R squared computation
  ##########################
  
  if(verbose){
    
    cat("\n")
    cat("R squared computation")
    cat("\n")
    cat("\n")
    
  }
  
  
  R2 <- QTL_R2(mppData = mppData, trait = trait, QTL = QTL, Q.eff = Q.eff)
  
  # save R2 results
  
  QTL.R2 <- data.frame(QTL[, 1:5], round(R2[[3]], 2), round(R2[[4]], 2),
                       round(R2[[5]], 2), round(R2[[6]], 2),
                       stringsAsFactors = FALSE)
  
  colnames(QTL.R2)[6:9] <- c("R2.diff", "adj.R2.diff", "R2.sg", "adj.R2.sg")
  
  write.table(QTL.R2, file = file.path(folder.loc, "QTL_R2.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # 7. QTL effects estimation
  ###########################
  
  if(verbose){
    
    cat("\n")
    cat("QTL effects estimation")
    cat("\n")
    cat("\n")
    
  }
  
  QTL.effects <- QTL_gen_effects(mppData = mppData, trait = trait, QTL = QTL,
                                 Q.eff = Q.eff, ref.par = ref.par,
                                 sum_zero = sum_zero)
  
  # 8. CIM- and confidence interval computation
  #############################################
  
  if(CI | plot.Qprof | plot.gen.eff){
    
    # Build optional cluster
    if(n.cores > 1){
      
      parallel <- TRUE
      cluster <- makeCluster(n.cores)
      
    } else {
      
      parallel <- FALSE
      cluster <- NULL
      
    }
    
    if(CI){
      
      if(verbose){
        
        cat("\n")
        cat("CIM- and confidence intervals computation")
        cat("\n")
        cat("\n")
        
      }
      
      # determine larger chromosome distance
      map <- mppData$map
      chr.fact <- factor(x = map[, 2], levels = unique(map[, 2]))
      step.size <- max(tapply(X = map[, 4], INDEX = chr.fact, FUN = max)) + 100
      
      CIM.m <- mpp_CIM_clu(mppData = mppData, trait = trait, Q.eff = Q.eff,
                           cofactors = QTL[, 1], window = step.size,
                           plot.gen.eff = FALSE, parallel = parallel,
                           cluster = cluster)
      
      QTL.CI <- QTL_CI(QTL = QTL, Qprof = CIM.m, drop = drop)
      
      write.table(QTL.CI, file = file.path(folder.loc, "QTL_CI.txt"),
                  quote = FALSE, sep = "\t", row.names = FALSE)
      
    }
    
    
    if(plot.Qprof | plot.gen.eff) {
      
      CIM_p <- mpp_CIM_clu(mppData = mppData, trait = trait, Q.eff = Q.eff,
                           cofactors = QTL[, 1], window = window,
                           plot.gen.eff = plot.gen.eff, parallel = parallel,
                           cluster = cluster)
      
      class(CIM_p) <- c("QTLprof", "data.frame")
      
      if(plot.Qprof){
        
        main.cim <- paste("QTL profile", pop.name, trait.name, Q.eff)
        
        if (Q.eff == "biall") {
          
          pdf(file.path(folder.loc, "QTL_profile.pdf"), height = 10, width = 16)
          
          pl <- plot(x = CIM_p, type = "h", main = main.cim,
                     threshold = threshold, text.size = text.size)
          print(pl)
          
          dev.off()
          
        } else {
          
          # CIM profile
          
          pdf(file.path(folder.loc, "QTL_profile.pdf"), height = 10, width = 16)
          
          pl <- plot(x = CIM_p, QTL = QTL, type = "l", main = main.cim,
                             threshold = threshold, text.size = text.size)
          print(pl)
          
          dev.off()
          
        }
        
      }
      
      if(plot.gen.eff){
        
        main.Qeff <- paste("QTL gen. effects", pop.name, trait.name, Q.eff)
        
        pdf(file.path(folder.loc, "gen_eff.pdf"), height = 10, width = 16)
        
        pl <- plot(x = CIM_p, gen.eff = TRUE, mppData = mppData,  Q.eff = Q.eff,
                           QTL = QTL, main = main.Qeff, text.size = text.size)
        
        print(pl)
        
        dev.off()
        
      }
      
      
    }
    
    # stop the clusters
    
    if(n.cores > 1){stopCluster(cluster)}
    
  } else { QTL.CI <- NULL}
  
  
  
  # 9. Results processing
  #######################
  
  if(verbose){
    
    cat("\n")
    cat("Results processing")
    cat("\n")
    cat("\n")
    
  }
  
  ### 9.1: general results
  
  gen.res <- c(dim(QTL)[1], round(R2[[1]][1], 2), round(R2[[2]][1], 2))
  names(gen.res) <- c("nb.QTL", "glb.R2", "glb.adj.R2")
  
  write.table(gen.res, file = file.path(folder.loc, "QTL_genResults.txt"),
              quote = FALSE, sep = "\t", col.names = FALSE)
  
  
  ### 9.3: Report
  
  
  if(CI) {QTL.info <- data.frame(QTL[, c(1, 2, 4, 5)], QTL.CI[, 4:8],
                                 stringsAsFactors = FALSE)
  } else {QTL.info <-  QTL[, c(1, 2, 4, 5)]}
  
  QTL_report(out.file = file.path(folder.loc, "QTL_REPORT.txt"),
             main = paste(pop.name, trait.name, Q.eff), QTL.info = QTL.info,
             QTL.effects = QTL.effects[[1]], R2 = R2)
  
  
  ### 9.4: Return R object
  
  
  results <- list(n.QTL = dim(QTL)[1], QTL = QTL[, 1:5], R2 = R2,
                  QTL.effects = QTL.effects, QTL.CI = QTL.CI)
  
  return(results)
  
}
