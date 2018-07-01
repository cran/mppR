############
# MQE_proc #
############

#' Multi-QTL effect MPP analysis
#' 
#' Build multi-QTL effects (MQE) models in which different QTL effects
#' (cross-specific, parental, ancestral or bi-allelic) can be assumed at
#' different loci.
#' 
#' The possible QTL effect that the user wants to allow must be
#' specified in \code{Q.eff}. The procedure is the following:
#' 
#' \enumerate{
#' 
#' \item{Forward regression to determine a MQE model with different
#' possible assumptions for the QTL effect at different loci. The function
#' use.}
#' 
#' \item{Optional backward elimination (\code{backward = TRUE}) on the final
#' list of detected QTLs.} 
#' 
#' \item{Estimation of the QTL genetic effects and R squared statistics.}
#' 
#' \item{Optional plot (\code{plot.MQE = TRUE}) of the last CIM run of the
#' forward regression using the function.}
#' 
#' }
#' 
#' 
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP_MQE".
#' 
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#' 
#' @param mppData An object of class \code{mppData}.
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#' 
#' @param Q.eff \code{Character} vector of possible QTL effects the user want to
#' test. Elements of \code{Q.eff} can be "cr", "par", "anc" or "biall".
#' For details look at \code{\link{mpp_SIM}}.
#'
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
#' @param plot.MQE \code{Logical} value. If \code{plot.MQE = TRUE},
#' the function will plot the last run of the MQE model determination.
#' Default = FALSE.
#' 
#' @param n.cores \code{Numeric}. Specify here the number of cores you like to
#' use. Default = 1.
#' 
#' @param verbose \code{Logical} value indicating if the steps of MQE_proc should
#' be printed. Default = TRUE.
#' 
#' @param output.loc Path where a folder will be created to save the results.
#' 
#' @return Return:
#'
#' \code{List} containing the following items:
#' 
#' \item{n.QTL}{Number of detected QTLs.}
#' 
#' \item{QTL}{\code{Data.frame} with QTL positions.}
#' 
#' \item{R2}{\code{list} containing R squared statistics of the QTL effects.
#' for details see \code{\link{QTL_R2}}.}
#' 
#' \item{QTL.effects}{\code{List} of genetic effects per QTL.}
#' 
#'
#' Some output files are also saved at the location specified
#' (\code{output.loc}):
#' 
#' \enumerate{
#' 
#' \item{A QTL report (QTL_REPORT.txt) with: 1) the number of detected QTLs;
#' 2) the global R squared statistics; 3) for each QTL, position information
#'  and estimated QTL genetic effect per cross or parents.}
#' 
#' \item{The list of QTLs (QTL.txt).}
#' 
#' \item{The QTL R squared statistics (QTL_R2.txt) (for details see
#'  \code{\link{QTL_R2}}).}
#' 
#' \item{General results of the QTL detection process: Number of QTL and
#' global adjusted and non-adjusted R squared statistics. (QTL_genResults.txt).}
#' 
#' \item{if \code{plot.MQE = TRUE}, a plot of the last QTL detection run profile
#' (plot_MQE.pdf).}
#' 
#' 
#' }
#'
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mpp_perm}}, \code{\link{mpp_SIM}},
#' \code{\link{QTL_R2}}
#' 
#' @examples
#' 
#' data(mppData)
#' 
#' # Specify a location where your results will be saved
#' my.loc <- tempdir()
#' 
#' MQE <- MQE_proc(pop.name = "USNAM", trait.name = "ULA", mppData = mppData,
#'                 Q.eff = c("par", "biall"), verbose = FALSE,
#'                 output.loc = my.loc)
#'                  
#' 
#' @export
#' 



MQE_proc <- function(pop.name = "MPP_MQE", trait.name = "trait1",
                     mppData = NULL, trait = 1, Q.eff,
                     threshold = 4, window = 30, backward = TRUE,
                     alpha.bk = 0.05, plot.MQE = FALSE, n.cores = 1,
                     verbose = TRUE, output.loc) {
  
  # 1. check the format of the data
  #################################
  
  check.MQE(mppData = mppData, trait = trait, Q.eff = Q.eff,
            VCOV = 'h.err', n.cores = n.cores, output.loc = output.loc,
            fct = "proc")
  
  
  # 2. Create a directory to store the results
  ############################################
  
  # create a directory to store the results of the QTL analysis
  
  folder.loc <- file.path(output.loc, paste("MQE", pop.name, trait.name,
                                            sep = "_"))
  
  dir.create(folder.loc)
  
  # 3. forward QTL detection
  ##########################
  
  if(verbose){
    
    cat("\n")
    cat("QTL detection")
    cat("\n")
    cat("\n")
    
  }
    
    QTL <- MQE_forward(mppData = mppData, trait = trait, Q.eff = Q.eff,
                       VCOV = 'h.err', threshold = threshold, window = window,
                       n.cores = n.cores, verbose = verbose)
    
    if(is.null(QTL)){
      
      warning(paste("No position is above the threshold, the stepwise procedure",
                    "could not select any QTL position."))
      
      return(NULL)
      
    } else {
      
      if(backward){
        
        if(verbose){
          
          cat("\n")
          cat("Backward elimination")
          cat("\n")
          cat("\n")
          
        }
        
        QTL <- MQE_BackElim(mppData = mppData, trait = trait, QTL = QTL[, 1],
                            Q.eff = QTL[, 5], VCOV = 'h.err', alpha = alpha.bk)
        
        if (is.null(QTL)) { # test if QTL have been selected
          
          stop("no QTL position stayed in the model after the backward elimination")
          
        }
        
      }
      
      # save the list of QTL
      
      file.QTL <- paste(folder.loc, "/", "QTL.txt", sep = "")
      
      write.table(QTL, file = file.QTL, quote = FALSE, row.names = FALSE, 
                  sep = "\t")
      
      # 4. QTL effects and R2
      #######################
      
      if(verbose){
        
        cat("\n")
        cat("Computation QTL genetic effects and R2")
        cat("\n")
        cat("\n")
        
      }  
      
      QTL_effect <- MQE_genEffects(mppData = mppData, trait = trait,
                                   QTL = QTL[, 1], Q.eff = QTL[, 5],
                                   VCOV = 'h.err')
      
      R2 <- MQE_R2(mppData = mppData, trait = trait, QTL = QTL[, 1],
                   Q.eff = QTL[, 5], glb.only = FALSE)
      
      
      # save R2 results
      
      QTL.R2 <- data.frame(QTL[, 1:5], round(R2[[3]], 2), round(R2[[4]], 2),
                           round(R2[[5]], 2), round(R2[[6]], 2),
                           stringsAsFactors = FALSE)
      
      colnames(QTL.R2)[6:9] <- c("R2.diff", "adj.R2.diff", "R2.sg", "adj.R2.sg")
      
      write.table(QTL.R2, file = paste0(folder.loc, "/", "QTL_R2.txt"),
                  quote = FALSE, sep = "\t", row.names = FALSE)
      
      # 5. Optional plot
      ##################
      
      if(plot.MQE){
        
        if(verbose){
          
          cat("\n")
          cat("Plot MQE last run profile")
          cat("\n")
          cat("\n")
          
        }
        
        CIM <- MQE_CIM(mppData = mppData, trait = trait, VCOV = 'h.err',
                       cofactors = QTL[, 1], cof.Qeff = QTL[, 5],
                       chg.Qeff = TRUE, window = window, n.cores = n.cores)
        
        main.plot <- paste("MQE", pop.name, trait.name)
        
        pdf(paste0(folder.loc, "/", "plot_MQE.pdf"), height = 10, width = 16)
        
        print(MQE_plot(mppData = mppData, Qprof = CIM, QTL = QTL, window = window,
                       threshold = threshold, main = main.plot))
        
        dev.off()
        
      }
      
      # 6. results processing
      #######################
      
      if(verbose){
        
        cat("\n")
        cat("Results processing")
        cat("\n")
        cat("\n")
        
      }
      
      
      # QTL report
      
      QTL_report(out.file = paste0(folder.loc, "/", "QTL_REPORT.txt"),
                 main = paste(pop.name, trait.name), QTL.info = QTL,
                 QTL.effects = QTL_effect, R2 = R2)
      
      # save general results
      
      gen.res <- c(dim(QTL)[1], round(R2[[1]][1], 2), round(R2[[2]][1], 2))
      names(gen.res) <- c("nb.QTL", "glb.R2", "glb.adj.R2")
      
      write.table(gen.res, file = paste0(folder.loc, "/", "QTL_genResults.txt"),
                  quote = FALSE, sep = "\t", col.names = FALSE)
      
      
      
      # form the R object to be returned
      
      results <- list(n.QTL = dim(QTL)[1], QTL = QTL, R2 = R2,
                      QTL.effects = QTL_effect)
      
      return(results)
      
    }
  
}