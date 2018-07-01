##########
# mpp_CV #
##########

#' MPP cross-validation
#' 
#' Evaluation of MPP QTL detection procedure by cross-validation (CV).
#' 
#' For details on the MPP QTL detection models see \code{\link{mpp_SIM}}
#' documentation. The CV scheme is adapted from Utz et al. (2000) to the MPP
#' context. A single CV run works like that:
#' 
#' \enumerate{
#' 
#' \item{Generation of a k-fold partition of the data. The partition is done
#' within crosses. Each cross is divided into k subsets. Then for the kth
#' repetition, the kth subset is used as validation set, the rest goes into the
#' training set.}
#' 
#' \item{For the kth repetition, utilization of the training set for cofactor
#' selection and multi-QTL model determination (\code{\link{mpp_SIM}} and
#' \code{\link{mpp_CIM}}). If \code{backward = TRUE}, the final list of QTLs is
#' tested simultaneously using a backward elimination
#' (\code{\link{mpp_back_elim}}).}
#' 
#' \item{Use the list of detected QTLs in the training set to calculate
#' the proportion of genetic variance explained by all detected QTLs in the
#' training set (p.ts = R2.ts/h2). Where R2.ts is the adjusted
#' R squared and h2 is the average within cross heritability (\code{her}). By
#' default, her = 1, which mean that
#' 
#' For each single QTL effect, difference partial R squared are also
#' calculated. Difference R squared are computed by doing the difference between
#' a model with all QTLs and a model without the ith position. For details about R
#' squared computation and adjustment look at \code{\link{QTL_R2}}.}
#' 
#' \item{Use the estimates of the QTL effects in the training set (B.ts) to
#' predict the phenotypic values of the validation set. y.pred.vs = X.vs*B.ts.
#' Computes the predicted R squared  in the validation set using the squared
#' Pearson correlation coefficient between the real values (y.vs) and the
#' predicted values (y.pred.vs). R2.vs = cor(y.ts,y.pred.ts)^2. Then
#' the predicted genetic variance in the validation set (p.vs) is equal to
#' p.vs = R2.vs/h2. For heritability correction, the user can provide a single
#' value for the average within cross heritability or a vector specifying each
#' within cross heritability. By default, \code{her = 1}, which means that the
#' results represent the proportion of phenotypic variance explained (predicted)
#' in the training (validation) sets.
#' 
#' The predicted R squared is computed per cross and then averaged
#' at the population level (p.ts). Both results are returned. Partial QTL
#' predicted R squared are also calculated using the difference between the
#' predicted R squared using all QTL and the predicted R squared without QTL i.
#' The bias between p.ts and p.vs is calculated as bias = 1 - (p.vs/p.ts).
#' 
#'   }
#' 
#' }
#' 
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP_CV".
#' 
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#'
#' @param mppData An object of class \code{mppData}.
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#' 
#' @param her \code{Numeric} value between 0 and 1 representing the heritability
#' of the trait. \code{her} can be a single value or a vector specifying each
#' within cross heritability. Default = 1.
#' 
#' @param Rep \code{Numeric} value representing the number of repetition of the
#' k-fold procedure. Default = 10.
#' 
#' @param k \code{Numeric} value representing the number of folds for the within
#' cross partition of the population. Default = 5.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. For more details see
#' \code{\link{mpp_SIM}}. Default = "cr".
#' 
#' @param thre.cof \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be peaked as a cofactor. Default = 3.
#' 
#' @param win.cof \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected cofactors. Default = 50.
#' 
#' @param N.cim \code{Numeric} value specifying the number of time the CIM
#' analysis is repeated. Default = 1.
#' 
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 20.
#' 
#' @param thre.QTL \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be selected as QTL. Default = 3.
#' 
#' @param win.QTL \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected QTLs. Default = 20.
#' 
#' @param backward \code{Logical} value. If \code{backward = TRUE},
#' the function performs a backward elimination on the list of selected QTLs.
#' Default = TRUE.
#' 
#' @param alpha.bk \code{Numeric} value indicating the significance level for
#' the backward elimination. Terms with p-values above this value will
#' iteratively be removed. Default = 0.05.
#' 
#' @param n.cores \code{Numeric}. Specify here the number of cores you like to
#' use. Default = 1.
#'
#' @param verbose \code{Logical} value indicating if the progresses of the CV
#' should be printed. Default = TRUE.
#'
#' @param output.loc Path where a folder will be created to save the results.
#' 
#'    
#' @return 
#' 
#' \code{List} containing the following results items:
#' 
#' \item{CV_res}{\code{Data.frame} containing for each CV run: 1) the number
#' of detected QTL; 2) the proportion of explained genetic variance in the TS
#' (p.ts); 3) the proportion of predicted genetic variance in the VS (p.vs) at
#' the population level (average of within cross prediction); the bias between
#' p.ts and p.vs (bias = 1-(p.vs/p.ts)).}
#' 
#' \item{p.vs.cr}{\code{Matrix} containing the within cross p.vs for each CV run.}
#' 
#' \item{QTL}{\code{Data.frame} containing: 1) the list of QTL position detected
#' at least one time during the entire CV process; 2) the number of times
#' the position has been detected; 3) the average partial p.ts of the QTL
#' position; 4) the average partial p.vs of the QTL position; 5) the average
#' partial bias of the QTL position.}
#' 
#' \item{QTL.profiles}{\code{Data.frame} -log10(p-value) QTL profiles of the
#' different CV runs.}
#' 
#' 
#' The results elements return as R object are also saved as text
#' files at the specified output location (\code{output.loc}). A transparency
#' plot of the CV results (plot.pdf) is also saved.
#' 
#' 
#' @author Vincent Garin
#' 
#' @references
#' 
#' Utz, H. F., Melchinger, A. E., & Schon, C. C. (2000). Bias and sampling error
#' of the estimated proportion of genotypic variance explained by quantitative
#' trait loci determined from experimental data in maize using cross validation
#' and validation with independent samples. Genetics, 154(4), 1839-1849.
#' 
#' @seealso
#' 
#' \code{\link{mpp_back_elim}},
#' \code{\link{mpp_CIM}},
#' \code{\link{mpp_perm}},
#' \code{\link{mpp_SIM}},
#' \code{\link{QTL_R2}}
#' 
#' @examples
#' 
#' data(mppData)
#' 
#' # Specify a location where your results will be saved
#' my.loc <- tempdir()
#' 
#' CV <- mpp_CV(pop.name = "USNAM", trait.name = "ULA", mppData = mppData,
#' her = .4, Rep = 1, k = 3, verbose = FALSE, output.loc = my.loc)
#' 
#' @export
#'


mpp_CV <- function(pop.name = "MPP_CV", trait.name = "trait1",
                   mppData, trait = 1, her = 1, Rep = 10, k = 5, Q.eff = "cr",
                   thre.cof = 3, win.cof = 50, N.cim = 1, window = 20,
                   thre.QTL = 3, win.QTL = 20, backward = TRUE,
                   alpha.bk = 0.05, n.cores = 1, verbose = TRUE,
                   output.loc)
{
  
  # 1. Check the validity of the parameters that have been introduced
  ###################################################################
  
  check.mpp.cv(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = 'h.err',
               n.cores = n.cores, output.loc = output.loc, her = her)
  
  # 2. Create a directory to store the results
  ############################################
  
  # create a directory to store the results of the QTL analysis
  
  folder.loc <- file.path(output.loc, paste("CV", pop.name, trait.name,
                                            Q.eff, sep = "_"))
  
  dir.create(folder.loc)
  
  # Build optional cluster
  
  if(n.cores > 1){
    
    parallel <- TRUE
    cluster <- makeCluster(n.cores)
    
  } else {
    
    parallel <- FALSE
    cluster <- NULL
    
  }
  
  # 3. Create space to store the results
  ######################################
  
  # global results
  
  N.QTL <- p.ts <- p.vs <- bias <- rep(0, (k*Rep))
  
  # within cross pVS
  
  p.vs.cr <- matrix(NA, mppData$n.cr, (k*Rep))
  
  ind.res <- 1 # index to feed the results later
  
  # individual QTL position results
  
  QTL.positions <- rep(0, dim(mppData$map)[1])
  N.Qeff.est.ts <- rep(0, dim(mppData$map)[1])
  N.Qeff.est.vs <- rep(0, dim(mppData$map)[1])
  
  QTL.pts.d <- QTL.pvs.d <- rep(0, dim(mppData$map)[1])
  
  profiles <- c()
  
  # keep the marker and in between position full list
  
  mk.list <- mppData$map[, 1]
  
  
  # 4. start to loop from 1 to r replicates
  ########################################
  
  for (i in 1:Rep) {
    
    if(verbose){
      
      cat(paste("CV repetition", i))
      cat("\n")
      
    }
    
    
    ### 4.1 generate a CV partition
    
    folds <- CV_partition(cross.ind = mppData$cross.ind, k = k)
    
    ### 4.2 iterate through the CV partition
    
    for (j in 1:k) {
      
      if(verbose){
        
        cat(paste("fold", j))
        cat("\n")
        
      }
      
      # training set
      
      mppData.ts <- subset(x = mppData, gen.list = folds[[j]]$train.set)
      
      # validation set
      
      mppData.vs <- subset(x = mppData, gen.list = folds[[j]]$val.set)
      
      prob.prog <- FALSE # indicator variable for a programmation problem
      
      # 4.2.1 cofactors selection
      
      SIM <- mpp_SIM_clu(mppData = mppData.ts, trait = trait, Q.eff = Q.eff,
                        parallel = parallel, cluster = cluster)
      
      if(sum(SIM$log10pval) == 0){prob.prog <- TRUE }
      
      cofactors <- QTL_select(Qprof = SIM, threshold = thre.cof,
                              window = win.cof, verbose = FALSE)
      
      
      # 4.2.2 multi-QTL model search
      
      # test if some cofactors have been selected
      
      if (!is.null(cofactors)) {
        
        # there are some cofactors
        
        CIM <- mpp_CIM_clu(mppData = mppData.ts, trait = trait, Q.eff = Q.eff,
                       cofactors = cofactors, window = window,
                       parallel = parallel, cluster = cluster)
        
        if(sum(CIM$log10pval) == 0){prob.prog <- TRUE }
        
        if (N.cim > 1) {
          
          for (l in 1:(N.cim - 1)) {
            
            # take the cofactors of the previous analysis
            
            cofactors <- QTL_select(Qprof = CIM, threshold = thre.cof,
                                    window = win.cof, verbose = FALSE)
            
            # test if some cofactors there before running next CIM
            
            if (!is.null(cofactors)) {
              
              CIM <- mpp_CIM_clu(mppData = mppData.ts, trait = trait,
                                 Q.eff = Q.eff, cofactors = cofactors,
                                 window = window, parallel = parallel,
                                 cluster = cluster)
              
              if(sum(CIM$log10pval) == 0){prob.prog <- TRUE }
              
              # else leave the loop
              
            } else { break }
            
          }
          
        }
        
        #### end multi QTL search
        
        # 4.2.3 QTL selection
        
        QTL <- QTL_select(Qprof = CIM, threshold = thre.QTL, window = win.QTL,
                          verbose = FALSE)
        
        # 4.2.4 Optional backward elimination
        
        if(!is.null(QTL) & backward){
          
          QTL.back <- mpp_back_elim(mppData = mppData.ts, trait = trait,
                                   QTL = QTL, Q.eff = Q.eff, alpha = alpha.bk)
          
          if(is.null(QTL.back)){ # If there was QTL position and backward return
            # no QTL it is (probably) due to programming error.
            
            prob.prog <- TRUE
            QTL <- NULL
            
          } else {QTL <- QTL.back}
          
        }
        
        
        if (!is.null(QTL)) {
          
          
          ### 4.3 compute the CV statistics (N.QTL, p.ts. etc.)
          
          # a) N.QTL
          
          N.QTL[ind.res] <- dim(QTL)[1]
          
          # b) QTL positions
          
          QTL.names <- QTL[, 1]
          
          QTL.positions <- QTL.positions + (is.element(mk.list, QTL.names) * 1)
          
          # store the profiles
          
          profiles <- cbind(profiles, CIM$log10pval)
          
          # c) p.ts (adjusted R2 / heritability)
          
          # get the R squared training set
          ################################
          
          
          R2.ts <- QTL_R2(mppData = mppData.ts, trait = trait, QTL = QTL,
                          Q.eff = Q.eff)
          
          
          # compute predicted R squared
          ##############################
          
          R2.vs <- QTL_pred_R2(mppData.ts = mppData.ts, mppData.vs = mppData.vs,
                               trait = trait, Q.eff = Q.eff, QTL = QTL,
                               her = her)
          
          
          # global results
          ################
          
          # non adjusted
          
          p.ts[ind.res] <- R2.ts$glb.R2/mean(her)
          p.vs[ind.res] <- R2.vs$glb.R2
          bias[ind.res] <- round(1 - (p.vs[ind.res]/p.ts[ind.res]), 2)
          
          
          # within cross results
          
          p.vs.cr[unique(mppData$cross.ind) %in%
                    names(R2.vs$R2.cr), ind.res] <- R2.vs$R2.cr
          
          
          # store partial R squared
          #########################
          
          # count the number of time partial QTL R2 could be estimated 
          
          if(!is.na(R2.vs[[1]])){ # validation set
            
            N.Qeff.est.vs <- N.Qeff.est.vs + (is.element(mk.list, QTL.names) * 1)
            
          }
          
          if(!is.na(R2.ts[[1]])){ # test set
            
            N.Qeff.est.ts <- N.Qeff.est.ts + (is.element(mk.list, QTL.names) * 1)
            
          }
          
          
          # general function to store individual QTL results
          
          update.res <- function(input, output, Q.names, mk.list){
            
            R.ij <- rep(0, length(mk.list))
            R.ij[match(Q.names, mk.list)] <- input
            rowSums(data.frame(output, R.ij), na.rm = TRUE)
            
          }
          
          
          # difference QTL R2 values
          
          QTL.pts.d <- update.res(input = R2.ts$part.R2.diff/mean(her),
                                  output = QTL.pts.d, Q.names = QTL.names,
                                  mk.list = mk.list)
          
          QTL.pvs.d <- update.res(input = R2.vs$part.R2.diff,
                                  output = QTL.pvs.d, Q.names = QTL.names,
                                  mk.list = mk.list)
          
        } else {
          
          
          if(prob.prog){
            
            N.QTL[ind.res] <- p.ts[ind.res] <- p.vs[ind.res] <- NA
            bias[ind.res] <- NA
            
          } else { # only no QTL significant enough.
            
            # N.QTL p.ts and p.vs stay at 0 as in the initialisation. Only need to
            # put the bias at NA.
            
            bias[ind.res] <- NA
            
            profiles <- cbind(profiles, CIM$log10pval)
            
          }
          
        }
        
      } else {
        
        
        if(prob.prog){
          
          N.QTL[ind.res] <- p.ts[ind.res] <- p.vs[ind.res] <- NA
          bias[ind.res] <- NA
          
        } else { # Only no QTL significant enough.
          
          # N.QTL p.ts and p.vs stay at 0 as in the initialisation. Only need to
          # put the bias at NA.
          
          bias[ind.res] <- NA
          
          profiles <- cbind(profiles, SIM$log10pval)
          
        }
        
      }
      
      ind.res <- ind.res + 1
      
    }  # end jth fold loop
    
    
  }  # end ith repetition loop
  
  # stop the clusters
  
  if(n.cores > 1){stopCluster(cluster)}
  
  # 5. format the results of the CV process
  #########################################
  
  ### 5.1 global results
  
  CV_res <- data.frame(N.QTL, p.ts, p.vs, bias)
  CV_res <- round(CV_res, 2)
  
  write.table(CV_res, paste0(folder.loc, "/", "CV_res.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)
  
  p.vs.cr <- round(p.vs.cr, 2)
  rownames(p.vs.cr) <- unique(mppData$cross.ind)
  colnames(p.vs.cr) <- paste0("CV.run", 1:(Rep*k))
  
  write.table(p.vs.cr, paste0(folder.loc, "/", "p_vs_cr.txt"), quote = FALSE,
              sep = "\t")
  
  ### 5.2 ind QTL position res
  
  av.pts.d <- round(QTL.pts.d/N.Qeff.est.ts, 1)
  av.pvs.d <- round(QTL.pvs.d/N.Qeff.est.vs, 1)
  bias.d <- round(1 - (av.pvs.d/av.pts.d), 1)
  
  QTL.sum <- data.frame(mppData$map, QTL.positions, av.pts.d, av.pvs.d,
                        bias.d, stringsAsFactors = FALSE)
  
  QTL.sum <- QTL.sum[QTL.positions > 0, ]
  
  colnames(QTL.sum)[5:8] <- c("N", "p.ts", "p.vs", "bias")
  
  write.table(QTL.sum, paste0(folder.loc, "/", "QTL.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)
  
  ### 5.3 profiles
  
  file <- paste0(folder.loc, "/", "QTL_profiles.txt")
  
  # combine the profiles with the map
  
  profiles2 <- data.frame(mppData$map, QTL.positions, profiles,
                          stringsAsFactors = FALSE)
  
  colnames(profiles2)[5:dim(profiles2)[2]] <- c("N.QTL", paste0("QTL.prof_",
                                                                1:dim(profiles)[2]))
  
  write.table(profiles2, file, quote = FALSE, sep = "\t", row.names = FALSE)
  
  # produce profile plot
  
  pdf(paste0(folder.loc, "/", "plot.pdf"), height = 10, width = 16)
  
  plot_CV(CV.res = list(QTL.profiles = profiles2),
          main = paste("CV", trait.name, Q.eff))
  
  dev.off()
  
  # return R object
  
  return(list(CV_res = CV_res, p.vs.cr = p.vs.cr, QTL = QTL.sum,
              QTL.profiles = profiles2))
  
  
} 
