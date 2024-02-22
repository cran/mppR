##############
# QTL_report #
##############

# function making a report with the results of the QTL detection.

# arguments

# out.file: output location

# main: title of the QTL analysis

# QTL.info: data.frame with QTL position, log10pval and optional CI information.

# QTL.effects : list of QTL genetic effects

# R2 : QTL R2 statistics

QTL_report <- function(out.file, main, QTL.info, QTL.effects, R2){
  
  # write the title
  
  cat(main, file = out.file, append = TRUE)
  cat("\n", file = out.file, append = TRUE)
  cat(rep("*", nchar(main)), sep = "", file = out.file, append = TRUE)
  cat("\n", file = out.file, append = TRUE)
  cat("\n", file = out.file, append = TRUE)
  
  # write the number of QTL and total R squared
  
  cat(paste("Number of QTL(s):", dim(QTL.info)[1]), file = out.file, 
      append = TRUE)
  cat("\n", file = out.file, append = TRUE)
  cat("\n", file = out.file, append = TRUE)
  
  # add global R squared measurments
  
  cat(paste("Global R2:", round(R2[[1]], 2)), file = out.file, append = TRUE)
  
  cat("\n", file = out.file, append = TRUE)
  cat("\n", file = out.file, append = TRUE)
  
  cat(paste("Global adjusted R2:", round(R2[[2]], 2)), file = out.file,
      append = TRUE)
  
  cat("\n", file = out.file, append = TRUE)
  cat("\n", file = out.file, append = TRUE)
  
  
  sink(file = out.file, append = TRUE)
  
  for (i in 1:length(QTL.effects)){
    
    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    cat(paste("QTL",i), file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    cat(rep("-",5), sep = "", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    cat("QTL position:", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    
    # QTL position information (with or without CI information)
    
    print.data.frame(QTL.info[i, ], row.names = FALSE)
    
    
    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    cat("QTL effect per cross or parent:", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    
    print(QTL.effects[[i]])
    
    # partial adjusted R squared
    
    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    
    cat(paste("Partial adj. R2 (%):", round(R2[[4]][i], 2)), file = out.file,
        append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    
  }
  
  sink()
  
}