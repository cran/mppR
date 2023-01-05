#################
# QTL_report_GE #
#################

# Function to produce a small report with the results of a QTL MPP GxE analysis.

QTL_report_GE <- function(out.file, main, QTL.info, QTL.effects,
                          R2 = NULL){

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

  if(!is.null(R2)){

    cat(paste("Global R2:", round(R2[[1]], 2)), file = out.file, append = TRUE)

    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)

    cat(paste("Global adjusted R2:", round(R2[[2]], 2)), file = out.file,
        append = TRUE)

    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)

  }

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
    cat("QTL effects:", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)

    print(QTL.effects[[i]])

    # partial adjusted R squared

    cat("\n", file = out.file, append = TRUE)
    cat("\n", file = out.file, append = TRUE)

    if(!is.null(R2)){

      cat(paste("Partial adj. R2 (%):", round(R2[[4]][i], 2)), file = out.file,
          append = TRUE)
      cat("\n", file = out.file, append = TRUE)
      cat("\n", file = out.file, append = TRUE)

    }

  }

  sink()

}
