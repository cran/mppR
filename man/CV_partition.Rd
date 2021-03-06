% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CV_partition.R
\name{CV_partition}
\alias{CV_partition}
\title{Cross validation partition}
\usage{
CV_partition(cross.ind, k = 5)
}
\arguments{
\item{cross.ind}{\code{Character} vector with the same length as the number
of genotypes which specifies to which cross each genotype belongs.}

\item{k}{\code{Numeric} value representing the number of subsets (fold) into
which data are spread within cross. Default = 5.}
}
\value{
Return:

\item{fold}{\code{List} of k lists (one for each fold). Each fold list
contains two vectors with genotypes indices of the training (\code{$train.set}) and
the validation set (\code{$val.set}).}
}
\description{
Partition the genotype indices into training and validation sets for
cross-validation (CV).
}
\details{
The genotype indices are randomly assigned within cross to k subsets (folds).
Then each subset is used once as validation set, the remaining data go in the
training set.
}
\examples{

data(mppData)

part.cv <- CV_partition(cross.ind = mppData$cross.ind, k = 5)

part.cv[[1]]$train.set
part.cv[[1]]$val.set

}
\seealso{
\code{\link{mpp_CV}}
}
\author{
Vincent Garin
}
