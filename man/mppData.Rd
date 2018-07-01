\name{mppData}
\alias{mppData}
\docType{data}
\title{
   Complete \code{mppData} object
}
\description{
  
  Complete \code{mppData} object made from a sample data of the maize US nested
  association mapping (NAM) population. This mppData object went through all
  the steps of the data processing: \code{\link{create.mppData}},
  \code{\link{QC.mppData}}, \code{\link{IBS.mppData}}, \code{\link{IBD.mppData}}, \code{\link{parent_cluster.mppData}}. The mppData contain all the data
  necessary for the QTL analysis procedures.
  
}
\usage{data(mppData)}
\format{
  \code{mppData}
}
\details{
  
  The complete \code{mppData} object is a \code{list} containing the following
  elements:

\enumerate{

\item{geno.IBS: IBS genotype marker matrix.}

\item{geno.IBD: R/qtl cross.object containing the genotype within cross IBD probabilities.}

\item{geno.id: List of genotypes.}

\item{allele.ref: Matrix containing for each marker the most and least frequent marker scores and the two heterozygous scores.}

\item{geno.par: Parents marker matrix.}

\item{geno.par.clu: Parent marker data used to cluster the parents.}

\item{par.clu: Parent clustering results.}

\item{mono.anc: Positions for which the ancestral clustering was monomorphic.}

\item{pheno: Phenotypic trait matrix.}

\item{map: Genetic map corresponding to the geno (IBS, IBD, par) arguments.}

\item{haplo.map: Genetic map corresponding to geno.par.clu.}

\item{cross.ind: Vector indicating to which cross each genotype belongs.}

\item{par.per.cross: Matrix with for each cross the parent 1 and 2.}

\item{parents: Vector of parents.}

\item{type: Type of population.}

\item{n.cr: Number of crosses.}

\item{n.par: Number of parents.}

\item{n.anc: Average number of ancestral group along the genome.}

\item{n.zigo: Number of possible allelic computations 2 (AA/BB) or 3 (AA/AB/BB).}

\item{rem.mk: Removed markers in the data processing.}

\item{rem.gen: Removed genotypes in the data processing.}

\item{status: Indicates the level of progression in the data processing.}

}
  
}

\examples{
  data(mppData)
}

\seealso{
  \code{\link{create.mppData}}, \code{\link{QC.mppData}},
  \code{\link{IBS.mppData}}, \code{\link{IBD.mppData}},
  \code{\link{parent_cluster.mppData}}
}

\keyword{datasets}