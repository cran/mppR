\name{USNAM_geno}
\alias{USNAM_geno}
\docType{data}
\title{
Reduced genotype data maize US-NAM population
}
\description{
Selection of markers and genotypes from the maize US nested association mapping
(NAM) population (McMullen et al., 2009).
}
\usage{data(USNAM_geno)}
\format{
 \code{data.frame}
}
\details{
Sample of the marker matrix of the US-NAM population. The selection correspond
to 102 markers coming from the two first
chromosomes present in \code{\link{USNAM_map}} and the 506 genotypes. These
genotypes correspond to the selected phenotypic values in
\code{\link{USNAM_pheno}}. The selected genotypes come from
the following crosses: (B73 x CML103), (B73 x CML322), (B73 x CML52),
(B73 x Hp301), (B73 x M37W). The data of the 6 parental lines are
also included. The data are available on www.panzea.org.
}
\source{
\url{https://www.panzea.org/}
}
\references{
McMullen, M. D., Kresovich, S., Villeda, H. S., Bradbury, P., Li, H., Sun, Q.,
... & Buckler, E. S. (2009). Genetic properties of the maize nested association
mapping population. Science, 325(5941), 737-740.


}
\examples{
data(USNAM_geno)
}

\seealso{
\code{\link{USNAM_pheno}}, \code{\link{USNAM_map}}
}

\keyword{datasets}
