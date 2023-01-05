\name{mppData_GE}
\alias{mppData_GE}
\docType{data}
\title{
  Example \code{mppData} object
}
\description{

Example \code{mppData} object representing a subset from the maize EU-NAM Flint
population (Bauer et al. 2013, Lehermeier et al. 2014, Giraud et al. 2014).


}
\usage{data(mppData_GE)}
\format{
  \code{mppData}
}
\details{

Sample data from the maize EU-NAM Flint population. The genotype data were
obtained from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50558.
The genotypes come from the five following crosses: (UH007 x D152),
(UH007 x F03802), (UH007 x F2), (UH007 x F283), and (UH007 x DK105). We selected
100 markers randomly spread on chromosomes five and six. The genetic map was
downloaded here http://maizegdb.org/data_center/reference?id=9024747.

The phenotypic data represent the within environment adjusted means for dry
matter yield (DMY) calculated at La Coruna (CIAM), at Roggenstein (TUM), at
Einbeck (KWS), and at Ploudaniel (INRA_P). The raw plot data were obtained
here: http://www.genetics.org/content/198/1/3/suppl/DC1.


}

\references{

Bauer, E., Falque, M., Walter, H., Bauland, C., Camisan, C., Campo, L., ...
  & Altmann, T. (2013). Intraspecific variation of recombination rate in maize.
  Genome biology, 14(9), R103.

Giraud, H., Lehermeier, C., Bauer, E., Falque, M., Segura, V., Bauland, C., ... & Schipprack, W. (2014). Linkage disequilibrium with linkage analysis of multiline
crosses reveals different multiallelic QTL for hybrid performance in the flint
and dent heterotic groups of maize. Genetics, 198(4), 1717-1734.

Lehermeier, C., Krämer, N., Bauer, E., Bauland, C., Camisan, C., Campo, L.,
  ... & Moreau, L. (2014). Usefulness of multiparental populations of maize
(Zea mays L.) for genome-based prediction. Genetics, 198(1), 3-16.

}

\examples{
  data(mppData_GE)
}

\keyword{datasets}
