\name{par_clu}
\alias{par_clu}
\docType{data}
\title{
  Parental clustering
}
\description{
  Example of parental clustering object.
}
\usage{data(par_clu)}

\details{
  The parent clustering matrix specifies at each genome position the results
  of a parent clustering into ancestral groups. The matrix rows
  represent the position and the columns correspond to each parent. For example,
  if we have at the ith row (1, 2, 3, 2, 1), this means that parents 1
  and 5 are in the same group, that 2 and 4 are in another one and that the
  third parent was assigned to any group.
  
}

\examples{
  data(par_clu)
}

\seealso{
  \code{\link{parent_cluster.mppData}}
}

\keyword{datasets}