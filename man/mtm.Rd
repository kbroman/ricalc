\name{mtm}
\alias{mtm}


\title{Symbolic transition matrices for meiosis}

\description{
  Transition matrices for meiosis: the probabilities of the possible
  haplotypes for a meiotic product from a given parental type.  These
  are symbolic versions.
}

\usage{
data(mtm)
}

\format{
  A list of lists.  Each list corresponds to a different condition (for
  example, 2 alleles for the X chromosome at 3 points).  Each list
  contains a set of vectors, one for each possible parental type,
  containing a vector of character strings specifying the probabilities
  of each possible meiotic haplotype.
}

\source{
  These were formed using the internal functions
  \code{write.meiosis.tm.symbolic} and \code{read.meiosis.tm.symbolic},
  using mathematica to simplify the expressions.
}

\examples{
data(mtm)
sapply(mtm[["2A3"]],length)
}

%\seealso{ \code{\link{write.meiosis.tm.symbolic}} }

\keyword{datasets}
