\name{fulltm}
\alias{fulltm}

\title{Symbolic generation-to-generation transition matrices}

\description{
  Symbolic generation-to-generation transition matrices for the
  process of creating multiple strain recombinant inbred lines.
}

\usage{data(fulltm)}

\format{
  A list of lists; each list contains the non-zero elements in the rows
  of the transition matrices.
}

\source{
  These were formed using the internal function
  \code{get.full.tm.symbolic}, with help from mathematica to simplify
  the expressions.
}

\examples{
data(fulltm)
fulltm[["2sibX2"]]["AA|BB x AA"]
}

\seealso{ \code{\link{get.full.tm}},
  \code{\link{get.full.tm.symbolic}},
  \code{\link{convert.full.tm}}
}

\keyword{datasets}
