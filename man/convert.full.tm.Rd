\name{convert.full.tm}
\alias{convert.full.tm}

\title{Convert sparse version of transition matrix to full form}

\description{
  Convert a sparse version of a transition matrix (with only the
  non-zero elements) into a full matrix.

}

\usage{
convert.full.tm(full.tm)
}

\arguments{
  \item{full.tm}{A list whose elements are the rows of the matrix, each
    being a vector with just the non-zero elements of the matrix.}
}

\value{
  A full matrix version of the input.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }


\seealso{ \code{\link{get.full.tm}},
  \code{\link{get.full.tm.symbolic}} }

\keyword{manip}
