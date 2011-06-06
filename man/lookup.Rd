\name{lookup}
\alias{lookup}

\title{Look-up tables for multi-locus parental types}

\description{ 
  Look-up tables for multi-locus parental types in the generation of
  2-way or 4-way recombinant inbred lines by selfing, sibling mating
  with the X chromosome, or sibling mating with an autosome, and for
  cases of two or three loci.  The tables connect all possible parental
  types to the prototypical versions after collapsing according to
  various symmetries.
} 

\usage{
data(lookup)
}

\format{
  A list of vectors, each corresponding to a different set of conditions
  (for example, 2-way selfing at 3 points).  Each vector has names
  corresponding to all possible parental types and with entries being
  the corresponding prototypes after collapsing according to various
  symmetries.
}

\source{
  These were formed using the internal function \code{gtypes}.
}

\references{
Broman, K. W. (2005) The genomes of recombinant inbred
lines. \emph{Genetics} \bold{169}, 1133--1146.
}

\examples{
data(lookup)
length(lookup[["2self3"]])
length(unique(lookup[["2self3"]]))
}

%\seealso{ \code{\link{sim.ri}} }

\keyword{datasets}
