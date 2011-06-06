\name{threept}
\alias{threept.AI}
\alias{threept.ANI}
\alias{threept.XI}
\alias{threept.XNI}

\title{Three-point probabilities for four-way RILs by sibling mating}

\description{
  Three-point haplotype probabilities for four-way recombinant inbred
  lines formed by sibling mating.  A and X indicate autosome and the X
  chromosome.  I and NI indicate strong crossover interference (gamma
  model with \eqn{\nu=11.3}{nu=11.3}) and no interference.
} 

\usage{
data(threept.AI)
data(threept.ANI)
data(threept.XI)
data(threept.XNI)
}

\format{
  Each is matrix; the first two columns indicate the recombination
  fraction between markers (at meiosis) and the three-point coincidence
  (at meiosis).  The remaining columns are the haplotype probabilities
  for the fixed RIL chromosome.  These are for the reduced set of
  states, taking account of many symmetries.

  Use \code{\link{count.absorb}} to get the numbers of full
  states in each reduced state.
}

\source{
  These were calculated using perl programs included with this package. 
}

\examples{
  data(threept.AI)
  plot(threept.AI[,c(1,3)])
}

\seealso{ \code{\link{get.ril.prob}} }

\references{
  Haldane, J. B. S. and Waddington, C. H. (1931) Inbreeding and
  linkage. \emph{Genetics} \bold{16}, 357--374.

Broman, K. W. (2005) The genomes of recombinant inbred
lines. \emph{Genetics} \bold{169}, 1133--1146.
}


\keyword{datasets}
