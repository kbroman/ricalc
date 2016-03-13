\name{count.absorb}
\alias{count.absorb}

\title{Count absorbing states}

\description{
  Determine the absorbing states of the Markov chain for the generation
  of recombinant inbred lines.
}

\usage{
count.absorb(n.strains=c("2","4"),type=c("selfing","sibmating"),
         chrtype=c("A","X"),n.loci=c("2","3"))
}

\arguments{
  \item{n.strains}{Number of parental strains.}
  \item{type}{Method for mating.}
  \item{chrtype}{Indicates autosome or X chromosome (when
    \code{type=sibmating} only).}
  \item{n.loci}{Number of loci.}
}

\value{
  A vector whose names indicate the minimal set of absorbing states, and
  whose elements are the number of the full states corresponding to each.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
  count.absorb()
}

\references{
  Haldane, J. B. S. and Waddington, C. H. (1931) Inbreeding and
  linkage. \emph{Genetics} \bold{16}, 357--374.

Broman, K. W. (2005) The genomes of recombinant inbred
lines. \emph{Genetics} \bold{169}, 1133--1146.
}

\keyword{datagen}
