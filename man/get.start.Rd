\name{get.start}
\alias{get.start}

\title{Determine the initial state in the generation of RILs}

\description{
  Determine the initial state in the Markov chain to generate 2- or
  4-way recombination inbred lines.
}

\usage{
get.start(n.strains=c("2","4"),type=c("selfing","sibmating"),
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
  A character string.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
get.start()
}
                    
\references{
  Haldane, J. B. S. and Waddington, C. H. (1931) Inbreeding and
  linkage. \emph{Genetics} \bold{16}, 357--374.

Broman, K. W. (2005) The genomes of recombinant inbred
lines. \emph{Genetics} \bold{169}, 1133--1146.
}

\seealso{ \code{\link{count.absorb}} }

\keyword{datagen}
