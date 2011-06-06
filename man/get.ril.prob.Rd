\name{get.ril.prob}
\alias{get.ril.prob}

\title{Calculate the haplotype probabilities for an RIL chromosome}

\description{
  Calculate (numerically) the haplotype probabilities on a 2- or 4-way
  RIL chromosome.
}

\usage{
get.ril.prob(r, coinc=1, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"), verbose=TRUE)}

\arguments{
  \item{r}{Recombination fraction between markers.}
  \item{coinc}{The 3-point coincidence (used only when
    \code{n.loci=3}).}
  \item{n.strains}{Number of parental strains.}
  \item{type}{Method for mating.}
  \item{chrtype}{Indicates autosome or X chromosome (when
    \code{type=sibmating} only).}
  \item{n.loci}{Number of loci.}
  \item{verbose}{If TRUE, print stuff as things progress.}
}

\value{
  A vector of probabilities, with absorbing states in their reduced
  form. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

                    
\references{
  Haldane, J. B. S. and Waddington, C. H. (1931) Inbreeding and
  linkage. \emph{Genetics} \bold{16}, 357--374.

Broman, K. W. (2005) The genomes of recombinant inbred
lines. \emph{Genetics} \bold{169}, 1133--1146.
}

\seealso{ \code{\link{get.ril.coinc}} }

\keyword{datagen}
