\name{get.full.tm}
\alias{get.full.tm}

\title{Calculate the transition matrix for the formation of RILs}

\description{
  Calculate the full generation-to-generation transition matrix for the
  formation of RILs.
}

\usage{
get.full.tm(r, coinc=1, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"))
}

\arguments{
  \item{r}{Recombination fraction between markers.}
  \item{coinc}{The 3-point coincidence (used only when
    \code{n.loci=3}).}
  \item{n.strains}{Number of parental strains.}
  \item{type}{Method for mating.}
  \item{chrtype}{Indicates autosome or X chromosome (when
    \code{type=sibmating} only).}
  \item{n.loci}{Number of loci.}
}

\value{
  The transition matrix as a list whose components are the rows of the
  matrix.  Only the non-zero elements of each row are returned.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
  fulltm2self2 <- get.full.tm(0.1)
}

\references{
  Haldane, J. B. S. and Waddington, C. H. (1931) Inbreeding and
  linkage. \emph{Genetics} \bold{16}, 357--374.

Broman, K. W. (2005) The genomes of recombinant inbred
lines. \emph{Genetics} \bold{169}, 1133--1146.
}

\seealso{ \code{\link{get.full.tm.symbolic}},
  \code{\link{fulltm}}, \code{\link{convert.full.tm}}
}

\keyword{datagen}
