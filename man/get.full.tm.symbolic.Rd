\name{get.full.tm.symbolic}
\alias{get.full.tm.symbolic}

\title{Calculate symbolic form of transition matrix for the formation of RILs}

\description{
  Calculate the full generation-to-generation transition matrix for the
  formation of RILs, in symbolic form.
}

\usage{
get.full.tm.symbolic(n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"))
}

\arguments{
  \item{n.strains}{Number of parental strains.}
  \item{type}{Method for mating.}
  \item{chrtype}{Indicates autosome or X chromosome (when
    \code{type=sibmating} only).}
  \item{n.loci}{Number of loci.}
}

\value{
  The transition matrix as a list whose components are the rows of the
  matrix.  Only the non-zero elements of each row are returned, and they
  are returned as character strings.
}


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
  fulltm2self2 <- get.full.tm.symbolic()
}

\references{
  Haldane, J. B. S. and Waddington, C. H. (1931) Inbreeding and
  linkage. \emph{Genetics} \bold{16}, 357--374.

  Broman,  K. W., Rowe,  L. B., Churchill,  G. A. and Paigen,  K. (2002)
  Crossover interference in the mouse. \emph{Genetics} \bold{160},
  1123--1131.
}

\seealso{ \code{\link{get.full.tm}}, \code{\link{convert.full.tm}} }

\keyword{datagen}
