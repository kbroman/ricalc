\name{get.ril.coinc}
\alias{get.ril.coinc}

\title{Calculate 3-point coincidence on an RIL chromosome}

\description{
  Calculate (numerically) the quantity analogous to a three-point
  coincidence on a 2- or 4-way RIL chromosome.
}

\usage{
get.ril.coinc(r, coinc=1, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), verbose=TRUE)
}

\arguments{
  \item{r}{Recombination fraction between markers.}
  \item{coinc}{The 3-point coincidence (used only when
    \code{n.loci=3}).}
  \item{n.strains}{Number of parental strains.}
  \item{type}{Method for mating.}
  \item{chrtype}{Indicates autosome or X chromosome (when
    \code{type=sibmating} only).}
  \item{verbose}{If TRUE, print stuff as things progress.}
}

\value{
  A single number.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
  get.ril.coinc(0.01)
}
                    
\references{
  Broman, K. W. The genomes of recombinant inbred lines.
  \emph{Genetics}, in press.
}

\seealso{ \code{\link{get.ril.prob}} }

\keyword{datagen}
