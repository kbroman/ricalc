\name{genAILped}
\alias{genAILped}

\title{Generate a pedigree matrix for AIL}

\description{
  Generate a pedigree matrix, for use with the simulation function
  \code{\link{sim.ped}}, for advanced intercross lines.
}

\usage{
genAILped(ngen=8, npairs=20, sibship.size=10)
}

\arguments{
 \item{ngen}{Number of generations of crosses.}
 \item{npairs}{Number of mating pairs per generation.}
 \item{sibship.size}{Size of sibships in the final generation.}
}

\details{
  We generate equal numbers of males and females at each generation.
  Matings are with random male/female pairs, with no attempt to avoid
  matings between siblings.
}  

\value{
  The output is a matrix with five columns: \code{"id"}, \code{"sex"}
  (coded as \code{0}=female, \code{1}=male), \code{"mom"}, \code{"dad"},
  and \code{"generation"}.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
ailped <- genAILped(ngen=5, npairs=5, sibship.size=2)
}
                    
\seealso{ \code{\link{sim.ped}} }

\keyword{datagen}
