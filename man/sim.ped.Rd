\name{sim.ped}
\alias{sim.ped}

\title{Simulate a general pedigree}

\description{
  Simulates alleles for a general pedigree, using the chi-square
  model for interference at meiosis.
}

\usage{
sim.ped(pedigree, L, sexsp=1, xchr=FALSE, m=10, obligate.chiasma=FALSE)
}

\arguments{
 \item{pedigree}{A matrix with rows corresponding to individuals and
 columns being individual identifiers (\code{"id"}), sex
 (\code{0}/\code{1} or \code{F}/\code{M}), mom, and dad.  Parents must
 precede their offspring.}
 \item{L}{Length of chromosome in cM.}
 \item{sexsp}{Female:male recombination rate (must be > 0)}
 \item{xchr}{If TRUE, simulate the X chromosome.}
 \item{m}{The interference parameter (a non-negative integer).
   \code{m=0} corresponds to no interference; \code{m>0} corresponds to
   positive crossover interference.}
 \item{obligate.chiasma}{Indicates whether there is an obligate chiasma
   on the four-strand bundle.}
}

\details{
  Meiosis is simulated by the chi-square model (see Zhao et al. 1995),
  with the possibility of requiring an obligate chiasma on the
  four-strand bundle at meiosis.

  The \code{pedigree} matrix must have columns named \code{"id"},
  \code{"sex"}, \code{"mom"}/\code{"dam"} and
  \code{"dad"}/\code{"sire"}.  Other columns are ignored.
  Founders should have 0 in the mom and dad columns; others should have
  the identifiers for the mom and dad.  Each individual must have either
  both parents in the matrix or neither.

}

\value{
  The output is a list with length equal to the number of rows in
  \code{pedigree}.  Each component is itself a list with two components,
  the maternal and paternal chromosomes.  These are matrices with two
  rows:  The first row consists of the locations of exchanges
  (0 and \code{L} are included); the second row
  contains the allele (coded 1, 2, \dots, no. founders) in the interval to the
  left.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
ailped <- genAILped(ngen=5, npairs=5, sibship.size=2)
aildat <- sim.ped(ailped, 100)
}

\references{
  Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis of crossover
  interference using the chi-square model. \emph{Genetics} \bold{139},
  1045--1056.
}

\seealso{ \code{\link{sim.ri}}, \code{\link{genAILped}}, \code{\link{convert2gen}} }

\keyword{datagen}
