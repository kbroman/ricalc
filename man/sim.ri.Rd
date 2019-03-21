\name{sim.ri}
\alias{sim.ri}

\title{Simulate a recombinant inbred line}

\description{
  Simulates alleles for a multiple-strain RI line, using the chi-square
  model for interference at meiosis, performing inbreeding by selfing or
  sibling mating, until complete fixation has occurred.
}

\usage{
sim.ri(L, sexsp=1, type=c("selfing","sibmating"),
       n.strains=c("2","4","8"),
       xchr=FALSE, m=10, obligate.chiasma=FALSE)
}

\arguments{
 \item{L}{Length of chromosome(s) in cM; either a vector or a single
   number.}
 \item{sexsp}{Female:male recombination rate (must be > 0)}
 \item{type}{Indicates whether the inbreeding is by selfing or by
   sibling mating.}
 \item{n.strains}{Number of initial parental strains.}
 \item{xchr}{If \code{length(L)==1} and \code{type="sibmating"}, this
     indicates to simulate the X chromosome.  If \code{length(L) > 1},
       the chromosomes named \code{"X"} or \code{"x"} are taken to be
       the X chromosome and all others are assumed to be autosomes.}
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

  Inbreeding is performed until complete fixation.  Of course, we assume
  no mutation and no selection.
}

\value{
  If \code{length(L)==1}, the output is a two-row matrix whose first
  row consists of the locations of exchanges along the recombinant
  inbred chromosome (0 and \code{L} are included).  The second row
  contains the allele (coded 1, 2, \dots, 8) in the interval to the
  left. Attributes \code{"prop.het"} and \code{"num.het"} are included
  which contain the proportion of the genome that is not yet fixed and
  the number of segments of such heterozygosity, respectively, after
  each generation of inbreeding.  The attribute \code{"nubreak"} gives
  the number of unique breakpoints at each generation.

  If \code{length(L)>1}, the output is a list of length
  \code{length(L)}, with each component being a two-row matrix, as
  above.  The overall list also has attributes \code{"prop.het"} and
  \code{"num.het"}, specifying the proportion of the entire genome that
  has not yet been fixed and the number of segments of such
  heterozygosity, respectively, after each generation of inbreeding.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(mouseL)
ri <- sim.ri(mouseL, type="selfing", n.str="8")
}

\references{
  Haldane, J. B. S. and Waddington, C. H. (1931) Inbreeding and
  linkage. \emph{Genetics} \bold{16}, 357--374.

Broman, K. W. (2005) The genomes of recombinant inbred
lines. \emph{Genetics} \bold{169}, 1133--1146.

  Broman,  K. W., Rowe,  L. B., Churchill,  G. A. and Paigen,  K. (2002)
  Crossover interference in the mouse. \emph{Genetics} \bold{160},
  1123--1131.

  Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis of crossover
  interference using the chi-square model. \emph{Genetics} \bold{139},
  1045--1056.
}

\seealso{ \code{\link{sim.ri}} }

\keyword{datagen}
