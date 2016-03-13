\name{mf.gam}
\alias{mf.gam}

\title{Map function for gamma renewal model}

\description{
  Calculates the recombination fraction corresponding to a particular
  cM interval length for the gamma renewal model for recombination at
  meiosis.
}

\usage{
mf.gam(d, nu=1, tol=1e-12)
}

\arguments{
 \item{d}{Interval length (in cM), potentially a vector.}
 \item{nu}{Interference parameter (\code{nu = m+1} in the chi-square
   model; \code{nu=0} is NI).}
 \item{tol}{Tolerance for doing numerical integration.}
}

\details{
  Uses the function \code{\link[stats]{integrate}} to do numerical
  integration.
}

\value{
  Vector of recombination fractions, of same length as the input \code{d}.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
mf.gam(0:25, 11.3)
}

\references{
Zhao, H. and Speed, T. P. (1996) On genetic map
functions. \emph{Genetics} \bold{142}, 1369--1377.

  Broman,  K. W., Rowe,  L. B., Churchill,  G. A. and Paigen,  K. (2002)
  Crossover interference in the mouse. \emph{Genetics} \bold{160},
  1123--1131.
}

\seealso{ \code{\link{imf.gam}}, \code{\link{coinc.gam}} }

\keyword{utilities}
