\name{coinc.gam}
\alias{coinc.gam}

\title{Three-point coincidence under the gamma renewal model}

\description{
  Calculates the three-point coincidence under the gamma renewal model,
  under the assumption that the two defined intervals have the same
  recombination fraction.
}

\usage{
coinc.gam(r, nu=1, tol=1e-12)
}

\arguments{
 \item{r}{Recombination fraction for each interval (can be a vector).}
 \item{nu}{Interference parameter (\code{nu = m+1} in the chi-square
   model; \code{nu=0} is NI).}
 \item{tol}{Tolerance for doing numerical integration.}
}

\details{
  Uses the functions \code{\link{mf.gam}} and
  \code{\link{imf.gam}}.
}

\value{
  Vector of three-point coincidences, of same length as the input
  \code{r}.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
coinc.gam(seq(0.01,0.49,by=0.01), 11.3)
}

\references{
Zhao, H. and Speed, T. P. (1996) On genetic map
functions. \emph{Genetics} \bold{142}, 1369--1377.

  Broman,  K. W., Rowe,  L. B., Churchill,  G. A. and Paigen,  K. (2002)
  Crossover interference in the mouse. \emph{Genetics} \bold{160},
  1123--1131.
}

\seealso{ \code{\link{mf.gam}}, \code{\link{imf.gam}} }

\keyword{utilities}
