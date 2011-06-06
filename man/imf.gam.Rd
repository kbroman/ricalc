\name{imf.gam}
\alias{imf.gam}

\title{Inverse map function for gamma renewal model}

\description{
  Calculates the genetic length (in cM) corresponding to a particular
  recombination fraction for the gamma renewal model for recombination at
  meiosis.
}

\usage{
imf.gam(r, nu=1, tol=1e-12)
}

\arguments{
 \item{r}{Recombination fraction, potentially a vector.}
 \item{nu}{Interference parameter (\code{nu = m+1} in the chi-square
   model; \code{nu=0} is NI).}
 \item{tol}{Tolerance for doing numerical integration and finding root.}
}

\details{
  Uses \code{\link[stats]{uniroot}} to solve the \code{r = mf.gam(d,
    nu)}.  
}  

\value{
  Vector of interval lengths (in cM), of same length as the input \code{r}.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
imf.gam(seq(0,0.49,by=0.01), 11.3)
}
                    
\references{
Zhao, H. and Speed, T. P. (1996) On genetic map
functions. \emph{Genetics} \bold{142}, 1369--1377.

  Broman,  K. W., Rowe,  L. B., Churchill,  G. A. and Paigen,  K. (2002)
  Crossover interference in the mouse. \emph{Genetics} \bold{160},
  1123--1131.
}

\seealso{ \code{\link{mf.gam}}, \code{\link{coinc.gam}} }

\keyword{utilities}
