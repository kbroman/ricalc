\name{mouseL}
\alias{mouseL}

\title{Estimated genetic lengths of mouse chromosomes}

\description{
  Estimated genetic lengths (in cM) of mouse chromosomes.
}

\usage{data(mouseL)}

\format{
  An vector of length 20, each element being the length (in cM) of the
  corresponding chromosome.
}

\source{
  Mouse Genome Informatics,
  \url{http://www.informatics.jax.org/mgihome/other/mouse_facts1.shtml}
}

\examples{
data(mouseL)
ri <- sim.ri(mouseL, type="selfing", n.str="8")
}

\seealso{ \code{\link{sim.ri}} }

\keyword{datasets}
