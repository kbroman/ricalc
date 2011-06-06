\name{ricalc-package}
\docType{package}
\alias{ricalc-package}

\title{Overview of R/ricalc package}
\description{
  A brief introduction to the contents of the R/ricalc package.
}

\section{Stuff for the gamma model}
{
  \itemize{
    \item \code{\link{coinc.gam}}:  Calculate the three-point
    coincidence function for the gamma model.

    \item \code{\link{mf.gam}}: Map function for the gamma model.

    \item \code{\link{imf.gam}}: Inverse of map function for the gamma model.
  }
}

\section{Data}
{
  \itemize{
     \item \code{\link{mouseL}}: Estimated genetic lengths of mouse chromosomes.

     \item \code{\link{threept.AI}}: Three-point haplotype probabilities
     for 8-way RIL, autosome, under interference.

     \item \code{\link{threept.ANI}}: Three-point haplotype probabilities
     for 8-way RIL, autosome, with no interference.

     \item \code{\link{threept.XI}}: Three-point haplotype probabilities
     for 8-way RIL, X chromosome, under interference.

     \item \code{\link{threept.XNI}}: Three-point haplotype probabilities
     for 8-way RIL, X chromosome, with no interference.
  }
}

\section{Symbolic data}
{
  \itemize{
     \item \code{\link{fulltm}}: Generation-to-generation transition
     matrices for 2- and 4-way RIL by selfing and by sibling mating (in
     symbolic notation...lists of lists of character strings).

     \item \code{\link{mtm}}: Transition matrices for meiosis, again in
     symbols...lists of lists of character strings).
  }

}

\section{Counting and enumerating parental types}
{
  \itemize{
     \item \code{\link{lookup}}: Lookup tables for parental types in the
     generation of multiple-strain RIL, representing the 

     \item \code{\link{get.start}}: Symbolic notation for the starting
     state in the generation of multiple-strain RIL.

     \item \code{\link{count.absorb}}: Count the absorbing states in
     multiple-strain RIL.
  }
}   

\section{Calculations and simulations}
{
  \itemize{
     \item \code{\link{sim.ri}}: Simulate multiple-strain RIL.

     \item \code{\link{get.ril.coinc}}: Numerical calculation of the
     3-point coincidence-type quantity for 2- and 4-way RIL.

     \item \code{\link{get.ril.prob}}: Numerical calculation of
     haplotype probabilities on 2- or 4-way RIL chromosomes.

     \item \code{\link{get.full.tm}}: Calculate the
     generation-to-generation transition matrix for the formation of
     RIL. 

     \item \code{\link{get.full.tm.symbolic}}: Calculate the
     generation-to-generation transition matrix for the formation of
     RIL, in symbolic form.

     \item \code{\link{convert.full.tm}}: Convert a sparse version of a
     transition matrix into a full matrix.
  }
}



\references{
Broman, K. W. (2005) The genomes of recombinant inbred
lines. \emph{Genetics} \bold{169}, 1133--1146.

Teuscher, F. and Broman, K. W. (2007) Haplotype probabilities for
multiple-strain recombinant inbred lines. \emph{Genetics} \bold{175}, 1267--1274.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{models}
