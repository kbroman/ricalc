\name{ricalc-internal}
\alias{adjust.order}
\alias{alternates}
\alias{create.parent}
\alias{cross}
\alias{exchangeg}
\alias{get.meiosis.tm}
\alias{get.meiosis.tm.symbolic}
\alias{gtypes}
\alias{lengthHet}
\alias{lennum.het}
\alias{meiosis}
\alias{meiosis.sub}
\alias{modifyL}
\alias{my.kronecker}
\alias{my.kronecker.symbolic}
\alias{nubreak}
\alias{read.full.tm.symbolic}
\alias{read.meiosis.tm.symbolic}
\alias{reverseg}
\alias{where.het}
\alias{write.full.tm.symbolic}
\alias{write.meiosis.tm.symbolic}
\alias{write.sparse}

\title{Internal ricalc functions}

\description{
  Internal ricalc functions.
}

\usage{
adjust.order(string, chrtype=c("A","X"))
alternates(pat, n.strains=c("2","4"), chrtype=c("A","X"))
create.parent(L, allele=1)
cross(mom, dad, m=10, obligate.chiasma=FALSE, xchr=FALSE, male=FALSE, sexsp=1)
exchangeg(string,alle=c("A","B"))
get.meiosis.tm(r, coinc=1, n.strains=c("2","4"), chrtype=c("A","X"),
         n.loci=c("2","3"))
get.meiosis.tm.symbolic(n.strains=c("2","4"), chrtype=c("A","X"),
         n.loci=c("2","3"))
gtypes(n.strains=c("2","4"), type=c("sibpair","ind","hap"),
         chrtype=c("A","X"), n.loci=c("2","3"), use.symmetry=TRUE,
         verbose=TRUE)
lengthHet(chr)
lennum.het(chr)
meiosis(parent, m=10, obligate.chiasma=FALSE)
meiosis.sub(L, m=10, obligate.chiasma=FALSE)
modifyL(ind, factor)
my.kronecker(a, b, sep="|")
my.kronecker.symbolic(a, b, sep="|")
nubreak(ind)
read.full.tm.symbolic(file, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"))
read.meiosis.tm.symbolic(file, n.strains=c("2","4"), chrtype=c("A","X"),
         n.loci=c("2","3"))
reverseg(string)
where.het(ind)
write.full.tm.symbolic(file, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"))
write.meiosis.tm.symbolic(file, n.strains=c("2","4"), chrtype=c("A","X"),
         n.loci=c("2","3"))
write.sparse(file=NULL,n.strains=c("2","4"),type=c("selfing","sibmating"),
         chrtype=c("A","X"),n.loci=c("2","3"), where=c("math","perl"))
}

\details{
  These are generally not to be called by the user.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{internal}
