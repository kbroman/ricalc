#####################################################################
#
# meiosis_tm.R
#
# copyright (c) 2004-7, Karl W Broman
# last modified Oct, 2007
# first written May, 2004
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
# Part of the R/ricalc package
# Contains:  get.meiosis.tm, get.meiosis.tm.symbolic,
#            write.meiosis.tm.symbolic, read.meiosis.tm.symbolic
#
######################################################################

######################################################################
# get.meiosis.tm: get the transition matrix for one meiosis
#                 parent genotype -> haplotype in meiotic product
#
# r = recombination fraction
# coinc = 3-pt coincidence (needed only if n.loci=3)
# n.strains = number of strains (2 or 4)
# chrtype = autosomal or X-chromosome
# n.loci = number of loci (2 or 3)
#
######################################################################
get.meiosis.tm <- 
function(r, coinc=1, n.strains=c("2","4"), chrtype=c("A","X"),
         n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  n.loci <- match.arg(n.loci)
  chrtype <- match.arg(chrtype)

  ind <- gtypes(n.strains,"ind",chrtype,n.loci,FALSE,FALSE)
  chr <- lapply(strsplit(ind,"\\|"),strsplit,"")

  probs <-
  function(x, r, coinc, n.loci=c("2","3"))
  {
    n.loci <- match.arg(n.loci)
    y <- x[[1]]
    z <- x[[2]]
    if(n.loci=="2") {
      nam <- list(y, z, c(y[1],z[2]), c(z[1],y[2]))
      res <- c((1-r),(1-r),r,r)/2
    }
    else {
      nam <- list(y, z, c(y[1],z[2:3]),
                  c(z[1],y[2:3]), c(y[1:2],z[3]),
                  c(z[1:2],y[3]), c(y[1],z[2],y[3]),
                  c(z[1],y[2],z[3]))

      res <- rep(c(1-r-r+coinc*r*r,
                   r*(1-coinc*r), r*(1-coinc*r),
                   coinc*r*r)/2, rep(2,4))
    }
    nam <- sapply(nam,paste,collapse="")

    u <- unique(nam)
    finalres <- rep(0,length(u))
    names(finalres) <- u
    for(i in seq(along=u))
      finalres[i] <- sum(res[nam==u[i]])
    finalres
  }

  output <- lapply(chr, probs, r, coinc, n.loci)
  names(output) <- ind
  output

}

######################################################################
# get.meiosis.tm.symbolic:
#     just like get.meiosis.tm, but giving symbolic results
#     (as character strings)
#
# n.strains = number of strains (2 or 4)
# chrtype = autosomal or X-chromosome
# n.loci = number of loci (2 or 3)
#
######################################################################
get.meiosis.tm.symbolic <- 
function(n.strains=c("2","4"), chrtype=c("A","X"),
         n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  n.loci <- match.arg(n.loci)
  chrtype <- match.arg(chrtype)

  ind <- gtypes(n.strains,"ind",chrtype,n.loci,FALSE,FALSE)
  chr <- lapply(strsplit(ind,"\\|"),strsplit,"")

  probs <-
  function(x, n.loci=c("2","3"))
  {
    n.loci <- match.arg(n.loci)
    y <- x[[1]]
    z <- x[[2]]
    if(n.loci=="2") {
      nam <- list(y, z, c(y[1],z[2]), c(z[1],y[2]))
      res <- c("((1-r)/2)","((1-r)/2)","(r/2)","(r/2)")
    }
    else {
      nam <- list(y, z, c(y[1],z[2:3]),
                  c(z[1],y[2:3]), c(y[1:2],z[3]),
                  c(z[1:2],y[3]), c(y[1],z[2],y[3]),
                  c(z[1],y[2],z[3]))

      res <- rep(c("((1-2*r+c*r*r)/2)",
                   "(r*(1-c*r)/2)", "(r*(1-c*r)/2)",
                   "(c*r^2/2)"), rep(2,4))
    }
    nam <- sapply(nam,paste,collapse="")

    u <- unique(nam)
    finalres <- rep(0,length(u))
    names(finalres) <- u
    for(i in seq(along=u))
      finalres[i] <- paste(res[nam==u[i]],collapse="+")
    finalres
  }

  output <- lapply(chr, probs, n.loci)
  names(output) <- ind
  output

}

######################################################################
# write.meiosis.tm.symbolic:
#     write the results of get.meiosis.tm.symbolic to a file,
#     to be read into mathematica for simplicification.
#
# file = character string for output file
# n.strains = number of strains (2 or 4)
# chrtype = autosomal or X-chromosome
# n.loci = number of loci (2 or 3)
#
######################################################################
write.meiosis.tm.symbolic <- 
function(file, n.strains=c("2","4"), chrtype=c("A","X"),
         n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  n.loci <- match.arg(n.loci)
  chrtype <- match.arg(chrtype)

  output <- get.meiosis.tm.symbolic(n.strains, chrtype, n.loci)
  n.output <- length(output)

  write("tm = {",file=file,ncol=1,append=FALSE)
  for(i in 1:n.output) {
    string <- paste("{ ",paste(output[[i]],collapse=",")," }",sep="")
    if(i != n.output)
      string <- paste(string,",",sep="")
    write(string,file=file,ncol=1,append=TRUE)
  }
  write("}",file=file,ncol=1,append=TRUE)
}

######################################################################
# read.meiosis.tm.symbolic:
#     read back in the simplified results from mathematica
#
# file = character string for output file
# n.strains = number of strains (2 or 4)
# chrtype = autosomal or X-chromosome
# n.loci = number of loci (2 or 3)
#
######################################################################
read.meiosis.tm.symbolic <- 
function(file, n.strains=c("2","4"), chrtype=c("A","X"),
         n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  n.loci <- match.arg(n.loci)
  chrtype <- match.arg(chrtype)

  ind <- gtypes(n.strains,"ind",chrtype,n.loci,FALSE,FALSE)

  x <- scan(file,what=character())
  x <- paste(x,collapse="")
  x <- strsplit(x,"[{}]")[[1]]
  x <- x[x!="" & x!=","]
  x <- strsplit(x,",")

  a <- get.meiosis.tm(0.1,1,n.strains,chrtype,n.loci)

  names(x) <- names(a)
  for(i in 1:length(x))
    names(x[[i]]) <- names(a[[i]])

  x
}

# end of meiosis_tm.R
