######################################################################
#
# ril_prob.R
#
# copyright (c) 2004, Karl W Broman
# last modified Nov, 2004
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
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/ricalc package
# Contains:  get.ril.prob, get.start, count.absorb, get.ril.coinc
#
######################################################################

######################################################################
# get.ril.prob
#
# Get the "absorption" probabilities" for multiway-RILs.
#
######################################################################
get.ril.prob <-
function(r, coinc=1, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"), verbose=TRUE)
{
  n.strains <- match.arg(n.strains)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  n.loci <- match.arg(n.loci)

  if(type=="selfing" && chrtype=="X") {
    chrtype <- "A"
    warning("With type 'selfing', chrtype must be 'A'")
  }

  if(verbose) cat(" -Calculating transition matrix\n")
  fulltm <- get.full.tm(r,coinc,n.strains,type,chrtype,n.loci)
  if(verbose) cat(" -Reforming transition matrix\n")
  fulltm <- convert.full.tm(fulltm)

  start0 <- get.start(n.strains,type,chrtype,n.loci)
  absorb <- count.absorb(n.strains,type,chrtype,n.loci)
  absorbnam <- names(absorb)

  p <- rep(0,length(absorbnam))
  wh <- match(absorbnam,rownames(fulltm))
  A <- fulltm[-wh,-wh]
  A <- A - diag(rep(1,nrow(A)))
  for(i in 1:(length(absorbnam)-1)) {
    if(verbose) cat(" -Solving equation",i,"\n")
    p[i] <- solve(A,-fulltm[-wh,wh[i]])[start0]
  }
  p[length(p)] <- 1-sum(p)
  nam <- absorbnam
  nam <- matrix(unlist(strsplit(nam,"\\|")),ncol=length(nam))[1,]

  names(absorb) <- names(p) <- nam
  attr(p,"counts") <- absorb

  p
}


######################################################################
# get.ril.coinc
#
# Get the "absorption" probabilities" for multiway-RILs.
#
######################################################################
get.ril.coinc <-
function(r, coinc=1, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), verbose=TRUE)
{
  n.strains <- match.arg(n.strains)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)

  if(type=="selfing" && chrtype=="X") {
    chrtype <- "A"
    warning("With type 'selfing', chrtype must be 'A'")
  }


  p2 <- get.ril.prob(r,coinc,n.strains,type,chrtype,"2",FALSE)
  p3 <- get.ril.prob(r,coinc,n.strains,type,chrtype,"3",verbose)

  absorbnam2 <- names(p2)
  absorbnam3 <- names(p3)
  x2 <- matrix(unlist(strsplit(absorbnam2,"")),nrow=2)
  x3 <- matrix(unlist(strsplit(absorbnam3,"")),nrow=3)

  pxo <- sum(p2[x2[1,] != x2[2,]])
  pdblxo <- sum(p3[x3[1,] != x3[2,] & x3[2,] != x3[3,]])

  pdblxo/pxo^2
}


######################################################################
# get.start
#
# determine starting state for the Markov chain
######################################################################
get.start <-
function(n.strains=c("2","4"),type=c("selfing","sibmating"),
         chrtype=c("A","X"),n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  n.loci <- match.arg(n.loci)
  if(type=="selfing" && chrtype=="X") {
    chrtype <- "A"
    warning("With type 'selfing', chrtype must be 'X'.")
  }
  if(type=="selfing" && n.strains=="4") {
    stop("No need to deal with 4 strains by selfing.")
  }

  if(type=="selfing") {
    if(n.loci=="2") { # 2-self-2
      start0 <- "AA|BB"
    }
    else { # 2-self-3
      start0 <- "AAA|BBB"
    }
  }
  else {
    if(n.loci=="2") {
      if(n.strains=="2") {
        if(chrtype=="A") { # 2-sibA-2
          start0 <- "AA|BB x AA|BB"
        }
        else { # 2-sibX-2
          start0 <- "AA|BB x AA"
        }
      }
      else {
        if(chrtype=="A") { # 4-sibA-2
          start0 <- "AA|BB x CC|DD"
        }
        else { # 4-sibX-2
          start0 <- "AA|BB x CC"
        }
      }
    }
    else {
      if(n.strains=="2") {
        if(chrtype=="A") { # 2-sibA-3
          start0 <- "AAA|BBB x AAA|BBB"
        }
        else { # 2-sibX-3
          start0 <- "AAA|BBB x AAA"
        }
      }
      else {
        if(chrtype=="A") { # 4-sibA-3
          start0 <- "AAA|BBB x CCC|DDD"
        }
        else { # 4-sibX-3
          start0 <- "AAA|BBB x CCC"
        }
      }
    }
  }
  start0
}

######################################################################
# count.absorb
#
# Identify all of the absorbing states from a lookup table,
# and count how many of the total states correspond to each
######################################################################
count.absorb <-
function(n.strains=c("2","4"),type=c("selfing","sibmating"),
         chrtype=c("A","X"),n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  n.loci <- match.arg(n.loci)

  data(lookup,envir=environment())
  if(type=="selfing")
    lookup <- lookup[[paste(n.strains,"self",n.loci,sep="")]]
  else
    lookup <- lookup[[paste(n.strains,"sib",chrtype,n.loci,sep="")]]

  z <- x <- unique(lookup)
  if(length(grep(" x ", x))>0) {
    x <- matrix(unlist(strsplit(x," x ")),ncol=2,byrow=TRUE)

    for(i in 1:2) {
      if(length(grep("\\|", x[,i]))>0) {
        y <- matrix(unlist(strsplit(x[,i],"\\|")),ncol=2,byrow=TRUE)
        x[,i] <- y[,1]
        x <- cbind(x,y[,2])
      }
    }
  }
  else {
    if(length(grep("\\|",x))>0) {
      x <- matrix(unlist(strsplit(x,"\\|")),ncol=2,byrow=TRUE)
    }
  }
  wh <- apply(x,1,function(a) length(unique(a))==1)

  sapply(z[wh],function(a,b) sum(a==b),lookup)
}



# end of ril_prob.R
