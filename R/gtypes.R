#####################################################################
#
# gtypes.R
#
# copyright (c) 2004, Karl W Broman
# last modified May, 2004
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
#     at https://www.r-project.org/Licenses/GPL-3
#
# Part of the R/ricalc package
# Contains: gtypes, alternates, reverseg, exchangeg, adjust.order
#
######################################################################

######################################################################
# gtypes: generate the possible genotype patterns
#
# n.strains = number of strains (2 or 4)
# type = sibpair: genotypes of a sib pair
#      = ind: individual genotypes
#      = hap: haplotypes
# chrtype = A or X (autosome or X chromosome)
# n.loci = number of loci (2 or 3)
#
######################################################################
gtypes <-
function(n.strains=c("2","4"), type=c("sibpair","ind","hap"),
         chrtype=c("A","X"), n.loci=c("2","3"), use.symmetry=TRUE,
         verbose=TRUE)
{
  n.strains <- match.arg(n.strains)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  n.loci <- match.arg(n.loci)

  if(chrtype=="X" && n.strains=="4")
    alleles <- LETTERS[1:3]
  else alleles <- LETTERS[1:as.numeric(n.strains)]

  na <- length(alleles)

  if(verbose) cat(" -Creating hap\n")
  # the possible haplotypes
  if(n.loci=="3") {
    hap <- as.vector(t(outer(alleles,alleles,paste,sep="")))
    hap <- as.vector(t(outer(hap,alleles,paste,sep="")))
  }
  else
    hap <- as.vector(t(outer(alleles,alleles,paste,sep="")))

  if(type=="hap") return(hap)


  if(verbose) cat(" -Creating ind\n")
  # the possible individuals
  ind <- outer(hap,hap,paste,sep="|")
  ind <- sort(ind[row(ind) <= col(ind)])

  if(type=="ind") {
    if(use.symmetry) {
      n.ind <- length(ind)
      uind <- rep("",n.ind)
      names(uind) <- ind
      for(i in 1:n.ind) {
        if(uind[i] == "") {
         temp <- alternates(ind[i], n.strains, chrtype)
         uind[temp] <- temp[1]
       }
      }
      ind <- uind
    }
    return(ind)
  }

  if(verbose) cat(" -Creating pairs\n")
  if(chrtype=="X")
    pair <- as.vector(t(outer(ind,hap,paste,sep=" x ")))
  else {
    if(n.strains=="2" || n.loci=="2") {
      pair <- outer(ind,ind,paste,sep=" x ")
      pair <- sort(pair[row(pair) <= col(pair)])
    }
    else { # in the 4-allele, 3 locus case, do this a bit at a time
      first.half <- 1:(length(ind)/2)
      pair <- outer(ind[first.half],ind[first.half],paste,sep=" x ")
      pair <- pair[row(pair) <= col(pair)]
      pair2 <- outer(ind[-first.half],ind[-first.half],paste,sep=" x ")
      pair2 <- pair2[row(pair2) <= col(pair2)]
      pair3 <- outer(ind[first.half],ind[-first.half],paste,sep=" x ")
      pair <- c(pair,pair2,pair3)
      rm(pair2,pair3,envir=environment())
      if(verbose) cat(" --Sorting\n")
      pair <- sort(pair)
    }
  }

  if(use.symmetry) {
    if(verbose) cat(" -Reducing pairs\n")
    n.pair <- length(pair)
    upair <- rep("",n.pair)
    names(upair) <- pair
    for(i in 1:n.pair) {
      if(verbose && i==round(i,-2)) cat(" --",i,"out of ",n.pair,"\n")
      if(upair[i] == "") {
        temp <- alternates(pair[i], n.strains, chrtype)
        upair[temp] <- temp[1]
      }
    }
    pair <- upair
  }

  pair
}


######################################################################
# alternates
#
# Uses symmetry to find all patterns equivalent to the input one.
######################################################################
alternates <-
function(pat, n.strains=c("2","4"), chrtype=c("A","X"))
{
  n.strains <- match.arg(n.strains)
  chrtype <- match.arg(chrtype)

  if(n.strains=="2" && chrtype=="X") {
    pat <- adjust.order(c(pat,reverseg(pat)),chrtype)
    return(sort(unique(pat)))
  }
  else if(n.strains=="2" || chrtype=="A") {
    pat2 <- exchangeg(pat)
    pat <- adjust.order(c(pat,pat2,reverseg(pat),reverseg(pat2)),chrtype)
    return(sort(unique(pat)))
  }
  else {
    pat2 <- exchangeg(pat)
    pat3 <- exchangeg(pat,c("C","D"))
    pat4 <- exchangeg(pat2,c("C","D"))
    pat5 <- exchangeg(exchangeg(pat,c("A","C")),c("B","D"))
    pat6 <- exchangeg(pat5)
    pat7 <- exchangeg(pat6,c("C","D"))
    pat8 <- exchangeg(pat5,c("C","D"))
    pat <- adjust.order(c(pat,pat2,pat3,pat4,pat5,pat6,pat7,pat8,
                          reverseg(pat),reverseg(pat2),reverseg(pat3),
                          reverseg(pat4),reverseg(pat5),reverseg(pat6),
                          reverseg(pat7),reverseg(pat8)), chrtype)
    return(sort(unique(pat)))
  }
}


######################################################################
# reverseg
#
# Reverse the locus order for a pattern
######################################################################
reverseg <-
function(string)
{
  if(length(grep(" x ",string))==0) { # one individual
    string <- unlist(strsplit(string,"\\|"))
    o <- strsplit(string,"")
    return(paste(sapply(o, function(a) paste(rev(a),collapse="")),
                 collapse="|"))
  }
  else { # a sib pair
    string <- unlist(strsplit(string," x "))
    for(i in 1:length(string)) {
      if(length(grep("\\|",string[i]))>0)
        string[i] <- paste(sapply(strsplit(unlist(strsplit(string[i],"\\|")),""),
                                  function(a) paste(rev(a),collapse="")),collapse="|")
      else
        string[i] <- paste(rev(unlist(strsplit(string[i],""))),collapse="")
    }
    return(paste(unlist(string),collapse=" x "))
  }
}

######################################################################
# exchangeg
#
# Exchange two alleles in a pattern
######################################################################
exchangeg <-
function(string,alle=c("A","B"))
{
  out <- o <- unlist(strsplit(string,""))
  out[o==alle[1]] <- alle[2]
  out[o==alle[2]] <- alle[1]
  paste(out,collapse="")
}

######################################################################
# adjust.order
#
# Adjust the order of haplotypes / individuals in a pattern
######################################################################
adjust.order <-
function(string, chrtype=c("A","X"))
{
  chrtype <- match.arg(chrtype)

  if(length(grep(" x ", string))==0) { # one individual
    z <- matrix(unlist(strsplit(string, "\\|")),ncol=2,byrow=TRUE)
    z <- apply(z, 1, function(a) paste(sort(a),collapse="|"))
  }
  else {
    z <- matrix(unlist(strsplit(string, " x ")),ncol=2,byrow=TRUE)

    zz1 <- matrix(unlist(strsplit(z[,1], "\\|")),ncol=2,byrow=TRUE)
    z[,1] <- apply(zz1, 1, function(a) paste(sort(a),collapse="|"))

    if(chrtype == "A") {
      zz2 <- matrix(unlist(strsplit(z[,2], "\\|")),ncol=2,byrow=TRUE)
      z[,2] <- apply(zz2, 1, function(a) paste(sort(a),collapse="|"))
    }

    if(chrtype == "A")
      z <- apply(z,1,function(a) paste(sort(a),collapse=" x "))
    else z <- apply(z,1,paste,collapse=" x ")
  }
  z
}

# end of gtypes.R
