#####################################################################
#
# full_tm.R
#
# copyright (c) 2004-2012, Karl W Broman
# last modified Oct, 2012
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
# Contains:  get.full.tm, my.kronecker, get.full.tm.symbolic,
#            my.kronecker.symbolic, write.sparse, convert.full.tm
#
######################################################################

######################################################################
# get.full.tm
#
# Get the full generation-to-generation transition matrix for the
# formation of RILs
#
# r = recombination fraction
# coinc = 3-pt coincidence (needed only if n.loci=3)
# n.strains = number of strains (2 or 4)
# type = method of inbreeding
# chrtype = autosomal or X-chromosome (for sib-mating only)
# n.loci = number of loci (2 or 3)
#
######################################################################
get.full.tm <-
function(r, coinc=1, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  n.loci <- match.arg(n.loci)

  if(type=="selfing" && chrtype=="X") {
    chrtype <- "A"
    warning("With type 'selfing', chrtype must be 'A'")
  }

  # per meiosis transition matrix
  mtm <- get.meiosis.tm(r, coinc, n.strains, chrtype, n.loci)

  # parental genotypes
  data(lookup,envir=environment()) # load lookup locally within function
  if(type=="selfing")
    lookup <- lookup[[paste(n.strains,"self",n.loci,sep="")]]
  else
    lookup <- lookup[[paste(n.strains,"sib",chrtype,n.loci,sep="")]]

  # unique parental types
  upg <- unique(lookup)
  n.upg <- length(upg)

  output <- vector("list",n.upg)
  names(output) <- upg

  for(i in 1:n.upg) {
    if(type=="selfing") {
      mp <- mtm[[upg[i]]]
      temp <- my.kronecker(mp,mp)
    }
    else { # sib-mating
      parents <- unlist(strsplit(upg[i]," x "))
      mp.mom <- mtm[[parents[1]]]
      if(chrtype=="A") {
        mp.dad <- mtm[[parents[2]]]
        temp <- my.kronecker(mp.mom,mp.dad,"|")
        temp <- my.kronecker(temp,temp," x ")
      }
      else {
        temp <- mp.mom
        names(temp) <- as.vector(outer(names(temp),parents[2],paste,sep="|"))
        temp <- my.kronecker(temp,mp.mom,sep=" x ")
      }
    }

    # collapse
    names(temp) <- adjust.order(names(temp),chrtype)
    names(temp) <- lookup[names(temp)]

    utempnam <- unique(names(temp))
    n.utemp <- length(utempnam)
    if(length(temp) != n.utemp) {
      utemp <- rep(0,n.utemp)
      names(utemp) <- utempnam
      for(j in 1:n.utemp)
        utemp[j] <- sum(temp[names(temp)==utempnam[j]])
      temp <- utemp
    }

    output[[i]] <- temp
  }

  output
}

######################################################################
# convert.full.tm
#
# convert the sparse versions of the full transition matrix
# (calculated by get.full.tm or get.full.tm.sparse) to the
# explicit matrix form
######################################################################
convert.full.tm <-
function(full.tm)
{
  nam <- names(full.tm)
  if(is.character(full.tm[[1]][1])) zero <- "0"
  else zero <- 0
  out <- lapply(full.tm,function(a,b,d) {
    x <- rep(d,length(b))
    names(x) <- b
    x[names(a)] <- a
    x }, nam,zero)
  out <- matrix(unlist(out),ncol=length(nam),byrow=TRUE)
  dimnames(out) <- list(nam,nam)
  out
}




######################################################################
# my.kronecker
# This is just like the function kronecker, but it ensures that
# the elements' names are correct (pasted with the argument sep
# as a spacer); the kronecker function seems to ignore the names
# completely.
######################################################################
my.kronecker <-
function(a, b, sep="|")
{
  na <- names(a)
  nb <- names(b)
  n.na <- length(a)
  n.nb <- length(b)
  k <- kronecker(a,b)
  names(k) <- t(outer(na,nb,paste,sep=sep))
  k
}

######################################################################
# get.full.tm.symbolic
#
# Symbolic version of get.full.tm
#
# n.strains = number of strains (2 or 4)
# type = method of inbreeding
# chrtype = autosomal or X-chromosome (for sib-mating only)
# n.loci = number of loci (2 or 3)
#
######################################################################
get.full.tm.symbolic <-
function(n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  n.loci <- match.arg(n.loci)
  if(type=="selfing" && chrtype=="X") {
    chrtype <- "A"
    warning("With type 'selfing', chrtype must be 'A'")
  }

  # per meiosis transition matrix
  data(mtm,envir=environment()) # load mtm locally within function
  mtm <- mtm[[paste(n.strains,chrtype,n.loci,sep="")]]

  # parental genotypes
  data(lookup,envir=environment()) # load lookup locally within function
  if(type=="selfing")
    lookup <- lookup[[paste(n.strains,"self",n.loci,sep="")]]
  else
    lookup <- lookup[[paste(n.strains,"sib",chrtype,n.loci,sep="")]]

  # unique parental types
  upg <- unique(lookup)
  n.upg <- length(upg)

  output <- vector("list",n.upg)
  names(output) <- upg

  for(i in 1:n.upg) {
    if(type=="selfing") {
      mp <- mtm[[upg[i]]]
      temp <- my.kronecker.symbolic(mp,mp)
    }
    else { # sib-mating
      parents <- unlist(strsplit(upg[i]," x "))
      mp.mom <- mtm[[parents[1]]]
      if(chrtype=="A") {
        mp.dad <- mtm[[parents[2]]]
        temp <- my.kronecker.symbolic(mp.mom,mp.dad,"|")
        temp <- my.kronecker.symbolic(temp,temp," x ")
      }
      else {
        temp <- mp.mom
        names(temp) <- as.vector(outer(names(temp),parents[2],paste,sep="|"))
        temp <- my.kronecker.symbolic(temp,mp.mom,sep=" x ")
      }
    }

    # collapse
    names(temp) <- adjust.order(names(temp),chrtype)
    names(temp) <- lookup[names(temp)]

    utempnam <- unique(names(temp))
    n.utemp <- length(utempnam)
    if(length(temp) != n.utemp) {
      utemp <- rep(0,n.utemp)
      names(utemp) <- utempnam
      for(j in 1:n.utemp)
        utemp[j] <- paste(temp[names(temp)==utempnam[j]],collapse=" + ")
      temp <- utemp
    }

    output[[i]] <- temp
  }

  output
}

######################################################################
# my.kronecker.symbolic
#
# Symbolic version of my.kronecker
#
######################################################################
my.kronecker.symbolic <-
function(a, b, sep="|")
{
  na <- names(a)
  nb <- names(b)
  n.na <- length(a)
  n.nb <- length(b)
  k <- kronecker(a,b,function(x,y) paste("((",x,")*(",y,"))",sep=" "))
  names(k) <- t(outer(na,nb,paste,sep=sep))
  k
}

######################################################################
# write.meiosis.tm.symbolic:
#     write the results of get.meiosis.tm.symbolic to a file,
#     to be read into mathematica for simplicification.
#
# file = character string for output file
# n.strains = number of strains (2 or 4)
# type = method of inbreeding
# chrtype = autosomal or X-chromosome (for sib-mating only)
# n.loci = number of loci (2 or 3)
#
######################################################################
write.full.tm.symbolic <-
function(file, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  n.loci <- match.arg(n.loci)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  if(type=="selfing" && chrtype=="X") {
    chrtype <- "A"
    warning("With type 'selfing', chrtype must be 'A'")
  }

  output <- get.full.tm.symbolic(n.strains, type, chrtype, n.loci)
  n.output <- length(output)

  write("tm = {",file=file,ncolumns=1,append=FALSE)
  for(i in 1:n.output) {
    string <- paste("{ ",paste(output[[i]],collapse=",")," }",sep="")
    if(i != n.output)
      string <- paste(string,",",sep="")
    write(string,file=file,ncolumns=1,append=TRUE)
  }
  write("}",file=file,ncolumns=1,append=TRUE)
}

######################################################################
# read.full.tm.symbolic:
#     read back in the simplified results from mathematica
#
# file = character string for output file
# n.strains = number of strains (2 or 4)
# type = method of inbreeding
# chrtype = autosomal or X-chromosome (for sib-mating only)
# n.loci = number of loci (2 or 3)
#
######################################################################
read.full.tm.symbolic <-
function(file, n.strains=c("2","4"), type=c("selfing","sibmating"),
         chrtype=c("A","X"), n.loci=c("2","3"))
{
  n.strains <- match.arg(n.strains)
  n.loci <- match.arg(n.loci)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  if(type=="selfing" && chrtype=="X") {
    chrtype <- "A"
    warning("With type 'selfing', chrtype must be 'A'")
  }

  ind <- gtypes(n.strains,"ind",chrtype,n.loci,FALSE,FALSE)

  x <- scan(file,what=character())
  x <- paste(x,collapse="")
  x <- strsplit(x,"[{}]")[[1]]
  x <- x[x!="" & x!=","]
  x <- strsplit(x,",")

  a <- get.full.tm(0.1,1,n.strains,type,chrtype,n.loci)

  names(x) <- names(a)
  for(i in 1:length(x))
    names(x[[i]]) <- names(a[[i]])

  x
}


######################################################################
# write.sparse
#
# write the transition matrix and such in a form suitable for
# import into mathematica, so that we can get symbolic solutions
# of the limiting RIL probabilities
######################################################################
write.sparse <-
function(file=NULL,n.strains=c("2","4"),type=c("selfing","sibmating"),
         chrtype=c("A","X"),n.loci=c("2","3"), where=c("math","perl"))
{
  n.strains <- match.arg(n.strains)
  type <- match.arg(type)
  chrtype <- match.arg(chrtype)
  n.loci <- match.arg(n.loci)
  where <- match.arg(where)
  if(type=="selfing" && chrtype=="X") {
    chrtype <- "A"
    warning("With type 'selfing', chrtype must be 'X'.")
  }
  if(type=="selfing" && n.strains=="4") {
    stop("No need to deal with 4 strains by selfing.")
  }

  data(lookup,envir=environment())
  if(type=="selfing") wh <- paste(n.strains,"self",n.loci,sep="")
  else wh <- paste(n.strains,"sib",chrtype,n.loci,sep="")
  lookup <- lookup[[wh]]
  states <- unique(lookup)

  data(fulltm,envir=environment())
  fulltm <- fulltm[[wh]]

  start0 <- get.start(n.strains,type,chrtype,n.loci)
  absorb <- names(count.absorb(n.strains,type,chrtype,n.loci))

  fulltm <- fulltm[-match(absorb,names(fulltm))]

  rhs <- matrix("",ncol=length(absorb),nrow=length(fulltm))

  names(rhs) <- absorb
  for(i in 1:ncol(rhs)) {
    rhs[,i] <-
      unlist(lapply(fulltm,function(a,b) {
        n <- names(a)
        if(!any(n==b)) return("0")
        else return(paste("-(",a[b],")",sep=""))
      }, absorb[i]))
  }
  fulltm <- lapply(fulltm, function(a,b) {
    x <- match(b,names(a))
    if(any(!is.na(x))) a <- a[-x[!is.na(x)]]
    return(a) }, absorb)
  names(rhs) <- NULL
  dimnames(rhs) <- list(names(fulltm),absorb)

  states <- states[-match(absorb,states)]

  for(i in 1:length(fulltm)) {
    wh <- which(names(fulltm[[i]]) == names(fulltm)[i])
    if(length(wh)==0) {
      fulltm[[i]] <- c(fulltm[[i]],-1)
      names(fulltm[[i]])[length(fulltm[[i]])] <- names(fulltm)[i]
    }
    else {
      fulltm[[i]][wh] <- paste("(",fulltm[[i]][wh],")-1",sep="")
    }
  }

  if(is.null(file))
    return(list(fulltm=fulltm,rhs=rhs,states=states,start0=start0,absorb=absorb))

  if(where == "math") {

    fulltm <- lapply(fulltm,function(a,b) { names(a) <- match(names(a),b)
                                            a}, states)

    write("a = SparseArray[{",file,append=FALSE)
    for(i in 1:length(fulltm)) {
      for(j in 1:length(fulltm[[i]])) {
        string <- paste("{",i,",",names(fulltm[[i]])[j],"} -> ",fulltm[[i]][j],sep="")
        if(i != length(fulltm) || j != length(fulltm[[i]]))
          string <- paste(string,",",sep="")
        write(string, file, append=TRUE)
      }
    }
    write(paste("}, {",length(fulltm),",",length(fulltm),"}]\n",sep=""),
          file, append=TRUE)

    write("b = {",file,append=TRUE)
    for(i in 1:nrow(rhs)) {
      string <- paste("{",paste(rhs[i,],collapse=","),"}",sep="")
      if(i != nrow(rhs)) string <- paste(string,",",sep="")
      write(string,file,append=TRUE)
    }
    write("}\n",file,append=TRUE)

    write(paste("start =", match(start0,names(fulltm))),file,append=TRUE)

  }
  else { # write to perl
    write(paste("start",start0,sep=","),file,append=FALSE)
    write(paste("absorb",paste(absorb,collapse=","),sep=","),file,append=TRUE)

    write(paste("trmatrix",length(fulltm),sep=","),file,append=TRUE)
    for(i in 1:length(fulltm)) {
      write(paste(names(fulltm)[i],length(fulltm[[i]]),sep=","),file,append=TRUE)

      for(j in 1:length(fulltm[[i]]))
        write(paste(names(fulltm[[i]])[j],fulltm[[i]][j],sep=","),file,append=TRUE)

    }

    write(paste("rhs", ncol(rhs),sep=","),file,append=TRUE)
    for(i in 1:ncol(rhs)) {
      write(paste(colnames(rhs)[i],nrow(rhs),sep=","),file,append=TRUE)
      for(j in 1:nrow(rhs))
        write(paste(rownames(rhs)[j],rhs[j,i],sep=","),file,append=TRUE)
    }
  }

  invisible(matrix(unlist(strsplit(absorb,"\\|")),ncol=length(absorb))[1,])
}

# end of full_tm.R
