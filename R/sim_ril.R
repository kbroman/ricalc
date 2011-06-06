#####################################################################
#
# sim_ril.R
#
# copyright (c) 2004-9, Karl W Broman
# last modified Apr, 2009
# first written May, 2004
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/ricalc package
# Contains: meiosis.sub, meiosis, create.parent, where.het,
#           lengthHet, sim.ri, cross, lennum.het, modifyL, nubreak
#
######################################################################

######################################################################
# meiosis.sub 
#
# Simulate the locations of crossovers on a meiotic product
# via the chi-square model (m=0 corresponds to no interference)
# and potentially with an obligate chiasma.
#
# L = chromosome length (in cM)
# m = interference parameter (m=0 is NI, m>0 is positive interference)
# obligate.chiasma = if TRUE, obligate chiasma on four-strand bundle
######################################################################
meiosis.sub <-
function(L, m=10, obligate.chiasma=FALSE)
{
  if(obligate.chiasma) { # adjust mean no. chiasmata
    if(L <= 50) stop("L must be > 50 cM")
    if(m==0) f <- function(Lstar,f.L,f.m=0) f.L-Lstar/(1-exp(-Lstar/50)) 
    else {
      f <- function(Lstar,f.L,f.m=0)
        {
          lambdastar <- Lstar/50*(f.m+1)
          temp <- lambdastar
          for(i in 1:length(temp))
            temp[i] <- sum(exp(-lambdastar[i] + (0:f.m)*log(lambdastar[i])-
                               lgamma((0:f.m)+1)) * (f.m+1-(0:f.m))/(f.m+1))
          f.L - Lstar/(1-temp)
        }
    }

    Lstar <- uniroot(f,c(1e-5,L+1),f.m=m,f.L=L)$root
  }
  else Lstar <- L

  if(m==0) { # no interference
    if(!obligate.chiasma)  # no obligate chiasma
      n.xo <- rpois(1,Lstar/100)
    else {
      up <- qpois(1e-14,Lstar/50,lower=FALSE)
      p <- dpois(1:up,Lstar/50)/ppois(0,Lstar/50)
      n.chi <- sample(1:up,1,prob=p)
      n.xo <- rbinom(1,n.chi,0.5)
    }
    if(n.xo==0) xo <- NULL
    else xo <- sort(runif(n.xo,0,L))
  }
  else { # chi-square model
    n.chi <- 0
    while(n.chi == 0) {
      n.pts <- rpois(1,Lstar/50*(m+1))
      first <- sample(1:(m+1),1)
      if(first <= n.pts || !obligate.chiasma) n.chi <- 1
    }
    if(first > n.pts)
      xo <- NULL
    else {
      pt.loc <- sort(runif(n.pts,0,L))
      chi.loc <- pt.loc[seq(first,length(pt.loc),by=m+1)]
      n.xo <- rbinom(1,length(chi.loc),0.5)
      if(n.xo==0) xo <- NULL
      else if(length(chi.loc)==1) xo <- chi.loc
      else xo <- sort(sample(chi.loc,n.xo,repl=FALSE))
    }
  }
    
  if(length(xo) == 0) xo <- NULL
  xo
}

######################################################################
# create.parent
#
# create a parental individual
# 
# L = chromosome length
# allele = vector of length 1 or two 
######################################################################
create.parent <-
function(L, allele=1)
{
  if(length(allele) == 1) allele <- rep(allele,2)
  if(length(allele) != 2)
    stop("allele should be of length 1 or 2")
  
  list(mat=rbind(c(0,L),allele[1]),
       pat=rbind(c(0,L),allele[2]))
}

######################################################################
# meiosis
#
# Output a random meiotic product from an input individual.
######################################################################
meiosis <-
function(parent, m=10, obligate.chiasma=FALSE)
{
  L <- parent$mat[1,ncol(parent$mat)]
  if(abs(parent$pat[1,ncol(parent$pat)] - L) > 1e-13)
    stop("There is a problem with the parent's data structure.")

  product <- meiosis.sub(L, m, obligate.chiasma)
  a <- sample(1:2,1)

  if(length(product)==0) return(parent[[a]])

  else {
    for(i in 1:length(product)) {
      if(i == 1) 
        result <- parent[[a]][,parent[[a]][1,]<product[1],drop=FALSE]
      else {
        temp <- parent[[a]][1,]>=product[i-1] & parent[[a]][1,]<product[i]
        result <- cbind(result,parent[[a]][,temp])
      }
      u <- parent[[a]][2,parent[[a]][1,]>=product[i]]
      result <- cbind(result,c(product[i],u[1]))
      a <- 3-a
    }
    temp <- parent[[a]][1,]>=product[length(product)]
    result <- cbind(result,parent[[a]][,temp])
  }

  # clean out excess stuff in the result
  if(ncol(result)>2) {
    keep <- rep(TRUE,ncol(result))
    for(i in 2:(ncol(result)-1)) 
      if(result[2,i] == result[2,i+1])
        keep[i] <- FALSE
  }
#  print(result)
  result[,keep,drop=FALSE]
}

######################################################################
# cross
#
# cross two individuals to create a single progeny
######################################################################
cross <-
function(mom, dad, m=10, obligate.chiasma=FALSE, xchr=FALSE, male=FALSE,
         sexsp=1)
{
  if(sexsp != 1 && !xchr) {
    mom <- modifyL(mom,sexsp*2/(1+sexsp))
    dad <- modifyL(dad,2/(1+sexsp))
  }
                   
  mat <- meiosis(mom,m,obligate.chiasma)

  if(!xchr) 
    pat <- meiosis(dad,m,obligate.chiasma)
  else {
    if(male) 
      pat <- dad$pat
    else
      pat <- dad$mat
  }

  if(sexsp != 1 && !xchr) {
    mat <- modifyL(mat,(1+sexsp)/(sexsp*2))
    pat <- modifyL(pat,(1+sexsp)/2)
  }

  list(mat=mat,pat=pat)
}

######################################################################
# sim.ri
#
# simulate a recombinant inbred line
#
# L = chromomosome length (can be a vector)
# sexsp = female:male recombination rate (number > 0)
# type = by selfing or sib-mating
# n.strains = 2, 4, or 8 initial parental strains
# xchr = simulate X chromosome (when type=sib only)
#        If L is a vector, this is ignored, and names(L) is used...
#        chr with name="X" or ="x" is the X, others are asssumed
#        to be autosomes
# m = interference parameter
# obligate.chiasma = if TRUE, obligate chiasma on 4-strand bundle
#
# OUTPUT: If length(L)=1, output is a 2-row matrix; first row
#         gives locations of crossovers and second row gives alleles
#         in each segment
#
#         If length(L)>1, output is a list of 2-row matrices, one
#         for each chromosome.
######################################################################
sim.ri <-
function(L, sexsp=1, type=c("selfing","sibmating"),
         n.strains=c("2","4","8"), 
         xchr=FALSE, m=10, obligate.chiasma=FALSE)
{  
  n.strains <- match.arg(n.strains)
  type <- match.arg(type)

  if(sexsp <= 0) stop("sexsp must be > 0")

  if(length(L) > 1) { # more than one chr: use recursion
    output <- vector("list",length(L))
    if(is.null(names(L))) names(L) <- 1:length(L)
    names(output) <- names(L)
    for(i in 1:length(L)) {
      if(names(L)[i]=="X" || names(L)[i]=="x") {
        output[[i]] <- sim.ri(L[i],sexsp,type,n.strains,
                              xchr=TRUE,m,obligate.chiasma)
      }
      else {
        output[[i]] <- sim.ri(L[i],sexsp,type,n.strains,
                              xchr=FALSE,m,obligate.chiasma)
      }
    }

    # get overall prop.het
    lhet <- lapply(output, attr, "prop.het")
    for(i in 1:length(lhet)) lhet[[i]] <- lhet[[i]]*L[i]
    maxgen <- max(sapply(lhet, length))
    lhet <- matrix(unlist(lapply(lhet, function(a,b) {
      if(length(a) < b) a <- c(a,rep(0,b-length(a))); a }, maxgen)),
                   ncol=length(lhet))
    lhet <- apply(lhet,1,sum)/sum(L)
    attr(output, "prop.het") <- lhet

    # get overall number of heterozygous segments
    nhet <- lapply(output, attr, "num.het")
    nhet <- matrix(unlist(lapply(nhet,function(a,b) {
      if(length(a) < b) a <- c(a,rep(0,b-length(a))); a }, maxgen)),
                   ncol=length(nhet))
    nhet <- apply(nhet,1,sum)
    attr(output,"num.het") <- nhet
      
    # get overall number of breaks at each generation
    nubreak <- lapply(output, attr, "nubreak")
    n <- max(sapply(nubreak, length))
    nubreak <- matrix(unlist(lapply(nubreak, function(a,b) {
      if(length(a) < b) a <- c(a,rep(a[length(a)],b-length(a))); a }, n)),
                      nrow=n)
    nubreak <- apply(nubreak, 1, sum)
    attr(output,"nubreak") <- nubreak

    return(output)
  }

  # first bits of inter-breeding
  if(n.strains=="2") {
    par1 <- create.parent(L,1)
    par2 <- create.parent(L,2)
  }
  else if(n.strains=="4") {
    par1 <- create.parent(L,1:2)
    par2 <- create.parent(L,3:4)
  }
  else {
    par1 <- cross(create.parent(L,1:2), create.parent(L,3:4),
                  m, obligate.chiasma, xchr, male=FALSE,sexsp)
    par2 <- cross(create.parent(L,5:6), create.parent(L,7:8),
                  m, obligate.chiasma, xchr, male=TRUE,sexsp)
  }

  thenubreak <- nubreak(par1)

  # final inter-breeding
  sib1 <- cross(par1, par2, m, obligate.chiasma, xchr, male=FALSE, sexsp)
  if(type=="sibmating") 
    sib2 <- cross(par1, par2, m, obligate.chiasma, xchr, male=FALSE, sexsp)
  else sib2 <- sib1
  par1 <- sib1; par2 <- sib2

  thenubreak <- c(thenubreak, nubreak(par1))
  nhet <- lhet <- NULL
  
  while(1) { # stop loop when chromosome is fixed
    sib1 <- cross(par1, par2, m, obligate.chiasma, xchr, male=FALSE, sexsp)
    if(type=="sibmating") 
      sib2 <- cross(par1, par2, m, obligate.chiasma, xchr, male=FALSE, sexsp)
    else sib2 <- sib1

    if(type=="sibmating") 
      temp <- lennum.het(c(sib1,sib2))
    else
      temp <- lennum.het(sib1)
    
    lhet <- c(lhet,temp[1])
    nhet <- c(nhet,temp[2])

    thenubreak <- c(thenubreak, nubreak(sib1))

    if(lhet[length(lhet)] == 0) break

    par1 <- sib1; par2 <- sib2
  }
  
  output <- sib1[[1]]
  attr(output,"prop.het") <- lhet/L
  attr(output,"num.het") <- nhet
  attr(output, "nubreak") <- thenubreak

  output
}

##############################
# where.het
#
# find regions of heterozygosity
# in an individual
##############################
where.het <-
function(ind)
{
  if(ncol(ind$mat)==ncol(ind$pat) && all(ind$mat == ind$pat)) {
#    cat(" --No regions of heterozygosity\n")
    return(NULL)
  }
  u <- sort(unique(c(ind$mat[1,],ind$pat[1,])))
  het <- NULL
  for(i in 2:length(u)) {
    mat <- ind$mat[,ind$mat[1,] >= u[i],drop=FALSE]
    mat <- mat[2,1]

    pat <- ind$pat[,ind$pat[1,] >= u[i],drop=FALSE]
    pat <- pat[2,1]

    if(mat!=pat) { # heterozygous
      if(is.null(het)) het <- cbind(u[i-1],u[i])
      else het <- rbind(het,c(u[i-1],u[i]))
    }
  }

  # clean up
  if(nrow(het) > 1) {
    keep <- rep(TRUE,nrow(het))
    for(j in 2:nrow(het)) {
      if(het[j,1] == het[j-1,2]) {
        het[j,1] <- het[j-1,1]
        keep[j-1] <- FALSE
      }
    }
    het <- het[keep,,drop=FALSE]
  }
  het
}



##############################
# lengthHet
#
# length of heterozygosity
# in a set of chromosomes
##############################
lengthHet <-
function(chr)
{
  if(length(chr) == 1) stop("Need more than one chromosome.\n")
  if(length(unique(sapply(chr,ncol)))==1) {
    flag <- 0
    for(i in 2:length(chr)) {
      if(any(chr[[1]] != chr[[i]])) {
        flag <- 1
        break
      }
    }
    if(flag==0) return(0)
  }

  # find all breakpoints
  u <- sort(unique(unlist(lapply(chr,function(a) a[1,]))))

  het <- 0
  for(i in 2:length(u)) {
    temp <- sapply(chr,function(a,b) a[2,a[1,]>=b][1],u[i])
    if(length(unique(temp)) > 1)
      het <- het + u[i]-u[i-1]
  }

  het
}

######################################################################
# lennum.het
#
# length of heterozygosity and the number of het segments
# in a set of chromosomes
######################################################################
lennum.het <-
function(chr)
{
  if(length(chr) == 1) stop("Need more than one chromosome.\n")
  if(length(unique(sapply(chr,ncol)))==1) {
    flag <- 0
    for(i in 2:length(chr)) {
      if(any(chr[[1]] != chr[[i]])) {
        flag <- 1
        break
      }
    }
    if(flag==0) return(c(0,0))
  }

  # find all breakpoints
  u <- sort(unique(unlist(lapply(chr,function(a) a[1,]))))

  nhet <- lhet <- 0
  for(i in 2:length(u)) {
    temp <- sapply(chr,function(a,b) a[2,a[1,]>=b][1],u[i])
    if(length(unique(temp)) > 1) {
      lhet <- lhet + u[i]-u[i-1]
      nhet <- nhet + 1
    }
  }

  c(lhet,nhet)
}

######################################################################
# modifyL
#
# inflate/deflate length of chromosome to accommodate sex-specific
# recombination rates
######################################################################
modifyL <-
function(ind, factor)
{
  if(is.list(ind)) 
    return(lapply(ind, function(a,b) { a[1,] <- a[1,]*b; a }, factor))

  ind[1,] <- ind[1,]*factor
  
  ind  
}

######################################################################
# number of unique breakpoints; assume individual has one chromosome
######################################################################
nubreak <-
function(ind)
  length(unique(sort(unlist(lapply(ind, function(a) a[1,])))))-2
    


# end of ricalc.R


