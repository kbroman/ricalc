#####################################################################
#
# sim_pedigree.R
#
# copyright (c) 2012, Karl W Broman
# last modified Oct, 2012
# first written Oct, 2012
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
# Contains: sim.ped, genAILped, convert2gen
#
######################################################################

######################################################################
# sim.ped: Simulate genotype, using continuous XO locations, for a general pedigree
#
# pedigree = matrix indicating pedigree structure
#     columns = id, sex sire/dad, dam/mom
#     sex must be coded F/M or 0/1 or Female/Male
#     parents must precede offspring
# L = length of chromosome, in cM
# sexsp = female/male recombination rate ratio, must be > 0
# xchr: if TRUE, simulate X chromosome
# m = interference parameter (for chi-square model)
# obligate.chiasma: if TRUE, assume obligate chiasma on 4-strand bundle
######################################################################

sim.ped <-
function(pedigree, L, sexsp=1, xchr=FALSE, m=10, obligate.chiasma=FALSE)
{
  if(L <= 0) stop("L must be > 0")
  if(m < 0) stop("m must be >= 0")
  if(sexsp <= 0) stop("sexsp must be > 0")

  momcol <- match(c("mom", "dam"), colnames(pedigree))
  momcol <- momcol[!is.na(momcol)]
  if(length(momcol) != 1)
    stop("There must be one column labeled \"mom\" or \"dam\".")

  dadcol <- match(c("dad", "sire"), colnames(pedigree))
  dadcol <- dadcol[!is.na(dadcol)]
  if(length(dadcol) != 1)
    stop("There must be one column labeled \"dad\" or \"sire\".")

  sexcol <- grep("^sex$", colnames(pedigree), ignore.case=TRUE)
  if(length(sexcol) != 1)
    stop("Can't find the \"sex\" column.")

  idcol <- grep("^id$", colnames(pedigree), ignore.case=TRUE)
  if(length(idcol) != 1)
    stop("Can't find the \"id\" column.")

  usex <- sort(unique(as.character(pedigree[,sexcol])))
  sex <- as.character(pedigree[,sexcol])
  if(all(usex == c("0", "1")))
    sex <- match(sex, c("0", "1"))-1
  else if(all(toupper(substr(usex, 1, 1)) == c("F", "M")))
    sex <- match(toupper(substr(sex, 1, 1)), c("F", "M")) - 1
  else
    stop("Cannot figure out sex codes: should be 0/1, F/M, female/male, or similar")

  output <- vector("list", nrow(pedigree))
  id <- pedigree[,idcol]
  names(output) <- as.character(id)

  # replace NA mom/dad with 0
  pedigree[is.na(pedigree[,momcol]),momcol] <- 0
  pedigree[is.na(pedigree[,dadcol]),dadcol] <- 0

  if(any((pedigree[,momcol]==0 & pedigree[,dadcol]!=0) |
         (pedigree[,momcol]!=0 & pedigree[,dadcol]==0)))
    stop("Need mom and dad to be both 0/NA or neither 0/NA")

  nextparent <- 1
  for(i in 1:nrow(pedigree)) {

    if(pedigree[i,momcol]==0) {
      output[[i]] <- create.parent(L, nextparent)
      nextparent <- nextparent + 1
    }
    else {
      output[[i]] <- cross(output[[which(id==pedigree[i,momcol])]],
                           output[[which(id==pedigree[i,dadcol])]],
                           m=m, obligate.chiasma=obligate.chiasma, xchr=xchr,
                           male=(sex[i]==1), sexsp=sexsp)
    }

  }

  output
}


######################################################################
# genAILped: generate an AIL pedigree, to F(ngen)
#
# ngen: number of generations (must be > 1)
# npairs: number of mating pairs per generation
# sibship.size: number of offspring per pair in last generation
######################################################################
genAILped <-
function(ngen=8, npairs=20, sibship.size=10)
{
  if(ngen < 2) stop("ngen must be >= 2")
  if(npairs < 1) stop("npairs must be > 0")
  if(sibship.size < 1) stop("sibship size must be > 0")

  output <- matrix(ncol=5, nrow=4+2*npairs*(ngen-2)+sibship.size*npairs)
  colnames(output) <- c("id", "sex", "mom", "dad", "generation")

  # parents
  output[1,] <- c(1, 0, 0, 0, 0)
  output[2,] <- c(2, 1, 0, 0, 0)

  # F1 generation
  output[3,] <- c(3, 0, 1, 2, 1)
  output[4,] <- c(4, 1, 1, 2, 1)

  # F2 generation
  ind <- 4 + (1:(npairs*2))
  output[ind, 1] <- ind
  output[ind, 2] <- rep(0:1, npairs)
  output[ind, 3] <- 3
  output[ind, 4] <- 4
  output[ind, 5] <- 2

  if(ngen==2) return(output)

  for(gen in 3:(ngen-1)) {
    moms <- ind[output[ind,2]==0]
    dads <- ind[output[ind,2]==1]
    if(npairs > 1) {
      moms <- sample(moms)
      dads <- sample(dads)
    }

    ind <- max(c(moms, dads)) + (1:(npairs*2))
    output[ind, 1] <- ind
    output[ind, 2] <- (ind - 1) %% 2
    for(i in seq(along=moms)) {
      output[ind[i*2-1], 3] <- moms[i]
      output[ind[i*2],   3] <- moms[i]
      output[ind[i*2-1], 4] <- dads[i]
      output[ind[i*2],   4] <- dads[i]
    }
    output[ind, 5] <- gen
  }

  moms <- ind[output[ind,2]==0]
  dads <- ind[output[ind,2]==1]
  if(npairs > 1) {
    moms <- sample(moms)
    dads <- sample(dads)
  }

  this <- max(c(moms, dads)) + 1
  for(i in seq(along=moms)) {
    these <- this + (1:sibship.size)-1
    output[these, 1] <- these
    output[these, 2] <- (these - 1) %% 2
    output[these, 3] <- moms[i]
    output[these, 4] <- dads[i]
    output[these, 5] <- ngen
    this <- max(these) + 1
  }

  output
}

######################################################################
# convert2gen: convert the sort of continuous XO location information
#              simulated by sim.ped to marker genotypes
#
# xodat: the sort of detailed genotype/XO data generated by sim.ped
# map:   vector of marker locations
######################################################################

convert2gen <-
function(xodat, map)
{
  if(map[1] != 0) map <- map - map[1]
  if(max(map) > max(xodat[[1]][[1]][1,]))
    warning("maximum simulated position is less than the length of the map.")

  output <- list(matrix(ncol=length(map), nrow=length(xodat)),
                 matrix(ncol=length(map), nrow=length(xodat)))
  for(i in seq(along=xodat)) {
    for(j in 1:2) {
      dat <- xodat[[i]][[j]]

      wh <- sapply(map, function(a,b) max(which(b <= a)), dat[1,])
      output[[j]][i,wh==ncol(dat)] <- dat[2,ncol(dat)]
      output[[j]][i,wh<ncol(dat)] <- dat[2,wh[wh <ncol(dat)]+1]
    }
  }

  if(max(unlist(output)) == 2) { # 2 alleles, so use 1/2/3 codes
    output <- output[[1]] + output[[2]] # becomes 2/3/4
    output[!is.na(output) & output==2] <- 1
    output[!is.na(output) & output==3] <- 2
    output[!is.na(output) & output==4] <- 3
  }
  else # otherwise, use binary codes
    output <- 2^output[[1]] + 2^output[[2]]

  dimnames(output) <- list(names(xodat), names(map))

  output
}

# end of sim_pedigree.R
