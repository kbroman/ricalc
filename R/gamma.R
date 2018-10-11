#####################################################################
#
# gamma.R
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
#     at https://www.r-project.org/Licenses/GPL-3
#
# Part of the R/ricalc package
# Contains: mf.gam, imf.gam, coinc.gam
#
######################################################################

######################################################################
# mf.gam
#
# Map function for the gamma renewal model (no obligate chiasma)
#
# d = interval length (in cM)
# nu = interference parameter (= m+1); nu=0 is NI
# tol = tolerance for doing numerical integration
######################################################################
mf.gam <-
function(d, nu=1, tol=1e-12)
{
  if(any(d < 0)) stop("Must have d >= 0")

  if(nu==1) return(0.5*(1-exp(-d/50)))

  d <- d/100
  mf.gam.sub <- function(y,NU)
    pgamma(y,shape=NU,rate=2*NU,lower.tail=FALSE)

  z <- d
  for(i in seq(along=d))
    z[i] <- integrate(mf.gam.sub,lower=0,upper=d[i],NU=nu,rel.tol=tol)$value

  z
}

######################################################################
# imf.gam
#
# Inverse map function for the gamma renewal model (no obl chiasma)
#
# r = recombination fraction
# nu = interference parameter
# tol = tolerance for doing numerical integral and finding root
######################################################################
imf.gam <-
function(r, nu=1, tol=1e-12)
 {
   if(any(r < 0 | r >= 0.5)) stop("Must have 0 <= r < 0.5")

   if(nu==1) return(-50*log(1-2*r))

   imf.gam.sub <- function(y,NU,R,TOL)
     mf.gam(y,NU,TOL)-R

   z <- r
   for(i in seq(along=r)) {
     if(r[i]-mf.gam(r[i]*100,nu) < tol) z[i] <- r[i]*100
     else  {
       up <- -50*log(1-2*r[i])
       lo <- r[i]*100
       z[i] <- uniroot(imf.gam.sub,upper=up,lower=lo,NU=nu,R=r[i],TOL=tol,tol=tol)$root
     }
   }

  z
}

######################################################################
# coinc.gam
#
# The three-point coindence for the gamma model, for the case
# that the two defined intervals have the same recombination
# fraction
######################################################################
coinc.gam <-
function(r, nu=1, tol=1e-12)
{
  if(any(r<=0 || r>0.5))
    stop("Must have 0 < r <= 0.5")

  out <- r
  wh <- (r==0.5)
  out[wh] <- 1
  out[!wh] <- (1-mf.gam(imf.gam(r[!wh],nu,tol)*2,nu,tol)/(2*r[!wh]))/r[!wh]

  out
}

# end of gamma.R
