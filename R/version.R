#####################################################################
#
# version.R
#
# copyright (c) 2011, Karl W Broman
# last modified June, 2011
# first written June, 2011
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
######################################################################
# print the installed version of R/ricalc
######################################################################
ricalcversion <-
function()
{
  u <- strsplit(library(help=ricalc)[[3]][[1]][2]," ")[[1]]
  u[length(u)]
}

# end of version.R
