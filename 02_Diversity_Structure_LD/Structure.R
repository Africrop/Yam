###################################################################################################################################
#
# Copyright 2017 IRD
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
# If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD 
#
# Written by Nora Scarcelli and Olivier Francois
#
###################################################################################################################################


###################################################################################################################################
## This script represents the q matrix created by Structure-like sofwtare as barplot and on a map
####################################################################################################################################

library("tess3r")
library("maps")


as.qmatrix<- function(Q){
  if (class(Q) != "matrix") Q <- as.matrix(Q)
  if (min(Q) < 0) stop("Q contains negative elements.")
  sumofq<- apply(Q, MARGIN = 1, sum)
  if ( sum(sumofq) != nrow(Q)) stop("Input matrix is not an ancestry matrix.")
  class(Q) = "tess3Q"
  return(Q)
}

# Color choice
my.colors =rainbow(3)
my.palette<- CreatePalette(my.colors, 10)

# Read files
coo=read.table("Coordinate.txt") # One row per sample with longitude and latitude
q=read.table("Ancestry.txt") # One row per sample with q matrix 
q <- q/apply(q, 1, sum)
q=as.qmatrix(q)

# Represent the q matrix as barplot
barplot(q, border = NA, space = 0, col.palette = my.palette, sort.by.Q = FALSE) ->bp

# Represent the q matrix on a map
plot(q, coo, method = "map.max", interpol = FieldsKrigModel(10), resolution = c(300,300), cex = .4, col.palette = my.palette, xlab="", ylab="")
map(add = T)

