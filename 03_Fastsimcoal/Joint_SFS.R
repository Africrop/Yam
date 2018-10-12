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
# Intellectual property belongs to IRD and Grenoble-Alpes University
#
# Written by Philippe Cubry and Roland Akakpo
#
###################################################################################################################################

###################################################################################################################################
# R script to produce multi-dimensional SFS 
# Uses file produced by SFS_by_pop.R
###################################################################################################################################

library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(cowplot)
library(LEA)
library(grid)

## Read data
sfs.bypop <- read.table("sfs.bypop.txt",header=TRUE)

## Getting population sizes
pop.size <- apply((na.exclude(sfs.bypop)),MARGIN = 2,max)
names(pop.size)<- colnames(sfs.bypop)

## Producing multi-dimensional SFS according to the number of subgroups
## For fastsimcoal, the pop "0" will be the first called in factors list, if n population were analysed, pop "n-1"  wil be the last one
## Here is the example of a 3D SFS calculated  between the 3 groups
## joint_SFS for K=3 pop. The population are pop1, pop2 and pop3
joint_SFS<- (as.data.frame(na.omit(table((sfs.bypop$pop1), 
                                                (sfs.bypop$pop2),
                                                (sfs.bypop$pop3),
                                                 dnn = c("pop1","pop2","pop3")
))))
                                               

write.table(joint_SFS , "joint_SFS")
write.table("1 observations. No. of demes and sample sizes are on next line","joint_SFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("3",pop.size["pop1"],pop.size["pop2"],pop.size["pop3"])),"joint_SFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(joint_SFS$Freq),"joint_SFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)




