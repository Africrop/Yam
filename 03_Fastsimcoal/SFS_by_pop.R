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
# R script to produce SFS per group
# Requires a 012 format file
###################################################################################################################################


library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(cowplot)
library(LEA)
library(grid)
library(data.table)

## Reading data file
## Data file column names = List\t  Name\t  ID\t  Species\t  Country
data=as.data.frame(t(read.table("In.012", sep="\t", row.names=1, h=T)))

## Randomly select 100000 snp
data=data[, sample(1:ncol(data), 100000, replace=FALSE)] 
unique(as.numeric(as.matrix(data)))
head(data[,1:6])

## Add information about groups
ind = as.data.frame(read.table("Individual_Information.txt", h = T, sep="\t", row.names=1))
data=merge(ind, data, by="row.names", all.x=TRUE)

## Convert data
snps.mat <- as.matrix(data[,5:ncol(data)])

## Retaining only polymorphic loci 
lst <- which(lapply(apply(snps.mat[,],2,function(x){x=x[x!=9] ; unique(x)}),length)!=1)

## Function to compute Derived Allele Frequency (DAF) wih regards to missing data for diploid data
daf = function(x){
  x = x[x != 9]
  f = sum(x)/(2*length(x))
}

## Function to compute DAF wih regards to mean missing data for diploid data
sfs.daf = function(x,N=2*nrow(x)){ y = apply(x,2,daf);
            sapply(y,function(y){round(as.integer(N)*y)})}

## Calculating group-specific SFS
pop.size <- by(data$Species,data$Species,function(x){as.integer(length(x))})
write.table(rbind(pop.size),"pop.size_CORRECTED.txt")
sfs.bypop <- by(snps.mat[,lst],data$Species,FUN = sfs.daf)
id <- names(sfs.bypop)
save(sfs.bypop,file="sfs.bypop_Rformat_CORRECTED")


## Transforming in a matrix format
sfs.bypop.b <- NULL
for(i in 1:length(sfs.bypop)){ sfs.bypop.b <- cbind(sfs.bypop.b,sfs.bypop[[i]])}
colnames(sfs.bypop.b) <- id
sfs.bypop.b=as.data.frame(sfs.bypop.b)

## Saving sfs by pop information
write.table(sfs.bypop.b,"sfs.bypop.txt")



