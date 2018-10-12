###################################################################################################################################
#
# Copyright 2017 IRD and Grenoble-Alpes University
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
# Written by Philippe Cubry, Yves Vigouroux and Olivier Francois, modified by Nora Scarcelli
#
###################################################################################################################################

###################################################################################################################################
# R script to compute set of frequencies and summary statistics on simulated datasets obtained by Splatche analysis  
# Needs files produced by splatche_launcher.R   
###################################################################################################################################

## Gather arguments from script call ##
args <- commandArgs(TRUE)


## Define useful functions to calculate summarized statistics ##
Nbsample = 160 # Number of samples

## Function to construct SFS ##
myprocess1 = function(
  input= string,
  Nbsample = 160, # Number of samples
  Nbloci = 100,
  Nbindiv = rep(2, 160)) # Number of samples
{
  ## process DNA file
  X <- scan(file = input, what = character(), sep = "\t", quiet = TRUE, skip = 0, nlines = 0, comment.char = "#")
  Totindiv <- sum(Nbindiv)
  lx <- length(X)
  X <- X[30:lx]
  M <- data.frame(matrix(NA, nrow = Totindiv, ncol = (2 + Nbloci)))
  for(k in 1:Nbsample){
    shift <- 1
    if (k == 1) k.s <- 0 else k.s <- sum(Nbindiv[1:(k-1)])
    for (i in 1:(Nbindiv[k])){
      X[shift + (i-1)*(2+Nbloci)] <- k
      M[k.s + i, ] <- X[ ( shift + (i-1)*(2+Nbloci) ):(shift -1 +i*(2+Nbloci) )]
      shift <- shift + 1
    }
    lx <- length(X)
    X<-X[(10 + shift+i*(2+Nbloci)):lx]
  }
  return(M[,-2])
}

## Function to calculate minimal allele frequency ##
maf = function(x){ f = sum(x) ; min(f , length(x) - f ) }

## Function to bin SFS ##
histo.bin = function(x){
  n = length(x)
  y = sapply(1:80, FUN = function(i) sum(x == i) ) #### MUST BE IDENTICAL to function in sum_stat_obs.R
  c(y[1], y[2], sum(y[3:4]), sum(y[5:11]), sum(y[12:21]), sum(y[22:35]), sum(y[36:55]), sum(y[56:80]) )/n #### MUST BE IDENTICAL to function in sum_stat_obs.R
}

## Function to compute mean per group ##
z.groupm = function(z){
  sapply(1:11, FUN = function(i) mean(z[group == i]) ) #### MUST BE IDENTICAL to function in sum_stat_obs.R
}

## Load group definition ##
group = as.numeric(as.vector(read.table("./groups_sim.txt")[,2]))

## Variables definition ##
sum.stat = NULL
param = NULL
spect.stat=NULL

## Loading simulations parameters ##
print("loading parameters")
job_id <- gsub(".tar.gz","",gsub(paste(args[1],"_",sep=""),"",list.files(pattern = args[1])))
for (i in 1:length(job_id)){
  param.simu <- read.table(file = paste("param_",job_id[i],".txt",sep=""),sep="",header=T,colClasses="character")
  param.simu$job_id <- job_id[i]
  param = rbind(param,param.simu)
}

## Loop for creating summarized statistics from simulated datasets ##
for (n in 1:nrow(param)){
  print(paste("loop",n, "of",nrow(param)))
    ## creating a snp-containing object ##
  print("extracting files from archive")
  untar(tarfile=paste(args[1],"_",param["job_id"][n,],".tar.gz",sep=""),
        files = c(paste(args[1],"_",param["job_id"][n,],"/GeneticsOutput/settings",param$generations[n],
                      "_",param$accroissement[n],"_",param$migration[n],"_",param$SizeBeforeExpansion[n],
                      "_",param$TimeOfBottleneck[n],"_",param$AncestralSize[n],"_",param$TimeForRecentBott[n],
                      "_",param$MutationRate[n],"_",param$MainCarryingCapacity[n],"_GeneSamples_",seq(1,10,1),".arp" ,sep=""),
                  paste(args[1],"_",param["job_id"][n,],"/Arrival_cell_",param$generations[n],
                              "_",param$accroissement[n],"_",param$migration[n],"_output.txt" ,sep="")),
        list = FALSE, exdir = ".",compressed = "gzip", verbose = FALSE)


  string.snp = paste("./",args[1],"_",param["job_id"][n,],"/GeneticsOutput/settings",
                     param$generations[n],"_",param$accroissement[n],"_",param$migration[n],"_",
                     param$SizeBeforeExpansion[n],"_",param$TimeOfBottleneck[n],"_",
                     param$AncestralSize[n],"_",param$TimeForRecentBott[n],
                     "_",param$MutationRate[n],"_",param$MainCarryingCapacity[n],"_GeneSamples_1.arp" ,sep="")
  genotype = myprocess1(string.snp, Nbsample = 160, Nbloci = 500, Nbindiv = rep(2, 160) )[,-1] # Change number of samples (twice)


  arrival = (read.table(paste(args[1],"_",param["job_id"][n,],"/Arrival_cell_",param$generations[n],
                  "_",param$accroissement[n],"_",param$migration[n],"_output.txt" ,sep=""),skip=1, sep=":"))
  row.names(arrival) = arrival[,1]; arrival[,1]<-NULL; arrival <- t(arrival)

  for (r in 2:10){
    string1 = paste("./",args[1],"_",param["job_id"][n,],"/GeneticsOutput/settings",param$generations[n],
                    "_",param$accroissement[n],"_",param$migration[n],"_",param$SizeBeforeExpansion[n],
                    "_",param$TimeOfBottleneck[n],"_",param$AncestralSize[n],"_",param$TimeForRecentBott[n],
                    "_",param$MutationRate[n],"_",param$MainCarryingCapacity[n],"_GeneSamples_",r,".arp" ,sep="")
    genotype = cbind(genotype, myprocess1(string1, Nbsample = 160, Nbloci = 500, Nbindiv = rep(2, 160) )[,-1]) # Change number of samples (twice)


  }
  genotype1 = genotype[ (1:320)%%2 == 1, ] # Number of samples x 2
  genotype2 = genotype[ (1:320)%%2 == 0, ] # Number of samples x 2
  genotype = cbind(genotype1, genotype2)

  print("removing temp archive directory")
  unlink(x = paste("./",args[1],"_",param["job_id"][n,],sep=""),recursive = TRUE)

  print("computing summary statistics")


#### The following must be identical to sum_stat_obs.R

  # Load group definition from previous analysis

  group = (read.table("./groups_obs.txt", na.strings="-9"))
  group = as.numeric(as.vector(group[,2]))

  spect = apply(genotype, MARGIN = 2, FUN = function(x) {maf(as.numeric(as.factor(x)) - 1)}   )

  lst1 = which(spect == 1)
  geno1 = as.data.frame(genotype[,lst1])
  ind1 = apply(geno1, MARGIN = 2, FUN = function(x) {
    y = as.numeric(as.factor(x))
    i = which( y != y[1] )
    if (length(i) > 1) 1 else i
  })


  z = sapply(1:160, FUN = function(x) sum(ind1 == x) ) # Change number of samples

  stat.sim <- c((histo.bin(spect)/sum(histo.bin(spect))), z.groupm(z)/sum(z.groupm(z)))
  names(stat.sim) <- c(paste("SFS",seq(1,8,1), sep=""),paste("RareVariants",seq(1,11,1), sep=""))

  stat.temp <- as.data.frame(c(param[n,],stat.sim,as.data.frame(arrival)))
  spect.temp <- as.data.frame(c(param[n,],spect))

  sum.stat <- rbind(sum.stat, stat.temp )
  spect.stat <- rbind(spect.stat,spect.temp)
  cond <- c(rep("SFS",8),rep("RareVariants",11)) 

}
print("saving computed statistics to file")
write.table(file = paste(args[2],".splatche_sum.stat.txt",sep=""), sum.stat, row.names = F, quote=F)
