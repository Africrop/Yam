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
# Written by Philippe Cubry and Olivier Francois, modified by Nora Scarcelli
#
###################################################################################################################################

###################################################################################################################################
# This script launch Splatche2
# Usage : splatche_launcher.R <number of runs> <number of source populations>
# Needs splatche_input_files_creation.R
# Needs all the files within the folder "dataset_1_layer":
# Arrival_cell.col
# dens_init.txt
# dynamic_F.txt
# dynamic_K.txt
# GeneSamples.sam
# genetic_data_SEQ.par
# ppveg.asc
# rivers.asc
# roughness.asc
# settings.txt
# veg2F_pop1_time_1.txt
# veg2F_pop1_time_2.txt
# veg2K_pop1_time_1.txt
# veg2K_pop1_time_2.txt
###################################################################################################################################


#### getting Job_Id ####
job <- paste(Sys.getenv("JOB_ID"),"-",Sys.getenv("SGE_TASK_ID"),sep="")

#### getting call arguments (number of simulation = $1, number of source populations = $2) ####
args <- commandArgs(TRUE)

#### source useful file ####
source("./splatche_input_files_creation.R")

#### Definition of number of runs and source populations ($1 and $2 ou 1 and 1 by default) ####
if (is.na(args[1])==TRUE){
  nrun <- 1
  norigin <- 1
} else {nrun<-as.numeric(args[1]); norigin<-as.numeric(args[2])}


for (it in 1:nrun){
  print(paste(it,"over",nrun,sep=" "))

##### Random draw of coordinates #####
  coordo <- NULL

  for (ori in 1:norigin){
    ox = runif(1, -16, 40) # Longitude
    oy = runif(1, -5, 20) # Latitude

    while(!is.african(cbind(ox,oy))){
      ox = runif(1, -16, 40) # Longitude
      oy = runif(1, -5, 20) # Latitude
    }
    coordo = rbind(coordo,c(ox,oy))
  }

  ##### Random draw of variable #####
  ng <- sample(gen_min:gen_max,1) # Number of generation
  rate <- runif(1,acc_min,acc_max) # Growth rate
  mig <- runif(1,mig_min,mig_max) # Migration rate
  res <- sample(res_min:res_max,1) # Size during bottleneck
  Tres <- sample(Tres_min:Tres_max,1) # Bottleneck lenght
  resB <- sample(resB_min:resB_max,1) # Size before bottleneck
  MutRate<- sample(MutRate_prior,1) # Mutation rate
  K <- sample(K_min:K_max,1) # Carrying capacity

  TrecBott <- NULL
  TrecBott <- sample(TrecBott_min:ng,1)
  while(!TrecBott<= ng){
  TrecBott <- sample(TrecBott_min:TrecBott_max,1)
  }

  splatche(input = "./datasets_1layer/settings.txt", ng, rate, mig, norigin,res,Tres,resB,TrecBott, MutRate, K, coord.o = coordo) 

  
#### Launching Splatche2
 system(paste("splatche2-01 settings",ng,"_",rate,"_",mig,"_",res,"_",Tres,"_",resB,"_",TrecBott,"_",MutRate,"_",K, ".txt", sep=""))

 
##### Saving output files
  if(file.exists(paste("./datasets_1layer/GeneticsOutput/settings",ng,"_",rate,"_",mig,"_",res,"_",
                       Tres,"_",resB,"_",TrecBott,"_",MutRate,"_",K,"_GeneSamples_2.arp" ,sep=""))){

    coo1 <- NULL ; noms <- NULL
    for (npop in 1:norigin){
      coo1<-cbind(coo1,coordo[npop,1],coordo[npop,2])
      noms <-cbind(noms,paste("Long",npop,sep=""),paste("Lat",npop,sep=""))
    }
    param<- rbind(param,c(coo1,ng,rate,mig,res,Tres,resB,TrecBott,MutRate,K))
    colnames(param) <- c(noms,"generations","accroissement","migration","SizeBeforeExpansion",
                         "TimeOfBottleneck","AncestralSize","TimeForRecentBott","MutationRate","MainCarryingCapacity")
  }
  write.table(file = paste("param_",job,".txt",sep=""), param, row.names = F, quote=F)
}
