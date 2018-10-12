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
# This script create input files for Splatche2
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

#### Definition of prior distribution ####
gen_min <- 50 ; gen_max <- 2000 # min and max generations number
acc_min <- 0.1 ; acc_max <- 0.6 # min and max growth rate
mig_min <- 0.15 ; mig_max <- 0.7 # min and max migration rate
res_min <- 100 ; res_max <- 10000 # min and max size during bottleneck
resB_min <- 500 ; resB_max <- 100000 # min and max size before bottleneck
Tres_min <- 0 ; Tres_max <- 5000 # min and max duration of bottleneck parameter
TrecBott_min <- 1 ; TrecBott_max <- 2000 # min and max time in generations for the recent bottleneck
MutRate_prior <- 0.000001 # mutation rate
K_min <- 30 ; K_max <- 150 # min and max carrying capacity for each cluster

param <- NULL


#### Creating input files for Splatche ####

splatche = function(input = "settings.txt",
                    ng = 1200,
                    rate = 0.1,
                    mig = 0.03,
                    norigin = 1,
                    res = 0,
                    Tres = 0,
                    resB=0,
                    TrecBott=0,
                    MutRate= 0.000001,
                    K = 100,
                    coord.o = coordo) {
 
 XX = scan(file = input,
            what = character(),
            sep = "\t",
            quiet = TRUE,
            skip = 0,
            nlines = 0,
            comment.char = "#")
  XX[1] = paste("PopDensityFile=./datasets_1layer/dens_init_",job,".txt",sep = "")
  XX[8] = paste("EndTime=",ng,sep = "")
  XX[10] = paste("GrowthRate=",rate,sep = "")
  XX[11] = paste("MigrationRate=",mig,sep = "")
  XX[15] = paste("ArrivalCellFile=./datasets_1layer/Arrival_cell_",ng,"_",rate,"_",mig,".col",sep="")

  YY = NULL
  YY[1] = as.character(norigin)
  for (i in 1:norigin) {
    YY[i+1] = paste("source",i," ",  100," ",	coord.o[i,2]," ",	coord.o[i,1]," ",	res," ",	Tres," ",	resB," ",	0," ",	0," ",	0, sep = "")
  }
  YY[norigin + 2] = "Pop_name pop_size Lat Long ResizeA Time_ResizeB ResizeB MigrationFromCurrentSource NoLayer TimeOfExpansion"

  ZZ = NULL
  ZZ[1] = 2
  ZZ[2] = "1 ./datasets_1layer/veg2K_pop1_time_1.txt Time1"
  ZZ[3] = paste(TrecBott," ./datasets_1layer/veg2K_pop1_time_2.txt reducing by half the carrying capacity",sep="")

  AA = scan(file = "./datasets_1layer/genetic_data_SEQ.par",
            what = character(),
            sep = "\t",
            quiet = TRUE,
            skip = 0,
            nlines = 0,
            comment.char = "")
  AA[seq(4,3002,3)] = paste("DNA    1   0    ",MutRate,"   1",sep = "")

  BB = scan(file = "./datasets_1layer/veg2K_pop1_time_1.txt",
            what = character(),
            sep="\n",
            quiet = TRUE,
            skip = 0,
            nlines = 0,
            comment.char = "")
  BB[3] <- paste("3\t",K,"\tTropical woodland",sep="")
  BB[4] <- paste("4\t",K,"\tTropical thorn scrub and scrub woodland",sep="")
  BB[5] <- paste("5\t",K,"\tTropical semi-desert",sep="")
  BB[6] <- paste("6\t",K,"\tTropical grassland",sep="")
  BB[8] <- paste("8\t",K,"\tSavanna",sep="")
  BB[9] <- paste("9\t",K,"\tBroadleaved temperate evergreen forest",sep="")
  BB[12] <- paste("12\t",K,"\tTemperate deciduous broadleaved forest",sep="")

  CC = scan(file = "./datasets_1layer/veg2K_pop1_time_2.txt",
            what = character(),
            sep="\n",
            quiet = TRUE,
            skip = 0,
            nlines = 0,
            comment.char = "")
  CC[3] <- paste("3\t",K,"\tTropical woodland",sep="")
  CC[4] <- paste("4\t",K,"\tTropical thorn scrub and scrub woodland",sep="")
  CC[5] <- paste("5\t",K,"\tTropical semi-desert",sep="")
  CC[6] <- paste("6\t",K,"\tTropical grassland",sep="")
  CC[8] <- paste("8\t",K,"\tSavanna",sep="")
  CC[9] <- paste("9\t",K,"\tBroadleaved temperate evergreen forest",sep="")
  CC[12] <- paste("12\t",K,"\tTemperate deciduous broadleaved forest",sep="")


  AR = scan(file="./datasets_1layer/Arrival_cell.col", what = character(),  quiet = TRUE, skip = 0, nlines = 0, comment.char = "#",sep="\n")

  write(file = paste("settings",ng,"_",rate,"_",mig,"_",res,"_",Tres,"_",resB,"_",TrecBott,"_",MutRate,"_",K, ".txt", sep = ""),
        XX,
        sep = "\t")
  write(file = paste("./datasets_1layer/dens_init_",job,".txt",sep=""),
        YY,
        sep = "\t")
  write(file = paste("./datasets_1layer/dynamic_K.txt",sep=""),
        ZZ,
        sep = "\t")
  write(file = paste("./datasets_1layer/genetic_data_SEQ.par",sep=""),
        AA,
        sep = "\t")
  write(file = paste("./datasets_1layer/veg2K_pop1_time_1.txt",sep=""),
        BB,
        sep = "\t")
  write(file = paste("./datasets_1layer/veg2K_pop1_time_2.txt",sep=""),
        CC,
        sep = "\t")
  write(file = paste("./datasets_1layer/Arrival_cell_",ng,"_",rate,"_",mig,".col",sep=""),
        AR,
        sep = "\t")

}

africa = read.table("./datasets_1layer/ppveg.asc", skip = 6)

is.african = function(coord){
  i = 1 + as.integer((coord[,1]  + 16.73322)/0.83)
  j = 87 -  as.integer((coord[,2]  + 35.23718)/0.83)
  return( africa[j,i] >= 0)
}
