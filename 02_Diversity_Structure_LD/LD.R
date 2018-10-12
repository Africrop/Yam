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
# Written by Philippe Cubry, based on Benedicte Rhone, modified by Nora Scarcelli
#
###################################################################################################################################


###################################################################################################################################
##This script calculate mean r2 values for each position from vcftools generated files
# LD is then plotted 
####################################################################################################################################

library(data.table)
library(ggplot2)
library(cowplot)

# Initiate variable
DL <- NULL

# Read the vcftools produced table of r2 values
tmp <- fread("out.geno.ld")
names(tmp)=c("CHR","POS1","POS2","N_INDIV","R2")

# Subsampling of the matrix
tmp2=tmp[sample(nrow(tmp),size=500000,replace=FALSE),]

# Calculate distance between considered polymorphisms
tmp2$dist <- tmp2$POS2 - tmp2$POS1
DL <- cbind(tmp2$R2,tmp2$dist)
colnames(DL) <- c("R2","dist")

# Calculate mean r2 for each considered distance
DL <- as.data.frame(DL)
meanR2 <- NULL
meanR2 <- cbind(as.data.frame(tapply(DL$R2, DL$dist,mean)), as.data.frame(tapply(DL$R2, DL$dist,var)) ,as.data.frame(table(DL$dist)))
names(meanR2)=c("Mean","Var","BP","Freq")
# Save the resulting file
write.table(file="mean_R2.txt",x= meanR2)

##### Computation of Hill and Weir expectations - script adapted from Benedicte Rhone #####
n.sample<-100 # number of sample

dl <- read.table("mean_R2.txt",header=TRUE)
dl=dl[!is.na(dl$Mean), ]
distance<-dl$BP
LD.data<-dl$Mean

HW.st<-c(C=0.2)                   
HW.nonlinear<-nls(LD.data~(((10+C*distance)/
                                         ((2+C*distance)*(11+C*distance)))*
                                         (1+((3+C*distance)*(12+12*C*distance+
                                         (C*distance)^2))/(n.sample*(2+C*distance)*
                                         (11+C*distance)))), 
                        start=HW.st, control=nls.control(maxiter=100),
                        weights = dl$Freq)
HW.nonlinear.prediction <- predict(HW.nonlinear,list(distance = dl$Var1))

summary(HW.nonlinear)
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
new.rho

##### Plotting LD against distance in pb
smoothness <- 30

p <- ggplot() +
  geom_line(aes(dl$BP,HW.nonlinear.prediction,color="grey"),linetype = "solid")+
  stat_smooth(aes(dl$BP,dl$Mean,
                  weight=dl$Freq,color="black"),linetype = "solid",method = "gam",formula = y ~ s(x,bs="cs",k=smoothness),se = FALSE)+
  scale_color_manual(name="",values=c("black","grey"),labels=c("Smoothed mean r2","Hill and Weir (1988) expectation"))+
  labs(x="Distance (in bp)",y= "r2")
ggsave(filename = "LD_Decay.pdf",p)


