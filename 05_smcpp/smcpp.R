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
# Written by Philippe Cubry
#
###################################################################################################################################

####################################
# R script to represent smc++ results 
####################################


smc <- read.csv("output_unfold.csv")
smc.notk <- read.csv("output_unfold_notk.csv")

library(ggplot2)
library(cowplot)

smc_plot <- ggplot() +
  geom_line(aes(x,y,color="cult"),data=smc)+
  coord_cartesian(ylim=c(2000,150000),xlim=c(300,1000000))+
  scale_x_log10(breaks=c(250,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000),
                labels=c(200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,"100000","200000","300000","400000","500000","600000","700000","800000","900000","1000000"))+
  scale_y_log10(breaks=c(2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000),
                labels=c(2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,"1000000"))+
  xlab("Time in generations")+
  ylab("Population effective size")+theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1,size=7),
                                          axis.text.y = element_text(size=7))


save_plot(smc_plot,
          filename = "results.pdf",
          base_aspect_ratio = 2,
          base_height = 5)

smc_plot.notk <- ggplot() +
  geom_line(aes(x,y,color="cult"),data=smc.notk)+
  coord_cartesian(ylim=c(2000,150000),xlim=c(300,1000000))+
  scale_x_log10(breaks=c(250,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000),
                labels=c(200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,"100000","200000","300000","400000","500000","600000","700000","800000","900000","1000000"))+
  scale_y_log10(breaks=c(2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000),
                labels=c(2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,"1000000"))+
  xlab("Time in generations")+
  ylab("Population effective size")+theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1,size=7),
                                          axis.text.y = element_text(size=7))

