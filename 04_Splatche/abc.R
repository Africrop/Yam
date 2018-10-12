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
# Written by Philippe Cubry, Yves Vigouroux and Olivier Francois, modified by Nora Scracelli
#
###################################################################################################################################

###################################################################################################################################
# R script to compute abc analysis from Splatche2 output and represent the results on a map
# Needs files produced by Group.R, sum_stat_obs.R et sum_stat_sim.R
###################################################################################################################################

library(abc)
library(ggplot2)
library(raster)
library(rgdal)
library(cowplot)


# Load map elements - TO BE DOWNLOADED FROM NATURALEARTH website : http://www.naturalearthdata.com/
nat.earth <- stack("./naturalearth/HYP_50M_SR_W.tif")
lakes <- readOGR("./naturalearth",layer = "ne_10m_lakes")
rivers <- readOGR("./naturalearth",layer = "ne_50m_rivers_lake_centerlines")
playas <- readOGR("./naturalearth",layer = "ne_50m_playas")
ocean <- readOGR("./naturalearth",layer = "ne_10m_ocean")
shapefile <- readOGR("./naturalearth",layer = "ne_50m_admin_0_countries")
data <- fortify(shapefile)
data.lakes <- fortify(lakes)
data.playas <- fortify(playas)
data.ocean <- fortify(ocean)
data.rivers <- fortify(rivers)
# Define extent and load raster for the map
extent <- c(-25,57,-15,40) #### Coordinates for Africa
nat.crop <- crop(nat.earth, y=extent(extent))
rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop/255))
rast.table$rgb <- with(rast.table, rgb(HYP_50M_SR_W.1,
                                       HYP_50M_SR_W.2,
                                       HYP_50M_SR_W.3,
                                      1))

rm(nat.earth,nat.crop)

# functions definitions
densities_plot <- function(prior = prior, post = post, title = "") {ggplot()  +
            geom_density(aes(prior,fill="priors"),alpha=0.5) +
            geom_density(aes(post,fill="posteriors"),alpha=0.5) +
            scale_fill_manual(values=c('#E69F00','#999999')) +
            labs(title= title,x="") +
            theme(panel.grid=element_blank(),legend.title=element_blank())
}

#### Projection of longitude and latitude posteriors on a map
posteriorposition_map <- function(posteriorposition = posteriorposition,title="",n=50,h=25) {
  ggplot()+  
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill="lightblue", data = data.ocean, size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data,
                 colour ='antiquewhite4', fill = "white", size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data.lakes,
                 colour ='lightblue',fill="lightblue", size = .3) +
    stat_density2d(data = posteriorposition, 
                   aes(x=Long1,y=Lat1,  fill = ..level.., alpha = ..level..),
                   size = 0.1, geom = 'polygon',n = n,h=h) +
     geom_path(aes(x = long, y = lat, group = group), data = data.rivers,
                 colour ='lightblue') +
    scale_fill_gradient(low = "yellow", high = "red", guide = FALSE) +
    scale_alpha(range = c(0, 0.25), guide = FALSE) +
    ggtitle(title) +
    theme_bw() +
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "black",
                                      size = .5,
                                      linetype = "solid"
          )
    ) + coord_quickmap(xlim = c(-20,50),ylim =c(-10,35))
}



estimate_mode <- function(x) { 
  d <- density(x)
  d$x[which.max(d$y)]
}



#### Reading output created by sum_stat_obs.R and sum_stat_sim.R
list_input = list.files(path = ".", pattern="Splatche_sum.stat.txt")
param.moinsMR <- NULL
for (f in list_input) {
  dat <- read.table(f, header=T)
  param.moinsMR <- rbind(param.moinsMR, dat)
}

stat.obs_moinsMR = read.table("./sum_stats_obs.txt",header = T)

# boxplot of simulated summary statistics with comparison with observed summary statistics
#### Adjust SFS bin number, number of groups for rare variant calculation and corresponding column numbers for statistics and parameters
boxplot(param.moinsMR[,c(13:31)], #### Statistics column number (SFS + RV) from file XXX.Splatche_sum.stat.txt
        las=2, main = "Simulated vs Observed sum stat",xaxt="n", ylim=c(0,1))
axis(1, at= 1:19,c(paste("SFS",seq(1:8)),paste("RV",seq(1:11))), las=2) #### [at= 1:xx] = number of parameters; ["SFS",seq(1:xx)] = number of SFS bins; ["RV",seq(1:xx)] = number of groups for rare variant calculation
points(1:19,stat.obs_moinsMR,pch=19,col="red") #### [1:xx] = parameters column number from file sum.stats_obs.txt
legend("topright","Observed",pch=19,col="red")


# ABC analysis
#### Adjust SFS bin number, number of groups for rare variant calculation and corresponding column numbers for statistics and parameters
nn10pop <- abc(target = stat.obs_moinsMR[c(1:19)],sumstat = param.moinsMR[,c(13:31)],param = param.moinsMR[,c(1:8,11)],method = "neuralnet",tol = 0.01,numnet = 500) 

temp=nn10pop$unadj.values

# Plotting posteriors of the parameters
#### param.moinsMR$Long1 = prior
#### nn10pop$adj.values = posterior
long <- densities_plot(param.moinsMR$Long1,nn10pop$adj.values[,"Long1"],title = "Longitude")
lat <- densities_plot(param.moinsMR$Lat1,nn10pop$adj.values[,"Lat1"],title = "Latitude")
gen <- densities_plot(param.moinsMR$generations,nn10pop$adj.values[,"generations"],title = "n.o. generations")
acc <- densities_plot(param.moinsMR$accroissement,nn10pop$adj.values[,"accroissement"],title = "Growth rate")
mig <- densities_plot(param.moinsMR$migration,nn10pop$adj.values[,"migration"],title = "Migration rate")
sizeBeforeExp <- densities_plot(param.moinsMR$SizeBeforeExpansion,nn10pop$adj.values[,"SizeBeforeExpansion"],title = "Size before expansion")
TimeBott <- densities_plot(param.moinsMR$TimeOfBottleneck,nn10pop$adj.values[,"TimeOfBottleneck"],title = "Bottleneck length")
TimeRecBott <- densities_plot(param.moinsMR$TimeForRecentBott,nn10pop$adj.values[,"TimeForRecentBott"],title = "Rec_bott_time")
AncSize <- densities_plot(param.moinsMR$AncestralSize,nn10pop$adj.values[,"AncestralSize"],title = "Ancestral size")
K <- densities_plot(param.moinsMR$MainCarryingCapacity,nn10pop$adj.values[,"MainCarryingCapacity"],title = "K")

prow<-plot_grid(lat + theme(legend.position="none"),
                long + theme(legend.position="none"),
                gen + theme(legend.position="none"),
                mig + theme(legend.position="none"),
                acc + theme(legend.position="none"),
                sizeBeforeExp + theme(legend.position="none"),
                TimeBott + theme(legend.position="none"),
                TimeRecBott + theme(legend.position="none"),
                AncSize + theme(legend.position="none"),
                K + theme(legend.position="none")
                )
grobs <- ggplotGrob(lat)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid( prow, legend, rel_widths = c(3, .5))
dev.off()
save_plot(filename = "./posterior.pdf",plot_grid( prow, legend, rel_widths = c(3, .5)),base_height = 8,base_aspect_ratio = 2)
plot(nn10pop, param = param.moinsMR[,c(1:8,11)])

# Goodness of fit
param.moinsMR_NA = na.omit(param.moinsMR[,c(13:31)])
gfitpca(target = stat.obs_moinsMR[c(1:19)], sumstat = param.moinsMR_NA, cprob = 0.1, index=rep("90%",nrow(param.moinsMR_NA)))


# Plotting posterior on the map
map <- as.data.frame(nn10pop$adj.values[,1:2])
maps <- plot_grid(posteriorposition_map(map,n=100,h=NULL))
maps
save_plot(filename = "./map.pdf",maps)
