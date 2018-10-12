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
# Written by Philippe Cubry, adapted by Roland Akakpo
#
###################################################################################################################################

###################################################################################################################################
# R script to estimate the best scenario from FastSimCoal outputs
###################################################################################################################################

# get the best likelihoods results for each model and each run (100 runs were performed for each model)
results <- NULL
for(m in c("1","2","3")){
  for( i in list.files(paste0("Single",m))){
    if(file.exists(paste("./Single",m,"/",i,"/Single",m,".bestlhoods",sep=""))){
      results <- rbind(results,
                       cbind(m, i,
			     read.table(paste("./Single",m,"/",i,"/Single",m,".bestlhoods",sep=""),
                                             header=TRUE)))
    }
  }
}


# get the max likelihood for each model
ML <- NULL
for(m in levels(results$m)){
  ML <- rbind(ML, data.frame(m, results[which(results$MaxEstLhood==(max(results[results$m == m,]$MaxEstLhood))),]$i,max(results[results$m == m,]$MaxEstLhood)))
}
colnames(ML) <- c("model","run", "max_estimated_likelihood")
print(ML)



## Simulations for likelihood estimation
# We then used the Maximum Likelihood parameters estimates identified for each scenario in the previous section as an input to estimate more precisely the likelihood. We performed 100 runs of 1000000 simulations for each scenario.
# We then identified the best run among those 100 in order to identify the likelihood that will be used to compute AIC for model comparison. We also calculate the Akaikes weight of evidence in favor of each model as suggested by Excoffier et al, 2013.


# get the best likelihoods results for each model and each run (100 runs were performed for each model)
results <- NULL
for(m in c("1","2","3")){
  for( i in list.files(paste0("Single",m))){
    if(file.exists(paste("./Single",m,"/",i,"/Single",m,"_","maxL.lhoods",sep=""))){
      results <- rbind(results,
                       cbind(m, i,
                                  read.table(paste("./Single",m,"/",i,"/Single",m,"_","maxL.lhoods",sep=""),
                                             header=TRUE)))
    }
  }
}

# get the max likelihood for each model
ML <- NULL
for(m in levels(results$m)){
  ML <- rbind(ML, data.frame(m,max(results[results$m == m,]$DAFLHood_1)))
}
colnames(ML) <- c("model","max_likelihood")

# Compute the AIC = 2d - 2ln(lhood)  for each model taking the best run
# Assuming xx parameters for each model (xx=8 for model including 4 pops (cultivated + 3wild) and xx=7 for model including 3 wild pops)
# Caution : convert log10(likelihoods) from fsc to ln(likelihoods) by multiplying by log(10)
ML$AIC <- 2*8 - 2*(ML$max_likelihood*log(10))
# Compute exp(-0.5*DELTAi) with DELTAi = AICi - AICmin
ML$expDeltai <- exp(-0.5*(ML$AIC - min(ML$AIC)))
# Compute wi the Akaike's weight of evidence in favor of ith model
ML$wi <- ML$expDeltai/sum(ML$expDeltai)

print(ML)

