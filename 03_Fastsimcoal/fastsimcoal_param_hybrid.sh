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
# Written by Roland Akakpo, modified by Nora Scarcelli
#
###################################################################################################################################


###################################################################################################################################
# This script simulates the scenario described in Hybrid.est and Hybrid.tlp
# Usage: qsub fastsimcaol_param_hybrid.sh $1
# $1 = name of the model to test (example files: $1=Hybrid)
# Require in the same folder the files $1.obs, $1.est and $1.tpl (example files: Hybrid.obs, Hybrid.est, Hybrid.tlp)
####################################################################################################################################


#!/bin/bash

# check if option --removeZeroSFS
model=$1
minL=40
maxL=150
minSim=250000
maxSim=1000000

fsc25221 -t ${model}.tpl -e ${model}.est -n $minSim -N $maxSim -M 0.001  -l $minL -L $maxL --multiSFS -d -B 1 -C 10


