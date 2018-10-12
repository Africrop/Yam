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
# R script to convert vcf file into 012 format  
###################################################################################################################################


# Set the input vcf file
input <- "In.vcf"
#Define the output012 file
output <- "Out.012"
#### Open the connection and determine the length (in row) of the header ####
print("open the connection")
con_in <- file(input, open='r')
i <- 1
test <- readLines(con_in,n=1)
while(grepl(pattern="#", test)==TRUE){
  test <- readLines(con_in,n=1)
  i <- i+1
} 
close(con_in)

#### Open a new connection to read only the last line of the header ####
con_in <- file(input, open='r')
i <- i-1
header <- readLines(con_in,n=i)
close(con_in)
head <- unlist(strsplit(header[length(header)], split="\t"))

#### Open a new connection to read only the last line of the header ####
con_in <- file(input, open='r')
readLines(con_in,n=i)

#### Saving the resulting dataframe in a new 012 file ####
print("saving")
write.table(t(c("ID",head[10:length(head)])),output,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

j <- 0
##### now we read the file line by line and extract SNPs that can be polarized ####
while (length(oneLine <- readLines(con_in, n = 1, warn = FALSE)) > 0) {
  ###### Increment the counter and print it ######
  j = j + 1
  print(paste("Reading line ",j,sep=""))
  
  ###### splitting the raw data ######
  oneLine<- unlist(strsplit(oneLine, split="\t"))
  names(oneLine) <- head
  ID <- paste(oneLine["#CHROM"],oneLine["POS"],sep="_"); names(ID) <- "ID"
  oneLine <- oneLine[10:length(oneLine)]  
  oneLine[oneLine == "0/0"] <- "0"
    oneLine[oneLine == "1/1"] <- "2"
    oneLine[oneLine == "0/1"] <- "1"
    oneLine[oneLine == "1/0"] <- "1"
    oneLine[oneLine == "./."] <- "9"
    oneLine <- c(ID,oneLine)
    write.table(t(oneLine),output,append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
  } 

close(con_in)


