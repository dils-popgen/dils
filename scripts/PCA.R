#!/shared/software/miniconda/envs/r-3.5.1/bin/Rscript

#################################################################################################################################
#################################################################################################################################
#####                                                                                                                       #####
#####    This file is part of Demographic Inferences with Linked Selection : DILS.                                          #####
#####                                                                                                                       #####   
#####    DILS is free software: you can redistribute it and/or modify                                                       #####
#####    it under the terms of the GNU General Public License as published by                                               #####
#####    the Free Software Foundation, either version 3 of the License, or                                                  #####
#####    (at your option) any later version.                                                                                #####
#####                                                                                                                       #####    
#####    DILS is distributed in the hope that it will be useful,                                                            #####
#####    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                     #####
#####    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                      #####
#####    GNU General Public License for more details.                                                                       #####
#####                                                                                                                       #####    
#####    You should have received a copy of the GNU General Public License                                                  #####
#####    along with DILS.  If not, see <https://www.gnu.org/licenses/>.                                                     #####
#####                                                                                                                       #####    
#####    Please send bugreports with examples or suggestions to                                                             #####
#####    camille.roux@univ-lille.fr                                                                                         #####
#####                                                                                                                       #####    
#####    Or write a post on https://groups.google.com/forum/#!forum/dils---demographic-inferences-with-linked-selection     #####
#####                                                                                                                       #####
#################################################################################################################################
#################################################################################################################################

library(FactoMineR)
for(i in commandArgs()){
	tmp = strsplit(i, '=')
	if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
}


input = read.table( paste( timeStamp, '/distribution_PCA.txt', sep=''), h=T, sep='\t')

toRemove = c(1)
for(i in 1:(ncol(input)-1)){
	if(sd(input[,i]) < 0.00001){
		toRemove = c(toRemove, i)
	}
}

toRemove = unique(toRemove)
input = input[, -toRemove]

res.pca <- PCA(input[, -ncol(input)], graph = FALSE, ncp=3)

output_coord = data.frame(res.pca$ind$coord, origin = input$origin)
output_contrib = data.frame(res.pca$var$contrib)

# write the output file
write.table(output_coord, paste( timeStamp, '/table_coord_PCA_SS.txt', sep=''), col.names=T, row.names=F, quote=F, sep='\t') # coordinates to plot the 3D PCA
write.table(output_contrib, paste( timeStamp, '/table_contrib_PCA_SS.txt', sep=''), col.names=T, row.names=T, quote=F, sep='\t') # coordinates to plot the table of contributions
write.table(res.pca$eig, paste( timeStamp, '/table_eigenvalues_PCA_SS.txt', sep=''), col.names=T, row.names=T, quote=F, sep='\t') # coordinates to plot the table of contributions

