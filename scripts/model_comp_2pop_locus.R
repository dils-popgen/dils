#!/shared/software/miniconda/envs/r-3.5.1/bin/Rscript
# #!/shared/home/croux/.conda/envs/R_env/bin/Rscript
# #!/usr/bin/Rscript

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

library('abcrf')
library('viridis')
library('tidyverse')
# model_comp_2pop.R nameA=txn nameB=ama nSubdir=20  ntree=1000 
for(i in commandArgs()){
	tmp = strsplit(i, '=')
	if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
	if(tmp[[1]][1] == 'ncores'){ ncores = as.integer(tmp[[1]][2]) } # number of cores for the random forest
	if(tmp[[1]][1] == 'ntree'){ ntree = as.integer(tmp[[1]][2]) }
}
# rejected statistics
rejected_statistics = c('minDivAB_avg', 'minDivAB_std', 'maxDivAB_avg', 'maxDivAB_std', 'Gmin_avg', 'Gmin_std', 'Gmax_avg', 'Gmax_std') 

# colors
coul = c('#ffffcc', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84')
coul = colorRampPalette(coul)

# observed data
obs_ss = read.table(paste(timeStamp, '/ABCstat_global.txt', sep=''), h=T) # global statistics over all loci (avg and std)
stats_obs = c(3:15, 20:24)

obs_loci = read.table(paste(timeStamp, '/ABCstat_loci.txt', sep=''), h=T) # individual statistics for each locus

best_model = read.table(paste(timeStamp, '/modelComp/hierarchical_models.txt', sep=''), h=T, sep='\t')
model_demographic = best_model$migration.versus.isolation[1]
submodel_demographic = best_model$IM.versus.SC[1]

### LOCUS SPECIFIC MODEL COMPARISON
outfile = paste(timeStamp, '/locus_modelComp/locus_specific_modelComp.txt', sep='')
if(model_demographic == 'migration'){
	model_genomic = best_model$M.homo.versus.M.hetero[1]
	if(model_genomic == 'Mhetero'){
		mig_ss = read.table(paste(timeStamp, '/locus_modelComp/migration/ABCstat.txt', sep=''), h=T)
		iso_ss = read.table(paste(timeStamp, '/locus_modelComp/isolation/ABCstat.txt', sep=''), h=T)
		
		# locus_specific model comparison
		modIndexes = c(rep('migration', nrow(mig_ss)), rep('isolation', nrow(iso_ss)))

		data_ss = bind_rows(obs_ss, mig_ss, iso_ss)
		
		# remove unused statistics
		toRemove = c(1)
		for(i in 1:ncol(data_ss)){
			if( anyNA(data_ss[,i], recursive=T) || sd(data_ss[-1,i], na.rm=T)<1e-5 || colnames(data_ss)[i]%in%rejected_statistics  ){
				toRemove=c(toRemove, i)
			}
		}
		std_stats = grep('std', colnames(data_ss))
		if( length(std_stats)>0 ){
			toRemove = c(toRemove, std_stats)
		}
		pearson_stats = grep('pearson', colnames(data_ss))
		if( length(pearson_stats)>0 ){
			toRemove = c(toRemove, pearson_stats)
		}
		toRemove = unique(toRemove)

		# keep the same observed stats that the ones retained on simulations
		obs_stat=NULL
		for(stat in colnames(data_ss)[-toRemove]){
			obs_stat = c(obs_stat, which(colnames(obs_loci)==stat))
		}

		mod_iso_mig = abcrf(modIndexes~., data = data.frame(modIndexes, data_ss[-1, -toRemove]), ntree = ntree, paral = T, ncores = ncores)
		#predicted_model_iso_mig = predict(mod_iso_mig, data.frame(data_ss[1, -toRemove]), training=data.frame(modIndexes, data_ss[-1, -toRemove]), ntree = ntree, paral = T, ncores = ncores)
		predicted_model_iso_mig = predict(mod_iso_mig, data.frame(obs_loci[, obs_stat]), training=data.frame(modIndexes, data_ss[-1, -toRemove]), ntree = ntree, paral = T, ncores = ncores)

		allocation = predicted_model_iso_mig$allocation
		post_proba = predicted_model_iso_mig$post.prob
		res = data.frame(obs_loci[, c(1, stats_obs)], allocation, post_proba)
		
		write.table(res, outfile, col.names=T, row.names=F, quote=F, sep='\t', append=F)
	}else{ # if migration but homogeneous
		allocation = rep('migration', nrow(obs_loci))
		post_proba = rep('1', nrow(obs_loci))
		res = data.frame(obs_loci[, c(1, stats_obs)], allocation, post_proba)
		write.table(res, outfile, col.names=T, row.names=F, quote=F, sep='\t', append=F)
	}
}else{
	# change directory
	allocation = rep('isolation', nrow(obs_loci))
	post_proba = rep('1', nrow(obs_loci))
	res = data.frame(obs_loci[, c(1, stats_obs)], allocation, post_proba)
	write.table(res, outfile, col.names=T, row.names=F, quote=F, sep='\t', append=F)
}

