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
# model_comp_2pop.R nameA=txn nSubdir=20  ntree=1000 
for(i in commandArgs()){
	tmp = strsplit(i, '=')
	if(tmp[[1]][1] == 'nameA'){ nameA = tmp[[1]][2] }
	if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nMin'){ nMin = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'sub_dir_sim'){ sub_dir_sim = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nSubdir'){ nSubdir = as.integer(tmp[[1]][2]) } # number of subdirectories where simulations were ran
	if(tmp[[1]][1] == 'ncores'){ ncores = as.integer(tmp[[1]][2]) } # number of cores for the random forest
	if(tmp[[1]][1] == 'ntree'){ ntree = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'outgroup'){ outgroup = as.integer(tmp[[1]][2]) } # 0: no outgroup, no SFS used. 1: outgroup, SFS used
}

outfile = paste(timeStamp, '/', sub_dir_sim, '/report_', nameA, '.txt', sep='')
outfile_best = paste(timeStamp, '/', sub_dir_sim, '/best_model.txt', sep='')

# colors
coul = c('#ffffcc', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84')
coul = colorRampPalette(coul)

# observed data
obs_ss = read.table(paste(timeStamp, '/ABCstat_global.txt', sep=''), h=T) # global statistics over all loci (avg and std)

obs_loci = read.table(paste(timeStamp, '/ABCstat_loci.txt', sep=''), h=T) # individual statistics for each locus
obs_sfs = read.table(paste(timeStamp, '/ABCjsfs.txt', sep=''), h=T)

ss_obs = cbind(obs_ss, obs_sfs)

sfs=matrix(as.numeric(obs_sfs), byrow=T, nrow=1)
colnames(sfs) = paste('f_', nameA, '=', 0:(nMin-1), sep='')
sfs = sfs[1, -c(1, length(sfs))]

write.table(matrix(sfs, ncol=length(sfs), dimnames=list(NULL, names(sfs))), paste(timeStamp, '/sfs_table.txt', sep=''), col.names=T, row.names=T, sep='\t', quote=F)


pdf(paste(timeStamp, '/sfs_plot.pdf', sep=''), bg="white")
barplot(sfs, col=rev(viridis(1)), cex=1.25, cex.axis=1.25, cex.lab=1.25)
dev.off()


# simulated data
#models = c('Constant_1N', 'Constant_2N', 'Discrete_1N', 'Discrete_2N', 'Expo_1N', 'Expo_2N')
#growth = c('constant', 'constant', 'changes', 'changes', 'changes', 'changes')
models = c('Constant_1N', 'Constant_2N', 'Expansion_1N', 'Expansion_2N', 'Contraction_1N', 'Contraction_2N')
growth = c('Constant', 'Constant', 'Expansion', 'Expansion', 'Contraction', 'Contraction')
homoN = c('1N', '2N', '1N', '2N', '1N', '2N')
ss_sim = list()
params_sim = list()
all_models_sim = NULL

for(m in models){
	ss_sim_tmp = NULL
	params_sim_tmp = NULL
	for(rep in seq(0, nSubdir-1, 1)){
		# check if the simulation went well
		if(file.exists(paste(timeStamp, '/', sub_dir_sim, '/', m, '_', rep, '/priorfile.txt', sep='')) == TRUE){
			# statistics
			tmp_ss = read.table(paste(timeStamp, '/', sub_dir_sim, '/', m, '_', rep, '/ABCstat.txt', sep=''), h=T)
			
			tmp_sfs = read.table(paste(timeStamp, '/', sub_dir_sim, '/', m, '_', rep, '/ABCjsfs.txt', sep=''), h=T)
			tmp = cbind(tmp_ss, tmp_sfs)
			
			ss_sim_tmp = rbind(ss_sim_tmp, tmp)
			
			# params
			tmp_params = read.table(paste(timeStamp, '/', sub_dir_sim, '/', m, '_', rep, '/priorfile.txt', sep=''), h=T)
			params_sim_tmp = rbind(params_sim_tmp, tmp_params)
		}
	}
	# statistics
	ss_sim[[m]] = ss_sim_tmp 
	all_models_sim = rbind(all_models_sim, ss_sim_tmp)
	
	# params
	params_sim[[m]] = params_sim_tmp
}


# remove uninformative statistics: those with no variation
ss_2_remove = c(1)
for(i in 2:ncol(all_models_sim)){
	if( sd(all_models_sim[, i] )<1e-4 ){
		ss_2_remove = c(ss_2_remove, i)
	}
}

for(m in names(ss_sim)){
	for(i in 2:ncol(ss_sim[[m]])){
		if( sd(ss_sim[[m]][,i] ) <1e-4 ){
			ss_2_remove = c(ss_2_remove, i)
		}
	}
}

ss_2_remove = unique(ss_2_remove)


summary_modelComp = NULL # ex: c('migration_vs_isolation', 'IM_vs_AM', 'Mhomo_vs_Mhetero', 'Nhomo_vs_Nhetero')
summary_bestModel = NULL # ex: c('migration', 'IM', 'Mhetero', 'Nhetero')
summary_proba = NULL # ex: c(0.9, 0.6, 0.8, 0.7)

# model comparison #1 --> two models: isolation versus migration
modIndexes = NULL
for(i in 1:length(models)){
	modIndexes = c(modIndexes, rep(growth[i], nrow(ss_sim[[models[i]]])))
}

mod_cons_grow = abcrf(modIndexes~., data = data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
predicted_model_cons_grow = predict(mod_cons_grow, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)

write('\n#####\n\nMODEL COMPARISON #1: constant versus growth', outfile, append=F)
write('#confusion matrix:', outfile, append=T)
write.table(mod_cons_grow$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
write(paste('\n#best model: ', predicted_model_cons_grow$allocation, sep=''), outfile, append=T)
write(paste('#proba best model: ', predicted_model_cons_grow$post.prob, sep=''), outfile, append=T)
write('\n#votes:', outfile, append=T)
write.table(t(as.matrix(predicted_model_cons_grow$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

summary_modelComp = c(summary_modelComp, 'constant vs expansion vs contraction')
summary_bestModel = c(summary_bestModel, as.character(predicted_model_cons_grow$allocation))
summary_proba = c(summary_proba, predicted_model_cons_grow$post.prob)



# if ongoing migration
# model comparison #2 --> IM versus SC
if(predicted_model_cons_grow$allocation=='Constant'){
	homoN = 'Constant_1N'
	heteroN = 'Constant_2N'
}
if(predicted_model_cons_grow$allocation=='Expansion'){
	homoN = 'Expansion_1N'
	heteroN = 'Expansion_2N'
}
if(predicted_model_cons_grow$allocation=='Contraction'){
	homoN = 'Contraction_1N'
	heteroN = 'Contraction_2N'
}

all_models = NULL
modIndexes = NULL
for(i in homoN){
	all_models = rbind(all_models, ss_sim[[i]])
	modIndexes = c(modIndexes, rep('homoN', nrow(ss_sim[[i]])))
}
for(i in heteroN){
	all_models = rbind(all_models, ss_sim[[i]])
	modIndexes = c(modIndexes, rep('heteroN', nrow(ss_sim[[i]])))
}

mod = abcrf(modIndexes~., data = data.frame(modIndexes, all_models[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
predicted_model = predict(mod, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, all_models[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)

write('\n#####\n\nMODEL COMPARISON #2: Nhomo versus Nhetero', outfile, append=T)
write('#confusion matrix:', outfile, append=T)
write.table(mod$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
write(paste('\n#best model: ', predicted_model$allocation, sep=''), outfile, append=T)
write(paste('#proba best model: ', predicted_model$post.prob, sep=''), outfile, append=T)
write('\n#votes:', outfile, append=T)
write.table(t(as.matrix(predicted_model$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

write(paste(predicted_model$allocation, '\n', sep=''), outfile_best, append=F)

summary_modelComp = c(summary_modelComp, 'Nhomo versus Nhetero')
summary_bestModel = c(summary_bestModel, as.character(predicted_model$allocation))
summary_proba = c(summary_proba, predicted_model$post.prob)

# summarized output
summary_outfile = paste(timeStamp, '/', sub_dir_sim, '/hierarchical_models.txt', sep='')
write(paste(c(summary_modelComp), collapse='\t'), summary_outfile, append=F)
write(paste(c(summary_bestModel), collapse='\t'), summary_outfile, append=T)
write(paste(c(summary_proba), collapse='\t'), summary_outfile, append=T)

# best model
best = as.character(predicted_model_cons_grow$allocation)
if(predicted_model$allocation=='homoN'){
	best = paste(best, '1N', sep='_')
}else{
	best = paste(best, '2N', sep='_')
}

write(paste(best, '\n', sep=''), outfile_best, append=F)

