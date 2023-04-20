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
# model_comp_2pop.R nameA=txn nameB=ama nSubdir=20  ntree=1000 
for(i in commandArgs()){
	tmp = strsplit(i, '=')
	if(tmp[[1]][1] == 'nameA'){ nameA = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nameB'){ nameB = tmp[[1]][2] }
	if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nMin'){ nMin = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'sub_dir_sim'){ sub_dir_sim = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nSubdir'){ nSubdir = as.integer(tmp[[1]][2]) } # number of subdirectories where simulations were ran
	if(tmp[[1]][1] == 'ncores'){ ncores = as.integer(tmp[[1]][2]) } # number of cores for the random forest
	if(tmp[[1]][1] == 'ntree'){ ntree = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'outgroup'){ outgroup = as.integer(tmp[[1]][2]) } # 0: no outgroup, no SFS used. 1: outgroup, SFS used
	if(tmp[[1]][1] == 'population_growth'){ population_growth = tmp[[1]][2] } # constant; variable
	if(tmp[[1]][1] == 'modeBarrier'){ modeBarrier = tmp[[1]][2] } # beta; bimodal
	if(tmp[[1]][1] == 'binpath'){ binpath = tmp[[1]][2] } # path to the bin directory
	if(tmp[[1]][1] == 'posterior2use'){ posterior2use = tmp[[1]][2] } # path to the posterior file used for the locus specific model comp
}

outfile = paste(timeStamp, '/', sub_dir_sim, '/report_', nameA, '_', nameB, '.txt', sep='')
outfile_best = paste(timeStamp, '/', sub_dir_sim, '/best_model.txt', sep='')

# colors
coul = c('#ffffcc', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84')
coul = colorRampPalette(coul)

# observed data
obs_ss = read.table(paste(timeStamp, '/ABCstat_global.txt', sep=''), h=T) # global statistics over all loci (avg and std)
obs_ss = obs_ss[, -grep('min', colnames(obs_ss))]
obs_ss = obs_ss[, -grep('max', colnames(obs_ss))]

obs_loci = read.table(paste(timeStamp, '/ABCstat_loci.txt', sep=''), h=T) # individual statistics for each locus
obs_sfs = read.table(paste(timeStamp, '/ABCjsfs.txt', sep=''), h=T)

ss_obs = cbind(obs_ss, obs_sfs)

sfs=matrix(as.numeric(obs_sfs), byrow=T, ncol=nMin+1)
colnames(sfs) = paste('f', nameB, 0:nMin, sep='_')
rownames(sfs) = paste('f', nameA, 0:nMin, sep='_')
write.table(sfs, paste(timeStamp, '/sfs_table.txt', sep=''), col.names=T, row.names=T, sep='\t', quote=F)

sfs[1,2]=0; sfs[2,1]=0 # remove the singletons, ONLY FOR THE REPRESENTATION
pdf(paste(timeStamp, '/sfs_plot.pdf', sep=''), bg="white")
image(log10(sfs), col=rev(viridis(100)), xlab = paste('frequency in ', nameA, sep=''), ylab = paste('frequency in ', nameB, sep=''), cex=1.5, cex.axis=1.5, cex.lab=1.5)
dev.off()


# simulated data
models = c('SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N', 'PAN_1N', 'PAN_2N')
migration = c('migration', 'migration', 'migration', 'migration', 'isolation', 'isolation', 'isolation', 'isolation', 'migration', 'migration', 'migration', 'migration', 'isolation', 'isolation', 'migration', 'migration')
homoM = c('1M', '1M', '2M', '2M', 'null', 'null', 'null', 'null', '1M', '1M', '2M', '2M', 'null', 'null', 'null', 'null')
homoN = c('1N', '2N', '1N', '2N', '1N', '2N', '1N', '2N', '1N', '2N', '1N', '2N', '1N', '2N', '1N', '2N')
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
			tmp_ss = tmp_ss[, -grep('min', colnames(tmp_ss))]
			tmp_ss = tmp_ss[, -grep('max', colnames(tmp_ss))]
			
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
	modIndexes = c(modIndexes, rep(migration[i], nrow(ss_sim[[models[i]]])))
}

mod_iso_mig = abcrf(modIndexes~., data = data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
predicted_model_iso_mig = predict(mod_iso_mig, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)

write('\n#####\n\nMODEL COMPARISON #1: migration versus isolation', outfile, append=F)
write('#confusion matrix:', outfile, append=T)
write.table(mod_iso_mig$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
write(paste('\n#best model: ', predicted_model_iso_mig$allocation, sep=''), outfile, append=T)
write(paste('#proba best model: ', predicted_model_iso_mig$post.prob, sep=''), outfile, append=T)
write('\n#votes:', outfile, append=T)
write.table(t(as.matrix(predicted_model_iso_mig$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

summary_modelComp = c(summary_modelComp, 'migration versus isolation')
summary_bestModel = c(summary_bestModel, as.character(predicted_model_iso_mig$allocation))
summary_proba = c(summary_proba, predicted_model_iso_mig$post.prob)



# if ongoing migration
if(predicted_model_iso_mig$allocation=='migration'){
	# model comparison #2 --> IM versus SC
	IM = c('IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N')
	SC = c('SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N')
	PAN = c('PAN_1N', 'PAN_2N')
	all_models = NULL
	modIndexes = NULL
	for(i in IM){
		all_models = rbind(all_models, ss_sim[[i]])
		modIndexes = c(modIndexes, rep('IM', nrow(ss_sim[[i]])))
	}
	for(i in SC){
		all_models = rbind(all_models, ss_sim[[i]])
		modIndexes = c(modIndexes, rep('SC', nrow(ss_sim[[i]])))
	}
	for(i in PAN){
		all_models = rbind(all_models, ss_sim[[i]])
		modIndexes = c(modIndexes, rep('PAN', nrow(ss_sim[[i]])))
	}
	
	
	mod = abcrf(modIndexes~., data = data.frame(modIndexes, all_models[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
	predicted_model = predict(mod, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, all_models[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
	
	write('\n#####\n\nMODEL COMPARISON #2: IM versus SC versus PAN', outfile, append=T)
	write('#confusion matrix:', outfile, append=T)
	write.table(mod$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
	write(paste('\n#best model: ', predicted_model$allocation, sep=''), outfile, append=T)
	write(paste('#proba best model: ', predicted_model$post.prob, sep=''), outfile, append=T)
	write('\n#votes:', outfile, append=T)
	write.table(t(as.matrix(predicted_model$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

	
	summary_modelComp = c(summary_modelComp, 'IM versus SC versus PAN')
	summary_bestModel = c(summary_bestModel, as.character(predicted_model$allocation))
	summary_proba = c(summary_proba, predicted_model$post.prob)

	best_hierarchical = as.character(predicted_model$allocation)
	

	if(as.character(predicted_model$allocation) != 'PAN'){
		# model comparison #3 --> two models: mig homo versus mig hetero
		Mhomo = NULL
		for( i in c('SC_1M_1N', 'SC_1M_2N', 'IM_1M_1N', 'IM_1M_2N')){
			Mhomo = rbind(Mhomo, ss_sim[[i]])
		}

		Mhetero = NULL
		for( i in c('SC_2M_1N', 'SC_2M_2N', 'IM_2M_1N', 'IM_2M_2N')){
			Mhetero = rbind(Mhetero, ss_sim[[i]])
		}

		modIndexes = c(rep('Mhomo', nrow(Mhomo)), rep('Mhetero', nrow(Mhetero)))
		mod_Mhomo_Mhetero = abcrf(modIndexes~., data = data.frame(modIndexes, rbind(Mhomo, Mhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
		predicted_model_Mhomo_Mhetero = predict(mod_Mhomo_Mhetero, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, rbind(Mhomo, Mhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)

		write('\n#####\n\nMODEL COMPARISON #3: Mhomo versus Mhetero', outfile, append=T)
		write('#confusion matrix:', outfile, append=T)
		write.table(mod_Mhomo_Mhetero$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
		write(paste('\n#best model: ', predicted_model_Mhomo_Mhetero$allocation, sep=''), outfile, append=T)
		write(paste('#proba best model: ', predicted_model_Mhomo_Mhetero$post.prob, sep=''), outfile, append=T)
		write('\n#votes:', outfile, append=T)
		write.table(t(as.matrix(predicted_model_Mhomo_Mhetero$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)
		
		summary_modelComp = c(summary_modelComp, 'M-homo versus M-hetero')
		summary_bestModel = c(summary_bestModel, as.character(predicted_model_Mhomo_Mhetero$allocation))
		summary_proba = c(summary_proba, predicted_model_Mhomo_Mhetero$post.prob)

		if(as.character(predicted_model_Mhomo_Mhetero$allocation) == 'Mhomo'){
			best_hierarchical = paste(best_hierarchical, '1M', sep='_')
		}else{
			best_hierarchical = paste(best_hierarchical, '2M', sep='_')

		}
		
		# model comparison #4 --> Nhomo versus Nhetero
		Nhomo = NULL
		for( i in c('SC_1M_1N', 'SC_2M_1N', 'IM_1M_1N', 'IM_2M_1N')){
			Nhomo = rbind(Nhomo, ss_sim[[i]])
		}

		Nhetero = NULL
		for( i in c('SC_1M_2N', 'SC_2M_2N', 'IM_1M_2N', 'IM_2M_2N')){
			Nhetero = rbind(Nhetero, ss_sim[[i]])
		}

		modIndexes = c(rep('Nhomo', nrow(Nhomo)), rep('Nhetero', nrow(Nhetero)))
		mod_Nhomo_Nhetero = abcrf(modIndexes~., data = data.frame(modIndexes, rbind(Nhomo, Nhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
		predicted_model_Nhomo_Nhetero = predict(mod_Nhomo_Nhetero, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, rbind(Nhomo, Nhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
		
		
		write('\n#####\n\nMODEL COMPARISON #4: Nhomo versus Nhetero', outfile, append=T)
		write('#confusion matrix:', outfile, append=T)
		write.table(mod_Nhomo_Nhetero$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
		write(paste('\n#best model: ', predicted_model_Nhomo_Nhetero$allocation, sep=''), outfile, append=T)
		write(paste('#proba best model: ', predicted_model_Nhomo_Nhetero$post.prob, sep=''), outfile, append=T)
		write('\n#votes:', outfile, append=T)
		write.table(t(as.matrix(predicted_model_Nhomo_Nhetero$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

		summary_modelComp = c(summary_modelComp, 'N-homo versus N-hetero')
		summary_bestModel = c(summary_bestModel, as.character(predicted_model_Nhomo_Nhetero$allocation))
		summary_proba = c(summary_proba, predicted_model_Nhomo_Nhetero$post.prob)

		if(as.character(predicted_model_Nhomo_Nhetero$allocation) == 'Nhomo'){
			best_hierarchical = paste(best_hierarchical, '1N', sep='_')
		}else{
			best_hierarchical = paste(best_hierarchical, '2N', sep='_')

		}
	}else{ # if PAN model
		# model comparison #4 --> Nhomo versus Nhetero
		Nhomo = NULL
		for( i in c('PAN_1N')){
			Nhomo = rbind(Nhomo, ss_sim[[i]])
		}

		Nhetero = NULL
		for( i in c('PAN_2N')){
			Nhetero = rbind(Nhetero, ss_sim[[i]])
		}

		modIndexes = c(rep('Nhomo', nrow(Nhomo)), rep('Nhetero', nrow(Nhetero)))
		mod_Nhomo_Nhetero = abcrf(modIndexes~., data = data.frame(modIndexes, rbind(Nhomo, Nhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
		predicted_model_Nhomo_Nhetero = predict(mod_Nhomo_Nhetero, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, rbind(Nhomo, Nhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
		
		
		write('\n#####\n\nMODEL COMPARISON #2: Nhomo versus Nhetero', outfile, append=T)
		write('#confusion matrix:', outfile, append=T)
		write.table(mod_Nhomo_Nhetero$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
		write(paste('\n#best model: ', predicted_model_Nhomo_Nhetero$allocation, sep=''), outfile, append=T)
		write(paste('#proba best model: ', predicted_model_Nhomo_Nhetero$post.prob, sep=''), outfile, append=T)
		write('\n#votes:', outfile, append=T)
		write.table(t(as.matrix(predicted_model_Nhomo_Nhetero$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

		summary_modelComp = c(summary_modelComp, 'N-homo versus N-hetero')
		summary_bestModel = c(summary_bestModel, as.character(predicted_model_Nhomo_Nhetero$allocation))
		summary_proba = c(summary_proba, predicted_model_Nhomo_Nhetero$post.prob)

		if(as.character(predicted_model_Nhomo_Nhetero$allocation) == 'Nhomo'){
			best_hierarchical = paste(best_hierarchical, '1N', sep='_')
		}else{
			best_hierarchical = paste(best_hierarchical, '2N', sep='_')

		}
	}
	write(paste(best_hierarchical, '\n', sep=''), outfile_best, append=F)
}



# if ongoing migration
if(predicted_model_iso_mig$allocation=='isolation'){
	AM = c('AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N')
	SI = c('SI_1N', 'SI_2N')
	# model comparison #2 --> IM versus SC
	all_models = NULL
	modIndexes = NULL
	for(i in AM){
		all_models = rbind(all_models, ss_sim[[i]])
		modIndexes = c(modIndexes, rep('AM', nrow(ss_sim[[i]])))
	}
	for(i in SI){
		all_models = rbind(all_models, ss_sim[[i]])
		modIndexes = c(modIndexes, rep('SI', nrow(ss_sim[[i]])))
	}
	
	mod = abcrf(modIndexes~., data = data.frame(modIndexes, all_models[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
	predicted_model = predict(mod, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, all_models[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
	
	write('\n#####\n\nMODEL COMPARISON #2: AM versus SI', outfile, append=T)
	write('#confusion matrix:', outfile, append=T)
	write.table(mod$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
	write(paste('\n#best model: ', predicted_model$allocation, sep=''), outfile, append=T)
	write(paste('#proba best model: ', predicted_model$post.prob, sep=''), outfile, append=T)
	write('\n#votes:', outfile, append=T)
	write.table(t(as.matrix(predicted_model$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

	summary_modelComp = c(summary_modelComp, 'AM versus SI')
	summary_bestModel = c(summary_bestModel, as.character(predicted_model$allocation))
	summary_proba = c(summary_proba, predicted_model$post.prob)
	
	
	# model comparison #3 --> Nhomo versus Nhetero
	Nhomo = NULL
	for( i in c('AM_1M_1N', 'AM_2M_1N', 'SI_1N')){
		Nhomo = rbind(Nhomo, ss_sim[[i]])
	}

	Nhetero = NULL
	for( i in c('AM_1M_2N', 'AM_2M_2N', 'SI_2N')){
		Nhetero = rbind(Nhetero, ss_sim[[i]])
	}

	modIndexes = c(rep('Nhomo', nrow(Nhomo)), rep('Nhetero', nrow(Nhetero)))
	mod_Nhomo_Nhetero = abcrf(modIndexes~., data = data.frame(modIndexes, rbind(Nhomo, Nhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
	predicted_model_Nhomo_Nhetero = predict(mod_Nhomo_Nhetero, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, rbind(Nhomo, Nhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
	
	
	write('\n#####\n\nMODEL COMPARISON #3: Nhomo versus Nhetero', outfile, append=T)
	write('#confusion matrix:', outfile, append=T)
	write.table(mod_Nhomo_Nhetero$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
	write(paste('\n#best model: ', predicted_model_Nhomo_Nhetero$allocation, sep=''), outfile, append=T)
	write(paste('#proba best model: ', predicted_model_Nhomo_Nhetero$post.prob, sep=''), outfile, append=T)
	write('\n#votes:', outfile, append=T)
	write.table(t(as.matrix(predicted_model_Nhomo_Nhetero$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)
	
	summary_modelComp = c(summary_modelComp, 'N-homo versus N-hetero')
	summary_bestModel = c(summary_bestModel, as.character(predicted_model_Nhomo_Nhetero$allocation))
	summary_proba = c(summary_proba, predicted_model_Nhomo_Nhetero$post.prob)
	
	
	# model comparison #4 --> within demographic scenario
	if( predicted_model$allocation == 'SI' ){
		write('\n#####\n\nMODEL COMPARISON #4: within SI', outfile, append=T)
		SI_1N = ss_sim[['SI_1N']]
		SI_2N = ss_sim[['SI_2N']]
		
		modIndexes = c(rep('SI_1N', nrow(SI_1N)), rep('SI_2N', nrow(SI_2N)))
		mod_best = abcrf(modIndexes~., data = data.frame(modIndexes, rbind(SI_1N, SI_2N)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
		predicted_best = predict(mod_best, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, rbind(SI_1N, SI_2N)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
		}else if( predicted_model$allocation == 'AM' ){
		write('\n#####\n\nMODEL COMPARISON #4: within AM', outfile, append=T)
		AM_1M_1N = ss_sim[['AM_1M_1N']]
		AM_1M_2N = ss_sim[['AM_1M_2N']]
		AM_2M_1N = ss_sim[['AM_2M_1N']]
		AM_2M_2N = ss_sim[['AM_2M_2N']]
		
		modIndexes = c(rep('AM_1M_1N', nrow(AM_1M_1N)), rep('AM_1M_2N', nrow(AM_1M_2N)), rep('AM_2M_1N', nrow(AM_2M_1N)), rep('AM_2M_2N', nrow(AM_2M_2N)))
		mod_best = abcrf(modIndexes~., data = data.frame(modIndexes, rbind(AM_1M_1N, AM_1M_2N, AM_2M_1N, AM_2M_2N)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
		predicted_best = predict(mod_best, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, rbind(AM_1M_1N, AM_1M_2N, AM_2M_1N, AM_2M_2N)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
		}
		
	write('#confusion matrix:', outfile, append=T)
	write.table(mod_best$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
	write(paste('\n#best model: ', predicted_best$allocation, sep=''), outfile, append=T)
	write(paste('#proba best model: ', predicted_best$post.prob, sep=''), outfile, append=T)
	write('\n#votes:', outfile, append=T)
	write.table(t(as.matrix(predicted_best$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)
	write(paste(predicted_best$allocation, '\n', sep=''), outfile_best, append=F)
}


# summarized output
summary_outfile = paste(timeStamp, '/', sub_dir_sim, '/hierarchical_models.txt', sep='')
write(paste(c(summary_modelComp), collapse='\t'), summary_outfile, append=F)
write(paste(c(summary_bestModel), collapse='\t'), summary_outfile, append=T)
write(paste(c(summary_proba), collapse='\t'), summary_outfile, append=T)


#### LOCUS SPECIFIC MODEL COMPARISON
#nrep = 10
#path = getwd()
#outfile = 'locus_specific_modelComp.txt'
#if(predicted_model_iso_mig$allocation=='migration'){
#	if(predicted_model_Mhomo_Mhetero$allocation=='Mhetero'){
#		# get the posterior
#		#posterior_IM = read.table(paste(timeStamp, "/estim/posterior_IM_2M_2N.txt", sep=''), h=T)
#		posterior_IM = read.table(posterior2use, h=T) # usualy : posterior of an IM model
#		
#		# get informations about loci 
#		bpfile = read.table(paste(timeStamp, "/bpfile", sep=''), skip=1, h=F)
#		L = median(as.numeric(bpfile[1,]))
#		nA = median(as.numeric(bpfile[2,]))
#		nB = median(as.numeric(bpfile[3,]))
#		theta = median(as.numeric(bpfile[4,]))
#		rho = median(as.numeric(bpfile[5,]))
#
#		# change directory
#		setwd(paste(timeStamp, '/modelComp', sep=''))
#
#		# write the monolocus bpfile
#		write(paste('#locus specific model comparison\n', L, '\n', nA, '\n', nB, '\n', theta, '\n', rho, sep=''), 'bpfile', append=F)
#		
#		## simulations
#		N1_mig = NULL
#		N2_mig = NULL
#		Na_mig = NULL
#		Tsplit_mig = NULL
#		M12_mig = NULL
#		M21_mig = NULL
#		
#		N1_iso = NULL
#		N2_iso = NULL
#		Na_iso = NULL
#		Tsplit_iso = NULL
#		M12_iso = NULL
#		M21_iso = NULL
#		
#		Tsplit = NULL
#		
#		if(population_growth=='variable'){
#			Tdem1_mig = NULL
#			Tdem2_mig = NULL
#			founders1_mig = NULL
#			founders2_mig = NULL
#			
#			Tdem1_iso = NULL
#			Tdem2_iso = NULL
#			founders1_iso = NULL
#			founders2_iso = NULL
#		}
#		
#		for(i in 1:nrow(posterior_IM)){
#			# population growth
#			if(population_growth=='variable'){
#				Tdem1_mig = c(Tdem1_mig, rep(posterior_IM$Tdem1[i], nrep))
#				Tdem2_mig = c(Tdem2_mig, rep(posterior_IM$Tdem2[i], nrep))
#				founders1_mig = c(founders1_mig, rep(posterior_IM$founders1[i], nrep))
#				founders2_mig = c(founders2_mig, rep(posterior_IM$founders2[i], nrep))
#				
#				Tdem1_iso = c(Tdem1_iso, rep(posterior_IM$Tdem1[i], nrep))
#				Tdem2_iso = c(Tdem2_iso, rep(posterior_IM$Tdem2[i], nrep))
#				founders1_iso = c(founders1_iso, rep(posterior_IM$founders1[i], nrep))
#				founders2_iso = c(founders2_iso, rep(posterior_IM$founders2[i], nrep))
#			}
#			
#			Tsplit = c(Tsplit, rep(posterior_IM$Tsplit[i], nrep))
#			
#			# migration
#			a_N = posterior_IM$shape_N_a[i]
#			b_N = posterior_IM$shape_N_b[i]
#			scalar_N = rbeta(nrep, a_N, b_N) / (a_N / (a_N + b_N))
#			N1_mig = c(N1_mig, scalar_N * posterior_IM$N1[i])
#			N2_mig = c(N2_mig, scalar_N * posterior_IM$N2[i])
#			Na_mig = c(Na_mig, scalar_N * posterior_IM$Na[i])
#		
#			if(modeBarrier == "beta"){
#				a_M1 = posterior_IM$shape_M12_a[i]
#				b_M1 = posterior_IM$shape_M12_b[i]
#				a_M2 = posterior_IM$shape_M21_a[i]
#				b_M2 = posterior_IM$shape_M21_b[i]
#				scalar_M12 = rbeta(nrep, a_M1, b_M1) / (a_M1 / (a_M1 + b_M1))
#				scalar_M21 = rbeta(nrep, a_M2, b_M2) / (a_M2 / (a_M2 + b_M2))
#			}else{
#				# bimodal
#				scalar_M12 = rep(1, nrep)
#				scalar_M21 = rep(1, nrep)
#			}
#			M12_mig = c(M12_mig, scalar_M12 * posterior_IM$M12[i])
#			M21_mig = c(M21_mig, scalar_M21 * posterior_IM$M21[i])
#			
#			# isolation
#			a_N = posterior_IM$shape_N_a[i]
#			b_N = posterior_IM$shape_N_b[i]
#			scalar_N = rbeta(nrep, posterior_IM$shape_N_a[i], posterior_IM$shape_N_b[i]) / (a_N / (a_N + b_N))
#			N1_iso = c(N1_iso, scalar_N * posterior_IM$N1[i])
#			N2_iso = c(N2_iso, scalar_N * posterior_IM$N2[i])
#			Na_iso = c(Na_iso, scalar_N * posterior_IM$Na[i])
#			
#			M12_iso = c(M12_iso, rep(0, nrep))
#			M21_iso = c(M21_iso, rep(0, nrep))
#		}
#		
#		theta_prior = rep(theta, length(Tsplit))
#		rho_prior = rep(rho, length(Tsplit))
#		L_prior = rep(L, length(Tsplit))
#		nA_prior = rep(nA, length(Tsplit))
#		nB_prior = rep(nB, length(Tsplit))
#		ntot = nA_prior + nB_prior
#		
#		if(population_growth=='variable'){
#			# msnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs
#			prior_mig = cbind(ntot, theta_prior, rho_prior, L_prior, nA_prior, nB_prior, N1_mig, N2_mig, Tdem1_mig, founders1_mig*Na_mig, Tdem2_mig, founders2_mig*Na_mig, M12_mig, M21_mig, Tsplit, Tsplit, Na_mig)
#			prior_iso = cbind(ntot, theta_prior, rho_prior, L_prior, nA_prior, nB_prior, N1_iso, N2_iso, Tdem1_iso, founders1_iso*Na_iso, Tdem2_iso, founders2_iso*Na_iso, M12_iso, M21_iso, Tsplit, Tsplit, Na_iso)
#		}else{
#			prior_mig = cbind(ntot, theta_prior, rho_prior, L_prior, nA_prior, nB_prior, N1_mig, N2_mig, M12_mig, M21_mig, Tsplit, Tsplit, Na_mig)
#			prior_iso = cbind(ntot, theta_prior, rho_prior, L_prior, nA_prior, nB_prior, N1_iso, N2_iso, M12_iso, M21_iso, Tsplit, Tsplit, Na_iso)
#		}
#		
#		# simulations of migration
#		# "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs"
#		write.table(prior_mig, 'prior_mig.txt', sep='\t', col.names=F, row.names=F, quote=F)
#		if(population_growth=='variable'){
#			commande = paste('cat prior_mig.txt | ', binpath, '/msnsam tbs ', nrow(prior_mig), ' -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs | ', binpath, '/mscalc_2pop.py', sep='')
#		}else{
#			commande = paste('cat prior_mig.txt | ', binpath, '/msnsam tbs ', nrow(prior_mig), ' -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs | ', binpath, '/mscalc_2pop.py', sep='')
#		}
#		system(commande)
#		mig_ss = read.table('ABCstat.txt', h=T)
#		
#		# simulations of isolation
#		# "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs"
#		write.table(prior_iso, 'prior_iso.txt', sep='\t', col.names=F, row.names=F, quote=F)
#		if(population_growth=='variable'){
#			commande = paste('cat prior_iso.txt | ', binpath, '/msnsam tbs ', nrow(prior_iso), ' -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs | ', binpath, '/mscalc_2pop.py', sep='')
#		}else{
#			commande = paste('cat prior_iso.txt | ', binpath, '/msnsam tbs ', nrow(prior_iso), ' -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs | ', binpath, '/mscalc_2pop.py', sep='')
#		}
#		system(commande)
#		iso_ss = read.table('ABCstat.txt', h=T)
#		
#		# locus_specific model comparison
#		stats_obs = c(3:15, 20:24)
#		stats_sim = c(4, 6, 8, 10, 12, 14, 16, 19, 21, 24, 26, 28, 30, 40, 45, 46, 47, 48)
#		modIndexes = c(rep('migration', nrow(mig_ss)), rep('isolation', nrow(iso_ss)))
#
#		mod_iso_mig = abcrf(modIndexes~., data = data.frame(modIndexes, rbind(mig_ss, iso_ss)[, stats_sim]), ntree = ntree, paral = T, ncores = ncores)
#		predicted_model_iso_mig = predict(mod_iso_mig, data.frame(obs_loci[, stats_obs]), training=data.frame(modIndexes, rbind(mig_ss, iso_ss)[, stats_sim]), ntree = ntree, paral = T, ncores = ncores)
#
#		allocation = predicted_model_iso_mig$allocation
#		post_proba = predicted_model_iso_mig$post.prob
#		res = data.frame(obs_loci[, c(1, stats_obs)], allocation, post_proba)
#		
#		write.table(res, outfile, col.names=T, row.names=F, quote=F, sep='\t', append=F)
#	}else{ # if migration but homogeneous
#		# change directory
#		setwd(paste(timeStamp, '/modelComp', sep=''))
#		stats_obs = c(3:15, 20:24)
#		allocation = rep('migration', nrow(obs_loci))
#		post_proba = rep('1', nrow(obs_loci))
#		res = data.frame(obs_loci[, c(1, stats_obs)], allocation, post_proba)
#		write.table(res, outfile, col.names=T, row.names=F, quote=F, sep='\t', append=F)
#	}
#}else{
#	# change directory
#	setwd(paste(timeStamp, '/modelComp', sep=''))
#	stats_obs = c(3:15, 20:24)
#	allocation = rep('isolation', nrow(obs_loci))
#	post_proba = rep('1', nrow(obs_loci))
#	res = data.frame(obs_loci[, c(1, stats_obs)], allocation, post_proba)
#	write.table(res, outfile, col.names=T, row.names=F, quote=F, sep='\t', append=F)
#}
#
