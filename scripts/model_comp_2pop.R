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
}
#nameA = 'txn'
#nameB = 'ama'
#ntree = 1000
#nSubdir = 6
nsims_monolocus = 10000 # number of monolocus simulations

outfile = paste(timeStamp, '/', sub_dir_sim, '/report_', nameA, '_', nameB, '.txt', sep='')
outfile_best = paste(timeStamp, '/', sub_dir_sim, '/best_model.txt', sep='')

# colors
coul = c('#ffffcc', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84')
coul = colorRampPalette(coul)

# observed data
obs_ss = read.table(paste(timeStamp, '/ABCstat_global.txt', sep=''), h=T)
obs_ss = obs_ss[, -grep('min', colnames(obs_ss))]
obs_ss = obs_ss[, -grep('max', colnames(obs_ss))]
if( outgroup == 1 ){
	obs_sfs = read.table(paste(timeStamp, '/ABCjsfs.txt', sep=''), h=T)
	ss_obs = cbind(obs_ss, obs_sfs)
	
	sfs=matrix(as.numeric(obs_sfs), byrow=T, ncol=nMin+1)
	colnames(sfs) = paste('f', nameB, 0:nMin, sep='_')
	rownames(sfs) = paste('f', nameA, 0:nMin, sep='_')
	write.table(sfs, paste(timeStamp, '/sfs_table.txt', sep=''), col.names=T, row.names=T, sep='\t', quote=F)
	
	sfs[1,2]=0; sfs[2,1]=0 # remove the singletons, ONLY FOR THE REPRESENTATION
	pdf(paste(timeStamp, '/sfs_plot.pdf', sep=''), bg="white")
	image(log10(sfs), col=coul(100), xlab = paste('frequency in ', nameA, sep=''), ylab = paste('frequency in ', nameB, sep=''), cex=1.5, cex.axis=1.5, cex.lab=1.5)
	dev.off()
}else{
	ss_obs = obs_ss
}


# simulated data
models = c('SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N')
migration = c('migration', 'migration', 'migration', 'migration', 'isolation', 'isolation', 'isolation', 'isolation', 'migration', 'migration', 'migration', 'migration', 'isolation', 'isolation')
homoM = c('1M', '1M', '2M', '2M', 'null', 'null', 'null', 'null', '1M', '1M', '2M', '2M', 'null', 'null')
homoN = c('1N', '2N', '1N', '2N', '1N', '2N', '1N', '2N', '1N', '2N', '1N', '2N', '1N', '2N')
ss_sim = list()
params_sim = list()
all_models_sim = NULL

for(m in models){
	ss_sim_tmp = NULL
	params_sim_tmp = NULL
	for(rep in seq(0, nSubdir-1, 1)){
		# statistics
		tmp_ss = read.table(paste(timeStamp, '/', sub_dir_sim, '/', m, '_', rep, '/ABCstat.txt', sep=''), h=T)
		tmp_ss = tmp_ss[, -grep('min', colnames(tmp_ss))]
		tmp_ss = tmp_ss[, -grep('max', colnames(tmp_ss))]
		if( outgroup == 1 ){ tmp_sfs = read.table(paste(timeStamp, '/', sub_dir_sim, '/', m, '_', rep, '/ABCjsfs.txt', sep=''), h=T)
			tmp = cbind(tmp_ss, tmp_sfs)
			ss_sim_tmp = rbind(ss_sim_tmp, tmp)
		}else{
			ss_sim_tmp = rbind(ss_sim_tmp, tmp_ss)
		}
		
		# params
		tmp_params = read.table(paste(timeStamp, '/', sub_dir_sim, '/', m, '_', rep, '/priorfile.txt', sep=''), h=T)
		params_sim_tmp = rbind(params_sim_tmp, tmp_params)
	}
	# statistics
	ss_sim[[m]] = ss_sim_tmp 
	all_models_sim = rbind(all_models_sim, ss_sim_tmp)
	
	# params
	params_sim[[m]] = params_sim_tmp
}


# remove uninformative statistics: those with no variation
ss_2_remove = c(1)
for(m in models){
	for(i in 2:ncol(ss_obs)){
		if( sd(ss_sim[[m]][, i])==0 ){
			ss_2_remove = c(ss_2_remove, i)
		}
	}
}
ss_2_remove = unique(ss_2_remove)



# model comparison #1 --> all models
modIndexes = NULL
for(m in models){
	modIndexes = c(modIndexes, rep(m, nrow(ss_sim[[m]])))
}

mod = abcrf(modIndexes~., data = data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
predicted_model = predict(mod, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)

write('MODEL COMPARISON #1: 14 models', outfile, append=F)
write('#confusion matrix:', outfile, append=T)
write.table(mod$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
write(paste('\n#best model among 14 models: ', predicted_model$allocation, sep=''), outfile, append=T)
write(paste('#proba best model among 14 models: ', predicted_model$post.prob, sep=''), outfile, append=T)
write('\n#votes:', outfile, append=T)
write.table(t(as.matrix(predicted_model$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)

write(predicted_model$allocation, outfile_best, append=F)


# model comparison #2 --> two models: isolation versus migration
modIndexes = NULL
for(i in 1:length(models)){
	modIndexes = c(modIndexes, rep(migration[i], nrow(ss_sim[[models[i]]])))
}

mod_iso_mig = abcrf(modIndexes~., data = data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
predicted_model_iso_mig = predict(mod_iso_mig, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, all_models_sim[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)

write('\n#####\n\nMODEL COMPARISON #2: 2 models', outfile, append=T)
write('#confusion matrix:', outfile, append=T)
write.table(mod_iso_mig$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
write(paste('\n#best model between migration and isolation: ', predicted_model_iso_mig$allocation, sep=''), outfile, append=T)
write(paste('#proba best model between migration and isolation: ', predicted_model_iso_mig$post.prob, sep=''), outfile, append=T)
write('\n#votes:', outfile, append=T)
write.table(t(as.matrix(predicted_model_iso_mig$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)



# model comparison #3 --> two models: N homo versus N hetero
Nhomo = NULL
for( i in c('SC_1M_1N', 'SC_2M_1N', 'AM_1M_1N', 'AM_2M_1N', 'IM_1M_1N', 'IM_2M_1N', 'SI_1N')){
	Nhomo = rbind(Nhomo, ss_sim[[i]])
}

Nhetero = NULL
for( i in c('SC_1M_2N', 'SC_2M_2N', 'AM_1M_2N', 'AM_2M_2N', 'IM_1M_2N', 'IM_2M_2N', 'SI_2N')){
	Nhetero = rbind(Nhetero, ss_sim[[i]])
}

modIndexes = c(rep('Nhomo', nrow(Nhomo)), rep('Nhetero', nrow(Nhetero)))
mod_Nhomo_Nhetero = abcrf(modIndexes~., data = data.frame(modIndexes, rbind(Nhomo, Nhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)
predicted_model_Nhomo_Nhetero = predict(mod_Nhomo_Nhetero, data.frame(ss_obs[, -ss_2_remove]), training=data.frame(modIndexes, rbind(Nhomo, Nhetero)[, -ss_2_remove]), ntree = ntree, paral = T, ncores = ncores)

write('\n#####\n\nMODEL COMPARISON #3: 2 models (Nhomo versus Nhetero)', outfile, append=T)
write('#confusion matrix:', outfile, append=T)
write.table(mod_Nhomo_Nhetero$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
write(paste('\n#best model between Nhomo and Nhetero: ', predicted_model_Nhomo_Nhetero$allocation, sep=''), outfile, append=T)
write(paste('#proba best model between Nhomo and Nhetero: ', predicted_model_Nhomo_Nhetero$post.prob, sep=''), outfile, append=T)
write('\n#votes:', outfile, append=T)
write.table(t(as.matrix(predicted_model_Nhomo_Nhetero$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)



# model comparison #4 --> two models: mig homo versus mig hetero
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

write('\n#####\n\nMODEL COMPARISON #4: 2 models (Mhomo versus Mhetero)', outfile, append=T)
write('#confusion matrix:', outfile, append=T)
write.table(mod_Mhomo_Mhetero$model.rf$confusion.matrix, outfile, append=T, col.names=T, row.names=T, sep='\t', quote=F)
write(paste('\n#best model between Mhomo and Mhetero: ', predicted_model_Mhomo_Mhetero$allocation, sep=''), outfile, append=T)
write(paste('#proba best model between Mhomo and Mhetero: ', predicted_model_Mhomo_Mhetero$post.prob, sep=''), outfile, append=T)
write('\n#votes:', outfile, append=T)
write.table(t(as.matrix(predicted_model_Mhomo_Mhetero$vote, ncol=1)), outfile, append=T, col.names=F, row.names=T, sep='\t', quote=F)


#######################################################
## parameters of the best model, IM_2M_2N and SI_2N
#nSimulations = 980000
#source('/shared/home/croux/softwares/ABConline/2pops/get_parameters.R')
#best_model = predicted_model$allocation
#write(paste('\n#####\n\nparameters of the best model: ', best_model, sep=''), outfile, append=T)
#posterior = get_posterior(nameA, nameB, nSubdir, best_model, nSimulations=nSimulations)
#write('param\tHPD2.5%\tmedian\tHPD%97.5', outfile, append=T)
#for(i in 1:ncol(posterior)){
#	write(paste(colnames(posterior)[i], as.numeric(quantile(posterior[,i], 0.025)), as.numeric(quantile(posterior[,i], 0.5)), as.numeric(quantile(posterior[,i], 0.975)), sep='\t'), outfile, append=T)
#}
#
## IM_2M_2N
#model_tmp = 'IM_2M_2N'
#if(best_model!=model_tmp){
#	write(paste('\n#####\n\nparameters of model: ', model_tmp, sep=''), outfile, append=T)
#	posterior = get_posterior(nameA, nameB, nSubdir, model_tmp, nSimulations=nSimulations)
#	write('param\tHPD2.5%\tmedian\tHPD%97.5', outfile, append=T)
#	for(i in 1:ncol(posterior)){
#		write(paste(colnames(posterior)[i], as.numeric(quantile(posterior[,i], 0.025)), as.numeric(quantile(posterior[,i], 0.5)), as.numeric(quantile(posterior[,i], 0.975)), sep='\t'), outfile, append=T)
#	}
#}
#
## SI_2N
#model_tmp = 'SI_2N'
#if(best_model!=model_tmp){
#	write(paste('\n#####\n\nparameters of model: ', model_tmp, sep=''), outfile, append=T)
#	posterior = get_posterior(nameA, nameB, nSubdir, model_tmp, nSimulations=nSimulations)
#	write('param\tHPD2.5%\tmedian\tHPD%97.5', outfile, append=T)
#	for(i in 1:ncol(posterior)){
#		write(paste(colnames(posterior)[i], as.numeric(quantile(posterior[,i], 0.025)), as.numeric(quantile(posterior[,i], 0.5)), as.numeric(quantile(posterior[,i], 0.975)), sep='\t'), outfile, append=T)
#	}
#}
#
