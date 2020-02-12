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

for(i in commandArgs()){
        tmp = strsplit(i, '=')
        if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] } # name of the directory where the project is written 
        if(tmp[[1]][1] == 'sub_dir'){ sub_dir = tmp[[1]][2] } # name of the directory where the simulations were run
        if(tmp[[1]][1] == 'nIterations_gof'){ nIterations_gof = as.integer(tmp[[1]][2]) } # number of iterations : ~/timeStamp/sub_dir_{i} where i = range(nIterations_gof)
        if(tmp[[1]][1] == 'writeDistribution'){ writeDistribution = as.logical(tmp[[1]][2]) }
}

### Summary Stats
# simulations
x = NULL
for(i in 0:(nIterations_gof-1)){
        x = rbind(x, read.table(paste(timeStamp, '/', sub_dir, '/gof_', i, '/ABCstat.txt', sep=''), h=T))
}

# observation
y = read.table(paste(timeStamp, '/ABCstat_global.txt', sep=''), h=T)

# function to compute the pval
pvalue = function(distribution, obs){
	median_x = median(distribution)
	if(obs==median_x){
		pval=0.5
	}else if(as.numeric(obs)>median_x){
		pval = length(which(distribution>as.numeric(obs)))/length(distribution)
	}else{
		pval = length(which(distribution<as.numeric(obs)))/length(distribution)
	}
	return(pval)
}

if( writeDistribution==TRUE){

	### Summary Stats
	prior_ss = gof1_ss = gof2_ss = NULL
	prior_sfs = gof1_sfs = gof2_sfs = NULL
	
	# simulations
	for(i in 0:(nIterations_gof-1)){
		prior_ss = rbind(prior_ss, read.table(paste(timeStamp, '/best_model/best_model_', i, '/ABCstat.txt', sep=''), h=T))
		gof1_ss = rbind(gof1_ss, read.table(paste(timeStamp, '/gof/gof_', i, '/ABCstat.txt', sep=''), h=T))
		gof2_ss = rbind(gof2_ss, read.table(paste(timeStamp, '/gof_2/gof_', i, '/ABCstat.txt', sep=''), h=T))
		prior_sfs = rbind(prior_sfs, read.table(paste(timeStamp, '/best_model/best_model_', i, '/ABCjsfs.txt', sep=''), h=T))
		gof1_sfs = rbind(gof1_sfs, read.table(paste(timeStamp, '/gof/gof_', i, '/ABCjsfs.txt', sep=''), h=T))
		gof2_sfs = rbind(gof2_sfs, read.table(paste(timeStamp, '/gof_2/gof_', i, '/ABCjsfs.txt', sep=''), h=T))
	}

	# sub sample the stats in order to reduce the size
	finalSize = 2000
	if( nrow(prior_ss) > finalSize ){
		sub_prior = sample(1:nrow(prior_ss), finalSize, replace=F)
		prior_ss = prior_ss[sub_prior, ]
		prior_sfs = prior_sfs[sub_prior, ]
	}

	if( nrow(gof1_ss) > finalSize ){
		sub_gof1 = sample(1:nrow(gof1_ss), finalSize, replace=F)
		gof1_ss = gof1_ss[sub_gof1, ]
		gof1_sfs = gof1_sfs[sub_gof1, ]
	}
	
	if( nrow(gof2_ss) > finalSize ){
		sub_gof2 = sample(1:nrow(gof2_ss), finalSize, replace=F)
		gof2_ss = gof2_ss[sub_gof2, ]
		gof2_sfs = gof2_sfs[sub_gof2, ]
	}
	
	# observation
	obs_ss = read.table(paste(timeStamp, '/ABCstat_global.txt', sep=''), h=T)
	obs_sfs = read.table(paste( timeStamp, '/ABCjsfs.txt', sep=''), h=T)

	# output for the web interface
	origin = c('observed dataset', rep('prior', nrow(prior_ss)), rep('posterior', nrow(gof1_ss)), rep('optimized posterior', nrow(gof2_ss)))
	PCA = rbind(cbind(obs_ss, obs_sfs), cbind(prior_ss, prior_sfs), cbind(gof1_ss, gof1_sfs), cbind(gof2_ss, gof2_sfs))
	PCA = cbind(PCA, origin)
	write.table(PCA, paste(timeStamp, '/distribution_PCA.txt', sep=''), col.names=T, row.names=F, quote=F, sep='\t')
}


# measure the pval over all statistics
ss = c(2:31, 40:48)

stats = NULL
pvals = NULL
mean_exp = NULL
mean_obs = NULL
for(i in ss){
	stats = c(stats, colnames(x)[i])
	pvals = c(pvals, round(pvalue(x[,i], y[i]), 5))
	mean_exp = c(mean_exp, round(mean(x[,i]), 5))
	mean_obs = c(mean_obs, round(as.numeric(y[i]), 5))
}

pvals_fdr_corrected = round(p.adjust(pvals, "fdr"), 5)

res = data.frame(stats, mean_exp, mean_obs, pvals_fdr_corrected)

# outfile
outfile = paste(timeStamp, "/", sub_dir, "/goodness_of_fit_test.txt", sep='')
write.table(x=res, file=outfile, quote=FALSE, sep='\t', col.names=T, row.names=F)


### jSFS
# expected sfs
exp_sfs = NULL
for(i in 0:(nIterations_gof-1)){
        exp_sfs = rbind(exp_sfs, read.table(paste(timeStamp, '/', sub_dir, '/gof_', i, '/ABCjsfs.txt', sep=''), h=T))
}
exp_sfs_2 = apply(exp_sfs, MARGIN=2, FUN="median")

# observed sfs
obs_sfs = read.table(paste(timeStamp, '/ABCjsfs.txt', sep=''), h=T)

# compute the pvalue
tested_sfs = NULL
for(i in 1:length(exp_sfs)){
	tested_sfs = c(tested_sfs, pvalue(exp_sfs[,i], obs_sfs[i]))
}

tested_sfs = round(p.adjust(tested_sfs, "fdr"), 5)

# all matrixes 
## obs | exp | exp-obs | pval
sfs = rbind(obs_sfs, exp_sfs_2, exp_sfs_2-obs_sfs, tested_sfs)

outfile_sfs = paste(timeStamp, "/", sub_dir, "/gof_sfs.txt", sep='')
write.table(x=sfs, file=outfile_sfs, quote=FALSE, sep='\t', col.names=T, row.names=F)
 
