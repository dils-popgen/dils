#!/shared/software/miniconda/envs/r-3.5.1/bin/Rscript
# #!/shared/home/croux/.conda/envs/R_env/bin/Rscript
# #!/usr/bin/Rscript

writeDistribution = FALSE
nIterations_gof = 1

 for(i in commandArgs()){
        tmp = strsplit(i, '=')
        if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
        if(tmp[[1]][1] == 'sub_dir'){ sub_dir = tmp[[1]][2] }
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


if( writeDistribution==TRUE ){
	### Summary Stats
	prior_ss = gof1_ss = gof2_ss = NULL
	prior_sfs = gof1_sfs = gof2_sfs = NULL
	# simulations
	for(i in 0:(nIterations_gof-1)){
		prior_ss = rbind(prior_ss, read.table(paste(timeStamp, '/best_model/best_model_', i, '/ABCstat.txt', sep=''), h=T))
		gof1_ss = rbind(gof1_ss, read.table(paste(timeStamp, '/gof/gof_', i, '/ABCstat.txt', sep=''), h=T)) # posterior
		gof2_ss = rbind(gof2_ss, read.table(paste(timeStamp, '/gof_2/gof_', i, '/ABCstat.txt', sep=''), h=T)) # optimized posterior
		prior_sfs = rbind(prior_sfs, read.table(paste(timeStamp, '/best_model/best_model_', i, '/ABCjsfs.txt', sep=''), h=T))
		gof1_sfs = rbind(gof1_sfs, read.table(paste(timeStamp, '/gof/gof_', i, '/ABCjsfs.txt', sep=''), h=T))
		gof2_sfs = rbind(gof2_sfs, read.table(paste(timeStamp, '/gof_2/gof_', i, '/ABCjsfs.txt', sep=''), h=T))
	}
	
	# subselect
	npoints = 2000
	if(nrow(prior_ss)>npoints){
		select = sample(1:(nrow(prior_ss)), npoints, replace=F)
		prior_ss = prior_ss[select,]
		prior_sfs = prior_sfs[select,]
	}

	if(nrow(gof1_ss)>npoints){
		select = sample(1:(nrow(gof1_ss)), npoints, replace=F)
		gof1_ss = gof1_ss[select,]
		gof1_sfs = gof1_sfs[select,]
	}

	if(nrow(gof2_ss)>npoints){
		select = sample(1:(nrow(gof2_ss)), npoints, replace=F)
		gof2_ss = gof2_ss[select,]
		gof2_sfs = gof2_sfs[select,]
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
stats = NULL
pvals = NULL
mean_exp = NULL
mean_obs = NULL
for(i in 2:ncol(x)){
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

#exp_sfs = read.table(paste(timeStamp, '/', sub_dir, '/simulations_jsfs.txt', sep=''), h=T)
exp_sfs_2 = apply(exp_sfs, MARGIN=2, FUN="median")

# observed sfs
obs_sfs = read.table(paste(timeStamp, '/ABCjsfs.txt', sep=''), h=T)

# to remove
toRemove = apply(exp_sfs, MARGIN=2, FUN="sum")
toRemove = which(toRemove==0)

if(sum(toRemove) != 0){
	exp_sfs = exp_sfs[, -toRemove]
	obs_sfs = obs_sfs[, -toRemove]
	exp_sfs_2 = exp_sfs_2[-toRemove]
}

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
# barplot(as.matrix(sfs[1:2,]), beside=T, col=viridis(2)) 

