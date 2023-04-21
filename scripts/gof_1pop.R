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
        if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] }
        if(tmp[[1]][1] == 'sub_dir'){ sub_dir = tmp[[1]][2] }
}

### Summary Stats
# simulations
x = read.table(paste(timeStamp, '/', sub_dir, '/simulations.txt', sep=''), h=T)

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
exp_sfs = read.table(paste(timeStamp, '/', sub_dir, '/simulations_jsfs.txt', sep=''), h=T)
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
