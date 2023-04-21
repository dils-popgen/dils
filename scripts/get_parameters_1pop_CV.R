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

# function to estimate the parameters
abc_nnet_multivar <- function(target,x,sumstat,tol,gwt,rejmethod=F,noweight=F,transf="none",bb=c(0,0),nb.nnet=10,size.nnet=5,trace=T, MaxNWts=10000){
	options(digits=5)
	require(nnet)
	# target is the set of target summary stats
	# x is the parameter vector (long vector of numbers from the simulations) and is the dependent variable for the regression
	# x can also be a matrix for multi-dimensional models. Each column corresponds to a parameter and each row to a simulation.
	# sumstat is an array of simulated summary stats (i.e. independent variables).
	# tol is the required proportion of points nearest the target values
	# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
	# if noweight=T, no Epanechnikov weights are calculated
	# if rejmethod=T it doesn't bother with the regression, and just does rejection.
	# nb.nnet>1 is the number of trained neural networks, the more neural nets the more robust is the inference
	# size.nnet is the number of hidden network in the regression. Typically >= the number of parameters in the model
	# transf the vector of transformation for the parameter, ex:transf=c("none","logit","log")
	# bb the vector of bounds for the logit transformation ex:bb=cbind(c(0,0),c(0,2),c(0,0)) (The second column is the only one to be taken into account)
	# If trace=T print messages during the algorithm
	# If rejmethod=F it returns a list with the following components:-
	# $x regression adjusted values
	# $vals - unadjusted values in rejection region (i.e. normal rejection)
	# $wt - the regression weight (i.e. the Epanechnikov weight)
	# $ss - the sumstats corresponding to these points
	# $predmean - estimate of the posterior mean
	if(class(x)=="numeric")
	{
		bb<-cbind(bb)
		x<-cbind(x)
	}

	if(rejmethod)
		transf<-rep("none", dim(x)[2])
	normalise <- function(x,y){
	if(mad(y) == 0)
	return (x)
	else
	return (x/mad(y))
	}

	####Define the weight-decay paramaeter
	repet<-floor(nb.nnet/3)+1
	the_decay<-rep(c(10^(-4),10^(-3),10^(-2)),repet)[1:nb.nnet]
	lt<-dim(x)[2]
	nb_simu<-dim(sumstat)[1]
	for (i in 1:lt)
	{
	if(sum(transf[i] == c("none","log","logit")) == 0){
		stop("transf must be none, log, or logit")
	}
	if(transf[i]=="logit"){
		if(bb[1,i] >= bb[2,i]){
			stop("bounds wrong for logit")
		}
	}
	}

	if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))
	nss <- length(sumstat[1,])
	# scale everything
	    scaled.sumstat <- sumstat
	    for(j in 1:nss){
		scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
	    }
	    target.s.tmp <- target
	    for(j in 1:nss){
		target.s.tmp[,j] <- normalise(target[,j],sumstat[,j][gwt])
	    }

	#déplacé (origine : après calcul de abstol)
		    for (i in 1:lt)
		{
		    if(transf[i] == "log"){
			if(min(x[,i]) <= 0){
				print("log transform: val out of bounds - correcting")
				x.tmp <- ifelse(x[,i] <= 0,max(x[,i]),x[,i])
				x.tmp.min <- min(x.tmp)
				xx[,i] <- ifelse(x[,i] <= 0, x.tmp.min,x[,i])
			}
			x[,i] <- log(x[,i])
		    }
		    else if(transf[i] == "logit"){
			if(min(x[,i]) <= bb[1,i]){
				x.tmp <- ifelse(x[,i] <= bb[1,i],max(x[,i]),x[,i])
				x.tmp.min <- min(x.tmp)
				x[,i] <- ifelse(x[,i] <= bb[1,i], x.tmp.min,x[,i])
			}
			if(max(x[,i]) >= bb[2,i]){
				x.tmp <- ifelse(x[,i] >= bb[2,i],min(x[,i]),x[,i])
				x.tmp.max <- max(x.tmp)
				x[,i] <- ifelse(x[,i] >= bb[2,i], x.tmp.max,x[,i])
			}
			x[,i] <- (x[,i]-bb[1,i])/(bb[2,i]-bb[1,i])
			x[,i] <- log(x[,i]/(1-x[,i]))
		    }
		}

	#boucle le long des targets
	results = NULL
	for(L in 1:nrow(target.s.tmp)){
		print(paste(L, 'th dataset over ', nrow(target.s.tmp), sep=''))
		target.s=as.numeric(target.s.tmp[L,])
		sum1=dst=abstol=wt1=regwt=l1=fit1=ll=predmean=array_pred=mean_averaged=my_residuals=predvar=var_averaged=the_sd=predsd=x.tmp=x.tmp.min=xx=x.tmp.max=NULL
		# calc euclidean distance
		    sum1 <- 0
		    for(j in 1:nss){
			sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
		   }
		   dst <- sqrt(sum1)
		# includes the effect of gwt in the tolerance
		    dst[!gwt] <- floor(max(dst[gwt])+10)

		# wt1 defines the region we're interested in
		    abstol <- quantile(dst,tol)
		    wt1 <- dst <= abstol
		    if(rejmethod){
			regwt <- 1-dst[wt1]^2/abstol^2
			  l1 <- list(x=cbind(x[wt1,]),wt=regwt,ind=wt1,dst=dst)
		    }
		    else{
			  regwt <- 1-dst[wt1]^2/abstol^2
			if(noweight)
				regwt <- rep(1,length(regwt))
			ll<-NULL
			if(trace==TRUE)
				cat("Regression of the mean ")

			for (i in 1:nb.nnet)
			{
				if(trace==TRUE)
					cat(i," ")
				fit1 <- nnet(scaled.sumstat[wt1,],x[wt1,],weights=regwt,decay=the_decay[i],size=size.nnet,linout=T,maxit=500,trace=F, MaxNWts=MaxNWts)
				ll<-c(ll,list(fit1))
			}
			#Compute the residuals
			predmean<-NULL
			array_pred<-array(dim=c(nb.nnet,sum(wt1),lt))
			for (i in 1:nb.nnet)
			{
				array_pred[i,,]<-ll[[i]]$fitted.values
				predmean<-cbind(predmean,as.numeric(predict(ll[[i]],data.frame(rbind(target.s)))))	
			}
			mean_averaged<-NULL
			for (j in 1:lt)
				mean_averaged<-cbind(mean_averaged,apply(array_pred[,,j],FUN=median,MARGIN=2))
			  predmean<-apply(predmean,FUN=median,MARGIN=1)

			  my_residuals<-(x[wt1,]-mean_averaged)

			#Fit a neural network for predicting the conditional variance
			if(trace==TRUE)
				cat("\nRegression of the variance ")
			ll<-NULL
			for (i in 1:nb.nnet)
			{
				if(trace==TRUE)
					cat(i," ")
				fit2 <- nnet(scaled.sumstat[wt1,],log(my_residuals^2),weights=regwt,decay=the_decay[i],size=size.nnet,linout=T,maxit=500,trace=F, MaxNWts=MaxNWts)
				ll<-c(ll,list(fit2))
			}
			if(trace==TRUE)
				cat("\n")
			predvar<-NULL
			array_pred<-array(dim=c(nb.nnet,sum(wt1),lt))
			for (i in 1:nb.nnet){
				array_pred[i,,]<-ll[[i]]$fitted.values
				predvar<-cbind(predvar,as.numeric(predict(ll[[i]],data.frame(rbind(target.s)))))	
			}
			var_averaged<-NULL
			for (j in 1:lt)
				var_averaged<-cbind(var_averaged,apply(array_pred[,,j],FUN=median,MARGIN=2))
				the_sd<-sqrt(exp(var_averaged))
				predsd<-sqrt(exp(apply(predvar,FUN=median,MARGIN=1)))
				res_correc<-sapply(1:lt,FUN=function(i){predmean[i]+ ((predsd[i]*my_residuals[,i])/the_sd[,i]) })
				l1 <- list(x=cbind(res_correc),vals=cbind(x[wt1,]),wt=regwt,ss=sumstat[wt1,],predmean=predmean)
			}
		for (i in 1:lt){
			if(transf[i] == "log"){
				l1$x[,i] <- exp(l1$x[,i])
				l1$vals[,i] <- exp(l1$vals[,i])
			}
			if(transf[i] == "logit"){
				l1$x[,i] <- exp(l1$x[,i])/(1+exp(l1$x[,i]))
				l1$x[,i] <- l1$x[,i]*(bb[2,i]-bb[1,i])+bb[1,i]
				l1$vals[,i] <- exp(l1$vals[,i])/(1+exp(l1$vals[,i]))
				l1$vals[,i] <- l1$vals[,i]*(bb[2,i]-bb[1,i])+bb[1,i]
			}
		}
		results = rbind(results, apply(l1$x, MARGIN=2, FUN="median"))
#		    return(l1)
	#	    write.table(l1$x, col.names=F, row.names=F, file=paste(output,L,sep=""))
		}
		return(results)
}

# function to plot the prior and posterior
babar<-function(a,b,space=2,breaks="auto",AL=0.5,nameA="A",nameB="B",xl="",yl="",mn="",legx="topright", legende=TRUE){ 
       aprime=a;
       bprime=b;
       if(length(a)>length(b)){ bprime=b; aprime=sample(a,length(b),replace=F) }
       if(length(a)<length(b)){ aprime=a; bprime=sample(b,length(a),replace=F) }

       if(breaks=="auto"){
            bks=hist(c(aprime,bprime),plot=F)$breaks
            bklong=space*length(bks)
            bks=hist(c(aprime,bprime),plot=F,breaks=bklong)$breaks
       }
       else{
            bks=breaks
       }

       h1=hist(a,breaks=bks,plot=F)
       h2=hist(b,breaks=bks,plot=F)
       w1=sum(h1$density)
       w2=sum(h2$density)
       d1=max(h1$density)/w1
       d2=max(h2$density)/w2
       d=max(d1,d2)
       par(lwd=1)
       x=barplot(h1$density/w1,col="white",border=par("fg"),ylim=c(0,d),width=.8,space=.2,ylab=yl,xlab=xl,main=mn, cex.lab=1.2) 
       y=c(x,x[length(x)]+.96)-.5
       axis(side=1,at=y,labels=h1$breaks)
       par(lwd=2,lty=1)
       barplot(h2$density/w2,col=rgb(red=.25,blue=.25,green=.25,alpha=AL),border=NA,ylim=c(0,d/w2),width=.8,space=.2,add=T) 
       if(legende==TRUE){legend(legx,legend=c(nameA,nameB),fill=c("white",rgb(red=.25,blue=.25,green=.25,alpha=.5)),cex=1.5,bty="n")}
}


## get the arguments
#for(i in commandArgs()){
#	tmp = strsplit(i, '=')
#	if(tmp[[1]][1] == 'nameA'){ nameA = tmp[[1]][2] }
#	if(tmp[[1]][1] == 'nameB'){ nameB = tmp[[1]][2] }
#	if(tmp[[1]][1] == 'nCPU'){ nCPU = as.integer(tmp[[1]][2]) }
#	if(tmp[[1]][1] == 'model'){ model = tmp[[1]][2] } # model to simulate
#	if(tmp[[1]][1] == 'nMin'){ nMin = as.integer(tmp[[1]][2]) } # minimal number of sequences
#}


get_posterior<-function(nameA, nSubdir, sub_dir_sim, model, sub_dir_model, nPosterior, start, end){
	library(data.table)
	options(digits=5)
	###################
	# get observed data
	# observed data
	coul = c('#ffffcc', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84')
	coul = colorRampPalette(coul)

	
	# get the number of statistics, number of simulations and number of parameters
	tmp = read.table(paste(timeStamp, '/', sub_dir_sim, '/', model, '_0/ABCstat.txt', sep=''), h=T)[,-1]
	nSimulations = nrow(tmp)
	nStats = ncol(tmp)
	ss_sim_tmp = matrix(NA, nrow=nSimulations*nSubdir, ncol=nStats)
	colnames(ss_sim_tmp) = colnames(tmp)
	
	tmp = read.table(paste(timeStamp, '/', sub_dir_sim, '/', model, '_0/priorfile.txt', sep=''), h=T)
	nParams = ncol(tmp)
	params = matrix(NA, nrow=nSimulations*nSubdir, ncol=nParams)
	colnames(params) = colnames(tmp)
	
	tmp = read.table(paste(timeStamp, '/', sub_dir_sim, '/', model, '_0/ABCjsfs.txt', sep=''), h=T)[, -1]
	nBins = ncol(tmp)
	sfs_sim_tmp = matrix(NA, nrow=nSimulations*nSubdir, ncol=nBins)
	colnames(sfs_sim_tmp) = colnames(tmp)
	
	for(rep in seq(0, nSubdir-1, 1)){
		# statistics
		#tmp_ss = read.table(paste(timeStamp, '/', sub_dir_sim, '/', model, '_', rep, '/ABCstat.txt', sep=''), h=T)
		tmp_ss = as.matrix(fread(paste(timeStamp, '/', sub_dir_sim, '/', model, '_', rep, '/ABCstat.txt', sep=''), h=T))[,-1]
		ss_sim_tmp[(rep*nSimulations+1):((rep+1)*nSimulations),] = as.matrix(tmp_ss)
		
		# sfs
		tmp_sfs = as.matrix(fread(paste(timeStamp, '/', sub_dir_sim, '/', model, '_', rep, '/ABCjsfs.txt', sep=''), h=T))[, -1]
		sfs_sim_tmp[(rep*nSimulations+1):((rep+1)*nSimulations),] = as.matrix(tmp_sfs)

		# params
		#tmp_params = read.table(paste(timeStamp, '/', sub_dir_sim, '/', model, '_', rep, '/priorfile.txt', sep=''), h=T)
		tmp_params = as.matrix(fread(paste(timeStamp, '/', sub_dir_sim, '/', model, '_', rep, '/priorfile.txt', sep=''), h=T))
		params[(rep*nSimulations+1):((rep+1)*nSimulations),] = tmp_params
	}
	
	# tested data_set
	obs = start:end
	sim = (1:nrow(ss_sim_tmp))[-obs]
		
	# statistics
	ss_sim = cbind(ss_sim_tmp, sfs_sim_tmp)[sim,]  # with SFS
	ss_obs =  cbind(ss_sim_tmp, sfs_sim_tmp)[obs,]  # with SFS
	
	# parameters
	params_sim = params[sim,]
	params_obs = params[obs,]
	
	##############
	# inferences
	target_rf = data.frame(ss_obs)

	# RANDOM FOREST	
	library(abcrf)
	sim_training = 1:5000 # in case of debug
	params_model_rf = as.matrix(params_sim[sim_training,]) # in case of debug
	stats_model_rf = ss_sim[sim_training,] # in case of debug

	res_rf = list()
	for(i in 1:nParams){
		parameter = params_model_rf[,i]
		data = data.frame(parameter, stats_model_rf)
		mod = regAbcrf(parameter~., data, ntree=1000)
		estimate = predict(mod, target_rf, data)

		param_name = colnames(params_sim)[i]
		res_rf[[param_name]] = list()
		res_rf[[param_name]][['expectation']] = estimate$expectation
		res_rf[[param_name]][['variance']] = estimate$variance
		res_rf[[param_name]][['quantile025']] = estimate$quantiles[1]
		res_rf[[param_name]][['quantile975']] = estimate$quantiles[2] 
	}	
		

	# NEURAL NETWORK
	library('nnet')
	
	target = target_rf
	sumstat = ss_sim

	toRemove=NULL
	for(i in 1:ncol(sumstat)){
		if(sd(sumstat[,i])==0){
			toRemove = c(toRemove, i)
		}
	}
	
	if(is.null(toRemove)==FALSE){
		target = target[, -toRemove]
		sumstat = sumstat[, -toRemove]
	}
		
	x = matrix(as.numeric(unlist(params_sim)), byrow=F, ncol=ncol(params_sim))
	transf_obs = rep("logit", ncol(params_sim))
	bb = rbind(apply(x, MARGIN=2, FUN="min"), apply(x, MARGIN=2, FUN="max"))
	#res2 = abc_nnet_multivar(target=target, x=x, sumstat=sumstat, tol=1000/nrow(x), rejmethod=F, noweight=F, transf=transf_obs, bb=bb, nb.nnet=2*ncol(x), size.nnet=10*ncol(x), trace=T)
	res_nnet = abc_nnet_multivar(target=target, x=x, sumstat=sumstat, tol=nPosterior/nrow(x), rejmethod=F, noweight=F, transf=transf_obs, bb=bb, nb.nnet=2*ncol(x), size.nnet=2*ncol(x), trace=T)

	colnames(res_nnet) = colnames(params_sim)
	res_tot = list()
	res_tot[['random_forest']] = res_rf
	res_tot[['neural_network']] = res_nnet
	res_tot[['real_values']] = params_obs
	# retur inferences	
	return(res_tot)
}

