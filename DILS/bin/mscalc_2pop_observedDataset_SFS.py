#!/shared/mfs/data/software/miniconda/envs/pypy-2.7-5.10.0/bin/pypy
# #!/usr/local/bin/pypy
import os
import sys
from random import sample
project_name = sys.argv[1] # the name of the directory containing the project 
outgroup = int(sys.argv[2]) # if 0: no outgroup, and so, no SFS. If 1: outgroup, and so, SFS
#print("outgroup is {0}".format(outgroup))
def cr_sqrt(x):
	# returns the square root of a variable x
	if x == 0.0:
		return 0.0
	else:
		M = 1.0
		xN = x 
		while xN >= 2.0:
			xN = 0.25*xN
			M = 2.0*M
		while xN<0.5:
			xN = 4.0*xN
			M = 0.5*M
		A = xN
		B = 1.0-xN
		while 1==1:
			A = A*(1.0+0.5*B)
			B = 0.25*(3.0+B)*B*B
			if B<1.0E-15:
				return A*M

def cr_mean(x):
	# returns the mean of a list
	nElement = len(x)
	if nElement == 0:
		return(0)
	else:
		return(sum(x)/(1.0 * nElement))

def cr_std(x, exp_X):
	# returns the standard variation of a list
	nElement = len(x)
	if nElement == 0:
		return(0)
	else:
		if sum(x) == 0:
			return(0)
		else:
			A = sum( [ i**2 for i in x ] )
			A = A/(1.0 * nElement)
			return(cr_sqrt(A-exp_X**2))


def cr_pearsonR(x, y):
	# computes correlation between arrays x and y
	sumXi = 0.0
	sumYi = 0.0
	sumSquareX = 0.0
	sumSquareY = 0.0
	sumXiYi = 0.0
	
	nX = len(x)
	
	for i in range(nX):
		sumXi += x[i];
		sumYi += y[i];
		sumSquareX += x[i]*x[i];
		sumSquareY += y[i]*y[i];
		sumXiYi += x[i]*y[i];
	
	numerator = nX*sumXiYi - sumXi * sumYi
	denom1 = cr_sqrt(nX*sumSquareX - sumXi*sumXi)
	denom2 = cr_sqrt(nX*sumSquareY - sumYi*sumYi)
	if denom1 == 0 or denom2 == 0:
		return(0)
	else:
		return(numerator/(denom1*denom2))


def compFreq(sequences, segsites):
	# returns derived allele frequency for a number of 'segsites' positions
	nDerAll = []
	nInd = len(sequences)
	nPair = nInd*(nInd-1)/2.0
	for i in range(segsites):
		nDerAll.append(0)
		for j in sequences:
			if j[i] == "1":
				nDerAll[i] += 1.0
	pi = [ i * (nInd-i) / nPair for i in nDerAll ]
	freq = [ i/(1.0 * nInd) for i in nDerAll ]
	res = {}
	res['nDer'] = nDerAll
	res['pi_SNPs'] = pi
	res['pi'] = sum(pi)
	res['freq'] = freq
#			if i == 0:		# uncomment to print sequence j #1
#				print(j)	# uncomment to print sequence j #2
	return(res)

def piTot(nDerA, nDerB, nSamA, nSamB, segsites):
	# returns pi tot for the pooled populations from two vectors of allele count per locus
	piT = []
	nTot = nSamA + nSamB
	for i in range(segsites):
		tmp = nDerA[i] + nDerB[i]
		tmp = (nTot - tmp) * tmp / (nTot*(nTot-1)/2.0) # "n ancesral" x "n derived" / C(2,k)
		piT.append(tmp)
	return(piT)


def Fst(piA, piB, piT, segsites):
	# returns Fst as: 1 - mean_pi_over_populations / pi_in_the_whole_alignment
	if piT==0:
		res = 0.0
	else:
		res = 1.0 - (piA+piB)/(2.0*piT)
	return(res)


def sites(freqA, freqB, segsites):
	# test whether the SNP is a polymorphism specific to species A (sxA), species B (sxB), is found in both species (Ss) or differentialy fixed among species (Sf)
	sxA, sxB, ss, sf = 0, 0, 0, 0
	successive_ss = [0]
	previous_site = ""
	for i in range(segsites):
		if freqA[i] == 0:
			if freqB[i] == 1:
				sf += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sf"
			if freqB[i] < 1:
				sxB += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxB"
			continue
		if freqA[i] == 1:
			if freqB[i] == 0:
				sf += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sf"
			if freqB[i] > 0:
				sxB += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxB"
			continue
		else:
			if freqB[i] == 0 or freqB[i] == 1:
				sxA += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxA"
			else:
				ss += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss[len(successive_ss)-1] += 1
				previous_site = "ss"
			continue
	res = {'sxA':sxA, 'sxB':sxB, 'sf':sf, 'ss':ss, 'successive_ss':max(successive_ss)}
	return(res)


def tajD(pi, theta, n, S, a1, a2):
	# returns Tajima D
	# pi = pi for a locus
	# theta = theta for a locus
	# n = sample size
	# S = number of polymorphic positions
	# a1 = sum(1/i)
	# a2 = sum(1/i**2)
	b1 = (n + 1.0) / (3*(n-1))
	b2 = 2.0 * (n**2 + n + 3) / (9*n*(n-1.0))
	c1 = b1 - 1.0/a1
	c2 = b2 - (n + 2.0)/(a1*n) + a2/(a1**2)
	e1 = c1/a1
	e2 = c2/(a1**2 + a2)
	denom = cr_sqrt(e1*S + e2*S*(S-1))
	if denom != 0:
		D = (pi - theta) / denom
	else:
		D = 0.0
	return(D)

def compDiv(spA, spB, segsites):
	div = [] # vector of divergence between A and B
	nPair = 0	
	for i in spA:
		for j in spB:
			nPair += 1
			div.append(0)
			for k in range(segsites):
				if i[k]!=j[k]:
					div[nPair - 1] += 1
	res = {}
	res['divAB'] = cr_mean(div)
	res['minDivAB'] = min(div)
	res['maxDivAB'] = max(div)
	if res['divAB'] > 0:
		res['Gmin'] = res['minDivAB']/res['divAB']
		res['Gmax'] = res['maxDivAB']/res['divAB']
	else:
		res['Gmin'] = 0.0
		res['Gmax'] = 0.0
	return(res)	


def nSs_nSf(ss, sf):
	# returns the number of loci with ss and sf, with ss but no sf, etc ...
	nLoci = len(ss)
	ss_sf, noSs_sf, ss_noSf, noSs_noSf = 0, 0, 0, 0
	for i in range(nLoci):
		if ss[i] == 0:
			if sf[i] == 0:
				noSs_noSf += 1
			if sf[i] > 0:
				noSs_sf += 1
		else:
			if sf[i] == 0:
				ss_noSf += 1
			if sf[i] > 0:
				ss_sf += 1
	res = {'ss_sf': ss_sf, 'ss_noSf': ss_noSf, 'noSs_sf': noSs_sf, 'noSs_noSf': noSs_noSf}
	return(res)


# read information about loci from the bpfile 
if os.path.isfile('{0}/bpfile'.format(project_name)) == False:
	sys.exit("\n\t\033[1;31;40mERROR: bpfile was not found\n\033[0m\n")

infile = open('{0}/bpfile'.format(project_name), "r")
bpfile_header = infile.readline() # first empty line
L = [ float(i) for i in infile.readline().strip().split("\t") ]
nSamA = [ int(i) for i in infile.readline().strip().split("\t") ]
nSamB = [ int(i) for i in infile.readline().strip().split("\t") ]
theta = [ float(i) for i in infile.readline().strip().split("\t") ]
rho = [ float(i) for i in infile.readline().strip().split("\t") ]
infile.close()

nLoci = len(L)
threshold_sim = 1000
if nLoci < threshold_sim:
	# if less than 1000 of loci
	nLoci_sim = nLoci
else:
	# if more than 1000 of loci, then subsample 1000 for simulations only
	nLoci_sim = threshold_sim
	sampled_loci = sample(range(nLoci), threshold_sim)

a1_spA, a1_spB, a2_spA, a2_spB= [], [], [], []
for nsam in nSamA:
	a1_spA.append(sum([ 1.0/i for i in range(1, nsam) ]))
	a2_spA.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))
for nsam in nSamB:
	a1_spB.append(sum([ 1.0/i for i in range(1, nsam) ]))
	a2_spB.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))


# jSFS
nBins_spA = max(nSamA)+1
nBins_spB = max(nSamB)+1
sfs = []
for i in range(nBins_spA):
	sfs.append([0]*nBins_spB)


# ms' output file
outfile_jsfs = open("{0}/ABCjsfs.txt".format(project_name), "w")
header = ''
for i in range(nBins_spA):
	for j in range(nBins_spB):
		header += 'fA{0}_fB{1}\t'.format(i, j)
header = header.strip() + '\n'

outfile_jsfs.write(header)

outfile = open("{0}/ABCstat_global.txt".format(project_name), "w")
outfile_loci = open("{0}/ABCstat_loci.txt".format(project_name), "w")

res = "dataset\tbialsites_avg\tbialsites_std\t"
res += "sf_avg\tsf_std\t"
res += "sxA_avg\tsxA_std\t"
res += "sxB_avg\tsxB_std\t"
res += "ss_avg\tss_std\t"
res += "successive_ss_avg\tsuccessive_ss_std\t"
res += "piA_avg\tpiA_std\t"
res += "piB_avg\tpiB_std\t"
res += "pearson_r_pi\t"
res += "thetaA_avg\tthetaA_std\t"
res += "thetaB_avg\tthetaB_std\t"
res += "pearson_r_theta\t"
res += "DtajA_avg\tDtajA_std\t"
res += "DtajB_avg\tDtajB_std\t"
res += "divAB_avg\tdivAB_std\t"
res += "netdivAB_avg\tnetdivAB_std\t"
res += "minDivAB_avg\tminDivAB_std\t"
res += "maxDivAB_avg\tmaxDivAB_std\t"
res += "Gmin_avg\tGmin_std\t"
res += "Gmax_avg\tGmax_std\t"
res += "FST_avg\tFST_std\t"
res += "pearson_r_divAB_netDivAB\t"
res += "pearson_r_divAB_FST\t"
res += "pearson_r_netDivAB_FST\t"
res += "ss_sf\t" # number of loci with both ss and sf
res += "ss_noSf\t" # number of loci with ss but no sf
res += "noSs_sf\t" # number of loci without ss but with sf
res += "noSs_noSf\n" # number of loci with no Ss or Sf
outfile.write(res)

res = "dataset\tbialsites_avg\t"
res += "sf_avg\t"
res += "sxA_avg\t"
res += "sxB_avg\t"
res += "ss_avg\t"
res += "successive_ss_avg\t"
res += "piA_avg\t"
res += "piB_avg\t"
res += "thetaA_avg\t"
res += "thetaB_avg\t"
res += "DtajA_avg\t"
res += "DtajB_avg\t"
res += "divAB_avg\t"
res += "netdivAB_avg\t"
res += "minDivAB_avg\t"
res += "maxDivAB_avg\t"
res += "Gmin_avg\t"
res += "Gmax_avg\t"
res += "FST_avg\t"
res += "ss_sf\t" # number of loci with both ss and sf
res += "ss_noSf\t" # number of loci with ss but no sf
res += "noSs_sf\t" # number of loci without ss but with sf
res += "noSs_noSf\n" # number of loci with no Ss or Sf
outfile_loci.write(res)

#infile = open(msfile, "r")

test = 0 
nSim_cnt = 0 # count the number of treated multilocus simulations
nLoci_cnt = 0 # count the number of treated loci within a simulation

dataset = [] # vector containing the loci names
#for line in infile:
for line in sys.stdin: # read the ms's output from the stdin
	line = line.strip()
	if '//' in line:
		dataset.append(line.split('//')[1])
	if "segsites" in line:
		if nLoci_cnt == 0:
			# vector to record the sfs
			sfs = []
			for i in range(nBins_spA):
				sfs.append([0]*nBins_spB)
			# end of the vector recording the sfs
			ss_sf, noSs_sf, ss_noSf, noSs_noSf = 0, 0, 0, 0
			bialsites = []
			sf = []
			sxA = []
			sxB = []
			ss = []
			successive_ss = []
			piA = []
			piB = []
			thetaA = []
			thetaB = []
			DtajA = []
			DtajB = []
			divAB = []
			netdivAB = []
			minDivAB = []
			maxDivAB = []
			Gmin = []
			Gmax = []
			FST = []
		nLoci_cnt += 1
		nSam_cnt = 0 # count the number of treated individuals within a locus
		test = 1
		segsites = int(line.split(":")[1])
		bialsites.append(segsites)
		spA, spB = [], []
		continue
	if test == 1:
		if segsites == 0:
			test = 0
			sf.append(0)
			sxA.append(0)
			sxB.append(0)
			ss.append(0)
			successive_ss.append(0)
			piA.append(0)
			piB.append(0)
			thetaA.append(0)
			thetaB.append(0)
			DtajA.append(0)
			DtajB.append(0)
			divAB.append(0)
			netdivAB.append(0)
			minDivAB.append(0)
			maxDivAB.append(0)
			Gmin.append(0)
			Gmax.append(0)
			FST.append(0)
			noSs_noSf += 1
		if segsites != 0:
			if "positions" not in line and line!="\n":
				nSam_cnt += 1
				if nSam_cnt <= nSamA[nLoci_cnt - 1]:
					spA.append(line.strip())
				if nSam_cnt > nSamA[nLoci_cnt - 1] and nSam_cnt <= (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1]):
					spB.append(line.strip())
				if nSam_cnt == (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1]):
					tmpA = compFreq(spA, segsites)
					freqA = tmpA['freq']
					piA.append(tmpA['pi']/L[nLoci_cnt - 1])
					
					tmpB = compFreq(spB, segsites)
					freqB = tmpB['freq']
					piB.append(tmpB['pi']/L[nLoci_cnt - 1])
					
					tmp = sites(freqA, freqB, segsites) # determines the # of sxA, sxB, ss, sf
					sf.append(tmp['sf']/(1.0*L[nLoci_cnt - 1]))
					sxA.append(tmp['sxA']/(1.0*L[nLoci_cnt - 1]))
					sxB.append(tmp['sxB']/(1.0*L[nLoci_cnt - 1]))
					ss.append(tmp['ss']/(1.0*L[nLoci_cnt - 1]))
					successive_ss.append(tmp['successive_ss'])
					# test if loci has both ss, sf, or only one of two, etc ...
					if tmp['sf'] > 0:
						if tmp['ss'] > 0:
							ss_sf += 1
						if tmp['ss'] == 0:
							noSs_sf += 1
					if tmp['sf'] == 0:
						if tmp['ss'] > 0:
							ss_noSf += 1
						if tmp['ss'] == 0:
							noSs_noSf += 1
					
					# theta = sx + ss
					thetaA_locus = (tmp['sxA']+tmp['ss'])/a1_spA[nLoci_cnt - 1]
					thetaB_locus = (tmp['sxB']+tmp['ss'])/a1_spB[nLoci_cnt - 1]
					thetaA.append(thetaA_locus/L[nLoci_cnt - 1 ])
					thetaB.append(thetaB_locus/L[nLoci_cnt - 1 ])
					
					# Taj
					DtajA.append(tajD(tmpA['pi'], thetaA_locus, nSamA[nLoci_cnt - 1], tmp['sxA']+tmp['ss'], a1_spA[nLoci_cnt-1], a2_spA[nLoci_cnt-1]))
					DtajB.append(tajD(tmpB['pi'], thetaB_locus, nSamB[nLoci_cnt - 1], tmp['sxB']+tmp['ss'], a1_spB[nLoci_cnt-1], a2_spB[nLoci_cnt-1]))
					
					# divAB
					div = compDiv(spA, spB, segsites)
					divAB.append(div['divAB']/L[nLoci_cnt - 1 ])
					netdivAB.append(div['divAB']/L[nLoci_cnt - 1 ] - (tmpA['pi']/L[nLoci_cnt - 1] + tmpB['pi']/L[nLoci_cnt - 1]) / 2.0)
					minDivAB.append(div['minDivAB']/L[nLoci_cnt - 1 ])
					maxDivAB.append(div['maxDivAB']/L[nLoci_cnt - 1 ])
					Gmin.append(div['Gmin'])
					Gmax.append(div['Gmax'])
					
					# Fst
					# vector of piT over segsites
					piT = piTot(tmpA['nDer'], tmpB['nDer'], nSamA[nLoci_cnt - 1], nSamB[nLoci_cnt - 1], segsites)
					
					# mean Fst
					FST.append(Fst(sum(tmpA['pi_SNPs']), sum(tmpB['pi_SNPs']), sum(piT), segsites))
					
					#SFS
					if nLoci < threshold_sim:
						for SNPi in range(len(tmpA['nDer'])):
							nDerA = int(tmpA['nDer'][SNPi]) # number of copies of the derived allele in pop A
							nDerB = int(tmpB['nDer'][SNPi]) # number of copies of the derived allele in pop B
							nAncA = nSamA[nLoci_cnt - 1] - nDerA # number of copies of the ancestral allele in pop A
							nAncB = nSamB[nLoci_cnt - 1] - nDerB # number of copies of the ancestral allele in pop B
							if outgroup == 1:
								# unfolded sfs
								sfs[nDerA][nDerB] += 1 # can use the orientation
							else:
								# folded sfs : use the minor allele frequency
								nDerTot = nDerA + nDerB
								nAncTot = nAncA + nAncB
								if nDerTot < nAncTot: # if the allele labeled 0 is the major one
									sfs[nDerA][nDerB] += 1.0 # can use the orientation
								else:
									if nDerTot > nAncTot: # if the allele labeled 0 is the minor one
										sfs[nAncA][nAncB] += 1.0
									else: # if equal frequencies between alleles 0 and 1
										sfs[nAncA][nAncB] += 0.5 
										sfs[nDerA][nDerB] += 0.5
					else:
						if (nLoci_cnt - 1) in sampled_loci:
							for SNPi in range(len(tmpA['nDer'])):
								nDerA = int(tmpA['nDer'][SNPi]) # number of copies of the derived allele in pop A
								nDerB = int(tmpB['nDer'][SNPi]) # number of copies of the derived allele in pop B
								nAncA = nSamA[nLoci_cnt - 1] - nDerA # number of copies of the ancestral allele in pop A
								nAncB = nSamB[nLoci_cnt - 1] - nDerB # number of copies of the ancestral allele in pop B
								if outgroup == 1:
									# unfolded sfs
									sfs[nDerA][nDerB] += 1 # can use the orientation
								else:
									# folded sfs : use the minor allele frequency
									nDerTot = nDerA + nDerB
									nAncTot = nAncA + nAncB
									if nDerTot < nAncTot: # if the allele labeled 0 is the major one
										sfs[nDerA][nDerB] += 1.0 # can use the orientation
									else:
										if nDerTot > nAncTot: # if the allele labeled 0 is the minor one
											sfs[nAncA][nAncB] += 1.0
										else: # if equal frequencies between alleles 0 and 1
											sfs[nAncA][nAncB] += 0.5 
											sfs[nDerA][nDerB] += 0.5

					
	# compute average and std over of statistics over loci
	if nLoci_cnt != 0 and len(ss) == nLoci:
		test = 0
		nSim_cnt += 1
		nLoci_cnt = 0
		
		if nLoci < threshold_sim:
			# if there are less than **threshold_sim** loci, then use all of them 
			# statistics
			bialsites_avg = cr_mean(bialsites)
			bialsites_std = cr_std(bialsites, bialsites_avg)
			sf_avg = cr_mean(sf)
			sf_std = cr_std(sf, sf_avg)
			sxA_avg = cr_mean(sxA)
			sxA_std = cr_std(sxA, sxA_avg)
			sxB_avg = cr_mean(sxB)
			sxB_std = cr_std(sxB, sxB_avg)
			ss_avg = cr_mean(ss)
			ss_std = cr_std(ss, ss_avg)
			successive_ss_avg = cr_mean(successive_ss)
			successive_ss_std = cr_std(successive_ss, successive_ss_avg)
			piA_avg = cr_mean(piA)
			piA_std = cr_std(piA, piA_avg)
			piB_avg = cr_mean(piB)
			piB_std = cr_std(piB, piB_avg)
			pearson_r_pi = cr_pearsonR(piA, piB)
			thetaA_avg = cr_mean(thetaA)
			thetaA_std = cr_std(thetaA, thetaA_avg)
			thetaB_avg = cr_mean(thetaB)
			thetaB_std = cr_std(thetaB, thetaB_avg)
			pearson_r_theta = cr_pearsonR(thetaA, thetaB)
			DtajA_avg = cr_mean(DtajA)
			DtajA_std = cr_std(DtajA, DtajA_avg)
			DtajB_avg = cr_mean(DtajB)
			DtajB_std = cr_std(DtajB, DtajB_avg)
			divAB_avg = cr_mean(divAB)
			divAB_std = cr_std(divAB, divAB_avg)
			netdivAB_avg = cr_mean(netdivAB)
			netdivAB_std = cr_std(netdivAB, netdivAB_avg)
			minDivAB_avg = cr_mean(minDivAB)
			minDivAB_std = cr_std(minDivAB, minDivAB_avg)
			maxDivAB_avg = cr_mean(maxDivAB)
			maxDivAB_std = cr_std(maxDivAB, maxDivAB_avg)
			Gmin_avg = cr_mean(Gmin)
			Gmin_std = cr_std(Gmin, Gmin_avg)
			Gmax_avg = cr_mean(Gmax)
			Gmax_std = cr_std(Gmax, Gmax_avg)
			FST_avg = cr_mean(FST)
			FST_std = cr_std(FST, FST_avg)
			pearson_r_div_netDiv = cr_pearsonR(divAB, netdivAB)
			pearson_r_div_FST = cr_pearsonR(divAB, FST)
			pearson_r_netDiv_FST = cr_pearsonR(netdivAB, FST)
		else:
			# if there are >= **threshold_sim** loci, then use a sub sampling
			# statistics
			bialsites_avg = cr_mean([ bialsites[sampled] for sampled in sampled_loci ])
			bialsites_std = cr_std([ bialsites[sampled] for sampled in sampled_loci ], bialsites_avg)
			sf_avg = cr_mean([ sf[sampled] for sampled in sampled_loci ])
			sf_std = cr_std([ sf[sampled] for sampled in sampled_loci ], sf_avg)
			sxA_avg = cr_mean([ sxA[sampled] for sampled in sampled_loci ])
			sxA_std = cr_std([ sxA[sampled] for sampled in sampled_loci ], sxA_avg)
			sxB_avg = cr_mean([ sxB[sampled] for sampled in sampled_loci ])
			sxB_std = cr_std([ sxB[sampled] for sampled in sampled_loci ], sxB_avg)
			ss_avg = cr_mean([ ss[sampled] for sampled in sampled_loci ])
			ss_std = cr_std([ ss[sampled] for sampled in sampled_loci ], ss_avg)
			successive_ss_avg = cr_mean([ successive_ss[sampled] for sampled in sampled_loci ])
			successive_ss_std = cr_std([ successive_ss[sampled] for sampled in sampled_loci ], successive_ss_avg)
			piA_avg = cr_mean([ piA[sampled] for sampled in sampled_loci ])
			piA_std = cr_std([ piA[sampled] for sampled in sampled_loci ], piA_avg)
			piB_avg = cr_mean([ piB[sampled] for sampled in sampled_loci ])
			piB_std = cr_std([ piB[sampled] for sampled in sampled_loci ], piB_avg)
			pearson_r_pi = cr_pearsonR([ piA[sampled] for sampled in sampled_loci ], [ piB[sampled] for sampled in sampled_loci ])
			thetaA_avg = cr_mean([ thetaA[sampled] for sampled in sampled_loci ])
			thetaA_std = cr_std([ thetaA[sampled] for sampled in sampled_loci ], thetaA_avg)
			thetaB_avg = cr_mean([ thetaB[sampled] for sampled in sampled_loci ])
			thetaB_std = cr_std([ thetaB[sampled] for sampled in sampled_loci ], thetaB_avg)
			pearson_r_theta = cr_pearsonR([ thetaA[sampled] for sampled in sampled_loci ], [ thetaB[sampled] for sampled in sampled_loci ])
			DtajA_avg = cr_mean([ DtajA[sampled] for sampled in sampled_loci ])
			DtajA_std = cr_std([ DtajA[sampled] for sampled in sampled_loci ], DtajA_avg)
			DtajB_avg = cr_mean([ DtajB[sampled] for sampled in sampled_loci ])
			DtajB_std = cr_std([ DtajB[sampled] for sampled in sampled_loci ], DtajB_avg)
			divAB_avg = cr_mean([ divAB[sampled] for sampled in sampled_loci ])
			divAB_std = cr_std([ divAB[sampled] for sampled in sampled_loci ], divAB_avg)
			netdivAB_avg = cr_mean([ netdivAB[sampled] for sampled in sampled_loci ])
			netdivAB_std = cr_std([ netdivAB[sampled] for sampled in sampled_loci ], netdivAB_avg)
			minDivAB_avg = cr_mean([ minDivAB[sampled] for sampled in sampled_loci ])
			minDivAB_std = cr_std([ minDivAB[sampled] for sampled in sampled_loci ], minDivAB_avg)
			maxDivAB_avg = cr_mean([ maxDivAB[sampled] for sampled in sampled_loci ])
			maxDivAB_std = cr_std([ maxDivAB[sampled] for sampled in sampled_loci ], maxDivAB_avg)
			Gmin_avg = cr_mean([ Gmin[sampled] for sampled in sampled_loci ])
			Gmin_std = cr_std([ Gmin[sampled] for sampled in sampled_loci ], Gmin_avg)
			Gmax_avg = cr_mean([ Gmax[sampled] for sampled in sampled_loci ])
			Gmax_std = cr_std([ Gmax[sampled] for sampled in sampled_loci ], Gmax_avg)
			FST_avg = cr_mean([ FST[sampled] for sampled in sampled_loci ])
			FST_std = cr_std([ FST[sampled] for sampled in sampled_loci ], FST_avg)
			pearson_r_div_netDiv = cr_pearsonR([ divAB[sampled] for sampled in sampled_loci ], [ netdivAB[sampled] for sampled in sampled_loci ])
			pearson_r_div_FST = cr_pearsonR([ divAB[sampled] for sampled in sampled_loci ], [ FST[sampled] for sampled in sampled_loci ])
			pearson_r_netDiv_FST = cr_pearsonR([ netdivAB[sampled] for sampled in sampled_loci ], [ FST[sampled] for sampled in sampled_loci ])
			
			# reduce the size of the bpfile (i.e, the number of simulated loci)
			new_bpfile = bpfile_header
			new_bpfile += '\t'.join([ str(L[sampled]) for sampled in sampled_loci ]) + '\n'
			new_bpfile += '\t'.join([ str(nSamA[sampled]) for sampled in sampled_loci ]) + '\n'
			new_bpfile += '\t'.join([ str(nSamB[sampled]) for sampled in sampled_loci ]) + '\n'
			new_bpfile += '\t'.join([ str(theta[sampled]) for sampled in sampled_loci ]) + '\n'
			new_bpfile += '\t'.join([ str(rho[sampled]) for sampled in sampled_loci ]) + '\n'
			out_new_bpfile = open('{0}/bpfile'.format(project_name), "w")
			out_new_bpfile.write(new_bpfile)
			out_new_bpfile.close()
		
		#print("dataset {0}: {1} loci".format(nSim_cnt-1, len(ss)))
		res = ""
		res += "{0}\t{1:.5f}\t{2:.5f}\t".format(nSim_cnt-1, bialsites_avg, bialsites_std)
		res += "{0:.5f}\t{1:.5f}\t".format(sf_avg, sf_std)
		res += "{0:.5f}\t{1:.5f}\t".format(sxA_avg, sxA_std)
		res += "{0:.5f}\t{1:.5f}\t".format(sxB_avg, sxB_std)
		res += "{0:.5f}\t{1:.5f}\t".format(ss_avg, ss_std)
		res += "{0:.5f}\t{1:.5f}\t".format(successive_ss_avg, successive_ss_std)
		res += "{0:.5f}\t{1:.5f}\t".format(piA_avg, piA_std)
		res += "{0:.5f}\t{1:.5f}\t".format(piB_avg, piB_std)
		res += "{0:.5f}\t".format(pearson_r_pi)
		res += "{0:.5f}\t{1:.5f}\t".format(thetaA_avg, thetaA_std)
		res += "{0:.5f}\t{1:.5f}\t".format(thetaB_avg, thetaB_std)
		res += "{0:.5f}\t".format(pearson_r_theta)
		res += "{0:.5f}\t{1:.5f}\t".format(DtajA_avg, DtajA_std)
		res += "{0:.5f}\t{1:.5f}\t".format(DtajB_avg, DtajB_std)
		res += "{0:.5f}\t{1:.5f}\t".format(divAB_avg, divAB_std)
		res += "{0:.5f}\t{1:.5f}\t".format(netdivAB_avg, netdivAB_std)
		res += "{0:.5f}\t{1:.5f}\t".format(minDivAB_avg, minDivAB_std)
		res += "{0:.5f}\t{1:.5f}\t".format(maxDivAB_avg, maxDivAB_std)
		res += "{0:.5f}\t{1:.5f}\t".format(Gmin_avg, Gmin_std)
		res += "{0:.5f}\t{1:.5f}\t".format(Gmax_avg, Gmax_std)
		res += "{0:.5f}\t{1:.5f}\t".format(FST_avg, FST_std)
		res += "{0:.5f}\t".format(pearson_r_div_netDiv)
		res += "{0:.5f}\t".format(pearson_r_div_FST)
		res += "{0:.5f}\t".format(pearson_r_netDiv_FST)
		#res += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}".format(ss_sf, ss_noSf, noSs_sf, noSs_noSf) # total number of ss_sf, ss_noSf and noSs_sf noSs_noSf loci
		res += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}".format(ss_sf/(1.0*nLoci), ss_noSf/(1.0*nLoci), noSs_sf/(1.0*nLoci), noSs_noSf/(1.0*nLoci)) # proportion of ss_sf, ss_noSf, etc ... loci
		res += "\n"
		outfile.write(res)
		
		# write ABCstat for all loci
		for locus_tmp in range(nLoci):
			bialsites_avg = bialsites[locus_tmp]
			sf_avg = sf[locus_tmp]
			sxA_avg = sxA[locus_tmp]
			sxB_avg = sxB[locus_tmp]
			ss_avg = ss[locus_tmp]
			successive_ss_avg = successive_ss[locus_tmp]
			piA_avg = piA[locus_tmp]
			piB_avg = piB[locus_tmp]
			thetaA_avg = thetaA[locus_tmp]
			thetaB_avg = thetaB[locus_tmp]
			DtajA_avg = DtajA[locus_tmp]
			DtajB_avg = DtajB[locus_tmp]
			divAB_avg = divAB[locus_tmp]
			netdivAB_avg = netdivAB[locus_tmp]
			minDivAB_avg = minDivAB[locus_tmp]
			maxDivAB_avg = maxDivAB[locus_tmp]
			Gmin_avg = Gmin[locus_tmp]
			Gmax_avg = Gmax[locus_tmp]
			FST_avg = FST[locus_tmp]
			ss_sf_tmp = 0
			ss_noSf_tmp = 0
			noSs_sf_tmp = 0
			noSs_noSf_tmp = 0
			if ss_avg>0 and sf_avg>0:
				ss_sf_tmp = 1
			if ss_avg>0 and sf_avg==0:
				ss_noSf_tmp = 1
			if ss_avg==0 and sf_avg>0:
				noSs_sf_tmp = 1
			if ss_avg==0 and sf_avg==0:
				noSs_noSf_tmp = 1
			
			#print("dataset {0}: {1} loci".format(nSim_cnt-1, len(ss)))
			res = ""
			res += "{0}\t{1:.5f}\t".format(dataset[locus_tmp], bialsites_avg)
			res += "{0:.5f}\t".format(sf_avg)
			res += "{0:.5f}\t".format(sxA_avg)
			res += "{0:.5f}\t".format(sxB_avg)
			res += "{0:.5f}\t".format(ss_avg)
			res += "{0:.5f}\t".format(successive_ss_avg)
			res += "{0:.5f}\t".format(piA_avg)
			res += "{0:.5f}\t".format(piB_avg)
			res += "{0:.5f}\t".format(thetaA_avg)
			res += "{0:.5f}\t".format(thetaB_avg)
			res += "{0:.5f}\t".format(DtajA_avg)
			res += "{0:.5f}\t".format(DtajB_avg)
			res += "{0:.5f}\t".format(divAB_avg)
			res += "{0:.5f}\t".format(netdivAB_avg)
			res += "{0:.5f}\t".format(minDivAB_avg)
			res += "{0:.5f}\t".format(maxDivAB_avg)
			res += "{0:.5f}\t".format(Gmin_avg)
			res += "{0:.5f}\t".format(Gmax_avg)
			res += "{0:.5f}\t".format(FST_avg)
			res += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}".format(ss_sf_tmp, ss_noSf_tmp, noSs_sf_tmp, noSs_noSf_tmp) # proportion of ss_sf, ss_noSf, etc ... loci
			res += "\n"
			outfile_loci.write(res)
			
		# write sfs
		vector_sfs = []
		for i_spA in range(nBins_spA):
			for i_spB in range(nBins_spB):
				vector_sfs.append(sfs[i_spA][i_spB])
		vector_sfs = '\t'.join( [ str(fifj) for fifj in vector_sfs ]) + '\n'
		outfile_jsfs.write(vector_sfs)

outfile_loci.close()

outfile_jsfs.close()

infile.close()

outfile.close()

