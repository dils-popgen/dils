#!/shared/mfs/data/software/miniconda/envs/pypy-2.7-5.10.0/bin/pypy
# #!/usr/local/bin/pypy
# #!/usr/local/bin/pypy

# #!/gepv/home2/croux/bin/pypy

# #!/usr/bin/python


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


# read information about loci from the bpfile 
if os.path.isfile('{0}/bpfile'.format(project_name)) == False:
	sys.exit("\n\t\033[1;31;40mERROR: bpfile was not found\n\033[0m\n")

infile = open('{0}/bpfile'.format(project_name), "r")
bpfile_header = infile.readline() # first empty line
L = [ float(i) for i in infile.readline().strip().split("\t") ]
nSamA = [ int(i) for i in infile.readline().strip().split("\t") ]
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
 
a1_spA, a2_spA = [], []
for nsam in nSamA:
	a1_spA.append(sum([ 1.0/i for i in range(1, nsam) ]))
	a2_spA.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))


# SFS
nBins_spA = max(nSamA)
sfs = []
for i in range(nBins_spA):
	sfs.append(0)


# ms' output file
outfile_jsfs = open("{0}/ABCjsfs.txt".format(project_name), "w")
header = ''
for i in range(nBins_spA):
	header += 'fA{0}\t'.format(i)
header = header.strip() + '\n'

outfile_jsfs.write(header)

outfile = open("{0}/ABCstat_global.txt".format(project_name), "w")
outfile_loci = open("{0}/ABCstat_loci.txt".format(project_name), "w")

res = "dataset\tbialsites_avg\tbialsites_std\t"
res += "piA_avg\tpiA_std\t"
res += "thetaA_avg\tthetaA_std\t"
res += "DtajA_avg\tDtajA_std\n"
outfile.write(res)

res = "dataset\tbialsites_avg\t"
res += "piA_avg\t"
res += "thetaA_avg\t"
res += "DtajA_avg\n"
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
			# vector to recort the sfs
			sfs = []
			for i in range(nBins_spA):
				sfs.append(0)
			# end of the vector recording the sfs
			bialsites = []
			piA = []
			thetaA = []
			DtajA = []
		nLoci_cnt += 1
		nSam_cnt = 0 # count the number of treated individuals within a locus
		test = 1
		segsites = int(line.split(":")[1])
		bialsites.append(segsites)
		spA = []
		continue
	if test == 1:
		if segsites == 0:
			test = 0
			piA.append(0)
			thetaA.append(0)
			DtajA.append(0)
		if segsites != 0:
			if "positions" not in line and line!="\n":
				nSam_cnt += 1
				if nSam_cnt <= nSamA[nLoci_cnt - 1]:
					spA.append(line.strip())
				if nSam_cnt == nSamA[nLoci_cnt - 1]:
					tmpA = compFreq(spA, segsites)
					freqA = tmpA['freq']
					piA.append(tmpA['pi']/L[nLoci_cnt - 1])
					
					thetaA_locus = segsites/a1_spA[nLoci_cnt - 1]
					thetaA.append(thetaA_locus/L[nLoci_cnt - 1 ])
					
					# Taj
					DtajA.append(tajD(tmpA['pi'], thetaA_locus, nSamA[nLoci_cnt - 1], segsites, a1_spA[nLoci_cnt-1], a2_spA[nLoci_cnt-1]))
					
					# SFS
					if nLoci < threshold_sim:
						for SNPi in range(len(tmpA['nDer'])):
							nDerA = int(tmpA['nDer'][SNPi]) # number of copies of the derived allele in pop A
							nAncA = nSamA[nLoci_cnt - 1] - nDerA # number of copies of the ancestral allele in pop A
							if outgroup == 1:
								# unfolded sfs
								sfs[nDerA] += 1 # can use the orientation
							else:
								# folded sfs : use the minor allele frequency
								if nDerA <= nAncA: # if the allele labeled 0 is the major one
									sfs[nDerA] += 1.0 # can use the orientation
								else:
									sfs[nAncA] += 1.0
					else:
						if (nLoci_cnt - 1) in sampled_loci:
							for SNPi in range(len(tmpA['nDer'])):
								nDerA = int(tmpA['nDer'][SNPi]) # number of copies of the derived allele in pop A
								nAncA = nSamA[nLoci_cnt - 1] - nDerA # number of copies of the ancestral allele in pop A
								if outgroup == 1:
									# unfolded sfs
									sfs[nDerA] += 1 # can use the orientation
								else:
									# folded sfs : use the minor allele frequency
									if nDerA <= nAncA: # if the allele labeled 0 is the major one
										sfs[nDerA] += 1.0 # can use the orientation
									else:
										sfs[nAncA] += 1.0
					
	# compute average and std over of statistics over loci
	if nLoci_cnt != 0 and len(piA) == nLoci:
		test = 0
		nSim_cnt += 1
		nLoci_cnt = 0
		if nLoci < threshold_sim:
			# if there are less than **threshold_sim** loci, then use all of them 
			# statistics
			bialsites_avg = cr_mean(bialsites)
			bialsites_std = cr_std(bialsites, bialsites_avg)
			piA_avg = cr_mean(piA)
			piA_std = cr_std(piA, piA_avg)
			thetaA_avg = cr_mean(thetaA)
			thetaA_std = cr_std(thetaA, thetaA_avg)
			DtajA_avg = cr_mean(DtajA)
			DtajA_std = cr_std(DtajA, DtajA_avg)
		else:
			# statistics
			bialsites_avg = cr_mean([ bialsites[sampled] for sampled in sampled_loci ])
			bialsites_std = cr_std([ bialsites[sampled] for sampled in sampled_loci ], bialsites_avg)
			piA_avg = cr_mean([ piA[sampled] for sampled in sampled_loci ])
			piA_std = cr_std([ piA[sampled] for sampled in sampled_loci ], piA_avg)
			thetaA_avg = cr_mean([ thetaA[sampled] for sampled in sampled_loci ])
			thetaA_std = cr_std([ thetaA[sampled] for sampled in sampled_loci ], thetaA_avg)
			DtajA_avg = cr_mean([ DtajA[sampled] for sampled in sampled_loci ])
			DtajA_std = cr_std([ DtajA[sampled] for sampled in sampled_loci ], DtajA_avg)
		
			# reduce the size of the bpfile (i.e, the number of simulated loci)
			new_bpfile = bpfile_header
			new_bpfile += '\t'.join([ str(L[sampled]) for sampled in sampled_loci ]) + '\n'
			new_bpfile += '\t'.join([ str(nSamA[sampled]) for sampled in sampled_loci ]) + '\n'
			new_bpfile += '\t'.join([ str(theta[sampled]) for sampled in sampled_loci ]) + '\n'
			new_bpfile += '\t'.join([ str(rho[sampled]) for sampled in sampled_loci ]) + '\n'
			out_new_bpfile = open('{0}/bpfile'.format(project_name), "w")
			out_new_bpfile.write(new_bpfile)
			out_new_bpfile.close()
		res = ""
		res += "{0}\t{1:.5f}\t{2:.5f}\t".format(nSim_cnt-1, bialsites_avg, bialsites_std)
		res += "{0:.5f}\t{1:.5f}\t".format(piA_avg, piA_std)
		res += "{0:.5f}\t{1:.5f}\t".format(thetaA_avg, thetaA_std)
		res += "{0:.5f}\t{1:.5f}".format(DtajA_avg, DtajA_std)
		res += "\n"
		outfile.write(res)
		
		# write ABCstat for all loci
		for locus_tmp in range(nLoci):
			bialsites_avg = bialsites[locus_tmp]
			piA_avg = piA[locus_tmp]
			thetaA_avg = thetaA[locus_tmp]
			DtajA_avg = DtajA[locus_tmp]
			
			#print("dataset {0}: {1} loci".format(nSim_cnt-1, len(ss)))
			res = ""
			res += "{0}\t{1:.5f}\t".format(dataset[locus_tmp], bialsites_avg)
			res += "{0:.5f}\t".format(piA_avg)
			res += "{0:.5f}\t".format(thetaA_avg)
			res += "{0:.5f}\t".format(DtajA_avg)
			res += "\n"
			outfile_loci.write(res)
			
		# write sfs
		vector_sfs = []
		for i_spA in range(nBins_spA):
			vector_sfs.append(sfs[i_spA])
		vector_sfs = '\t'.join( [ str(fi) for fi in vector_sfs ]) + '\n'
		outfile_jsfs.write(vector_sfs)
		
outfile_loci.close()

outfile_jsfs.close()

infile.close()

outfile.close()




