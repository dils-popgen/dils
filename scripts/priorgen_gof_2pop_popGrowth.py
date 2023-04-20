#!/usr/bin/python
# #!/home/roux/python/Python-2.7.14/python
# -*- coding: utf-8 -*-

#################################################################################################################################
#################################################################################################################################
#####														       #####
#####    This file is part of Demographic Inferences with Linked Selection : DILS.					  #####
#####														       #####   
#####    DILS is free software: you can redistribute it and/or modify						       #####
#####    it under the terms of the GNU General Public License as published by					       #####
#####    the Free Software Foundation, either version 3 of the License, or						  #####
#####    (at your option) any later version.										#####
#####														       #####    
#####    DILS is distributed in the hope that it will be useful,							    #####
#####    but WITHOUT ANY WARRANTY; without even the implied warranty of						     #####
#####    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the						      #####
#####    GNU General Public License for more details.								       #####
#####														       #####    
#####    You should have received a copy of the GNU General Public License						  #####
#####    along with DILS.  If not, see <https://www.gnu.org/licenses/>.						     #####
#####														       #####    
#####    Please send bugreports with examples or suggestions to							     #####
#####    camille.roux@univ-lille.fr											 #####
#####														       #####    
#####    Or write a post on https://groups.google.com/forum/#!forum/dils---demographic-inferences-with-linked-selection     #####
#####														       #####
#################################################################################################################################
#################################################################################################################################

import sys
from numpy import median
from numpy.random import randint
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from random import shuffle
from random import randint
from random import choice

def provideTimes(Tsplit, Tsmall, nMultilocus):
	# produces a vector of length "nMultilocus" with times smaller than Tsplit
	res = []
	for i in range(nMultilocus):
		Tsmall_tmp = choice(Tsmall)
		if Tsmall_tmp < Tsplit:
			res.append(Tsmall_tmp)
		else:
			res.append(Tsplit[i] * 0.999)
	return(res)


def randomBeta(posterior, nMultilocus):
	estimate = median(posterior)
	a = 15.0
	b = 15.0
	scalar_tmp = beta(a=a, b=b, size=nMultilocus)
	scalar = [ i/(a/(a+b)) for i in scalar_tmp ]
	res = [ estimate * i for i in scalar ]
	return(res)


def produceBarriers(nLoci, nBarriers):
	# produces a vector of 0 (non barrier) or 1 (barrier), of size equal to the number of loci
	barriers = [0]*nBarriers + [1]*(nLoci-nBarriers)
	shuffle(barriers)
	return(barriers)


help = "\t\033[1;31;40mTakes one model specifier, a number of multilocus simulations and a config.yaml file containing prior boundaries as arguments:\033[0m\n\t\t"
help += "\n\t\t".join(["SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N", "AM_1M_1N", "AM_1M_2N", "AM_2M_1N", "AM_2M_2N", "IM_1M_1N", "IM_1M_2N", "IM_2M_1N", "IM_2M_2N", "SI_1N", "SI_2N", "PAN_1N", "PAN_2N"])
help += "\n\n"
help += "\t\033[1;32;40m#SI\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs -g 1 tbs -g 2 tbs\n"
help += "\t\033[1;32;40m#AM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs -g 1 tbs -g 2 tbs\n"
help += "\t\033[1;32;40m#IM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs -g 1 tbs -g 2 tbs\n"
help += "\t\033[1;32;40m#SC\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs -g 1 tbs -g 2 tbs\n"
help += "\t\033[1;32;40m\nExample: ./priorgen_gof_2pop.py [SI_1N; SI_2N; ...; SC_2M_2N] [number of multilocus simulations] [name of the posterior posterior_file] [barriers are modelled by a 'beta' or 'bimodal' distribution] [values are chosen in the posterior_file by using 'joint' values; 'disjoint' values; or 'randomBeta' centered around the median of each parameter]\033[0m\n\n"
help += "\t\033[1;32;40mExample: ./priorgen_gof_2pop.py SC_2M_2N 1000 posterior_file beta joint\033[0m\n"

if len(sys.argv) != 6:
	print(help)
	sys.exit()

# Configuration of the prior distribution
nMultilocus = int(sys.argv[2])
modeBarrier = sys.argv[4] # bimodal / beta
modePrior = sys.argv[5] # joint / disjoint / randomBeta

# read bpfile
infile = open("bpfile", "r")
tmp = infile.readline()
L = infile.readline().strip().split("\t")
nsamA = infile.readline().strip().split("\t")
nsamB = infile.readline().strip().split("\t")
theta = infile.readline().strip().split("\t")
rho = infile.readline().strip().split("\t")
infile.close()

# number of loci
nLoci = len(L)

# sum of nsamA + nsamB
nsam_tot = [ int(nsamA[i]) + int(nsamB[i]) for i in range(nLoci) ]

# get the posterior
infile = open(sys.argv[3], 'r')
posterior = {}
params_posterior = []
header = infile.readline().strip().split('\t')
for i in header:
	params_posterior.append(i)
	if i not in posterior:
		posterior[i] = []
nLines = 0
for line in infile:
	nLines += 1
	line = line.strip().split('\t')
	cnt=0
	for i in line:
		cnt += 1
		posterior_value = float(i)
		if(posterior_value)<0:
			posterior_value = 0.00001
		posterior[ params_posterior[cnt-1] ].append(posterior_value)
#	for i in line:
#		cnt += 1
#		posterior[ params_posterior[cnt-1] ].append(float(i))
infile.close()

migration_models = ['SC_1M_1N', 'SC_2M_1N', 'SC_1M_2N', 'SC_2M_2N', 'AM_1M_1N', 'AM_2M_1N', 'AM_1M_2N', 'AM_2M_2N', 'IM_1M_1N', 'IM_2M_1N', 'IM_1M_2N', 'IM_2M_2N']

# get the lines of the posterior used for the simulations: vector of length nMultilocus
#used_posterior = randint(cnt, size=nMultilocus)

if modePrior == "joint": # modePrior in {joint; disjoint; randomBeta}, where joint takes the exact joint values of the posterior as a prior; disjoint takes random associations of parameter values from posterior; randombeta simulates a beta distribution around the median of the posterior
	used_posterior = [ randint(0, nLines-1) if nLines>1 else 0 for i in range(nMultilocus) ]
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	if 'PAN' not in sys.argv[1]:
		N2 = [ posterior['N2'][i] for i in used_posterior ]
		Na = [ posterior['Na'][i] for i in used_posterior ]
		Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
		founders2 = [ posterior['founders2'][i] for i in used_posterior ]
		Tdem2 = [ posterior['Tdem2'][i] if posterior['Tdem2'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	Tdem1 = [ posterior['Tdem1'][i] if posterior['Tdem1'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]
	
	if '2N' in sys.argv[1]:
		shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
		shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	if sys.argv[1] in migration_models:
		M12 = [ posterior['M12'][i] for i in used_posterior ]
		M21 = [ posterior['M21'][i] for i in used_posterior ]
		if 'SC' in sys.argv[1]:
			Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]
		if 'AM' in sys.argv[1]:
			Tam = [ posterior['Tam'][i] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]
		if '2M' in sys.argv[1]:
			if modeBarrier == "beta":
				shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
				shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
				shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
				shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]
			else:
				nBarriersM12 = [ int(posterior['nBarriersM12'][i]) for i in used_posterior ]
				nBarriersM21 = [ int(posterior['nBarriersM21'][i]) for i in used_posterior ]
				nBarriersM12 = [ i if i <= nLoci else nLoci for i in nBarriersM12 ]
				nBarriersM21 = [ i if i <= nLoci else nLoci for i in nBarriersM21 ]
else:
	if modePrior == "disjoint":
		used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
		N1 = [ posterior['N1'][i] for i in used_posterior ]
		if 'PAN' not in sys.argv[1]:
			used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
			N2 = [ posterior['N2'][i] for i in used_posterior ]
			used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
			Na = [ posterior['Na'][i] for i in used_posterior ]
			used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
			Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
			used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
			founders2 = [ posterior['founders2'][i] for i in used_posterior ]
			Tdem2 = provideTimes(Tsplit=Tsplit, Tsmall=posterior['Tdem2'], nMultilocus=nMultilocus)
		
		used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
		founders1 = [ posterior['founders1'][i] for i in used_posterior ]	
		Tdem1 = provideTimes(Tsplit=Tsplit, Tsmall=posterior['Tdem1'], nMultilocus=nMultilocus)
		
		if '2N' in sys.argv[1]:
			used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
			shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
			used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
			shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]
		
		if sys.argv[1] in migration_models:
			used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
			M12 = [ posterior['M12'][i] for i in used_posterior ]
			used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
			M21 = [ posterior['M21'][i] for i in used_posterior ]
			if 'SC' in sys.argv[1]:
				Tsc = provideTimes(Tsplit=Tsplit, Tsmall=posterior['Tsc'], nMultilocus=nMultilocus)
			if 'AM' in sys.argv[1]:
				Tam = provideTimes(Tsplit=Tsplit, Tsmall=posterior['Tam'], nMultilocus=nMultilocus)
				#Tam = [ posterior['Tam'][randint(0, cnt-1)] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i]*0.999 for i in range(nMultilocus) ]
			if '2M' in sys.argv[1]:
				if modeBarrier == "beta":
					used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
					shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
					used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
					shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
					used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
					shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
					used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
					shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]
				else:
					used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
					nBarriersM12 = [ int(posterior['nBarriersM12'][i]) for i in used_posterior ]
					used_posterior = [ randint(0, nLines-1) if nLines>1 else 0  for i in range(nMultilocus) ]
					nBarriersM21 = [ int(posterior['nBarriersM21'][i]) for i in used_posterior ]
					nBarriersM12 = [ i if i <= nLoci else nLoci for i in nBarriersM12 ]
					nBarriersM21 = [ i if i <= nLoci else nLoci for i in nBarriersM21 ]
	else: # if modeprior == 'randombeta'
		N1 = randomBeta(posterior['N1'], nMultilocus)
		if 'PAN' not in sys.argv[1]:
			N2 = randomBeta(posterior['N2'], nMultilocus)
			Na = randomBeta(posterior['Na'], nMultilocus)
			Tsplit = randomBeta(posterior['Tsplit'], nMultilocus)
			founders2 = randomBeta(posterior['founders2'], nMultilocus)
			Tdem2 = randomBeta(posterior['Tdem2'], nMultilocus)
			Tdem2 = [ Tdem2[i] if Tdem2[i]<Tsplit[i] else 0.999*Tsplit[i] for i in range(nMultilocus) ]
		founders1 = randomBeta(posterior['founders1'], nMultilocus)
		Tdem1 = randomBeta(posterior['Tdem1'], nMultilocus)
		Tdem1 = [ Tdem1[i] if Tdem1[i]<Tsplit[i] else 0.999*Tsplit[i] for i in range(nMultilocus) ]
		
		if '2N' in sys.argv[1]:
			shape_N_a = randomBeta(posterior['shape_N_a'], nMultilocus)
			shape_N_b = randomBeta(posterior['shape_N_b'], nMultilocus)

		if sys.argv[1] in migration_models:
			M12 = randomBeta(posterior['M12'], nMultilocus)
			M21 = randomBeta(posterior['M21'], nMultilocus)
			if 'SC' in sys.argv[1]:
				Tsc = randomBeta(posterior['Tsc'], nMultilocus)
				Tsc = [ Tsc[i] if Tsc[i]<Tsplit[i] else 0.999*Tsplit[i] for i in range(nMultilocus) ]
			if 'AM' in sys.argv[1]:
				Tam = randomBeta(posterior['Tam'], nMultilocus)
				Tam = [ Tam[i] if Tam[i]<Tsplit[i] else 0.999*Tsplit[i] for i in range(nMultilocus) ]
			if '2M' in sys.argv[1]:
				if modeBarrier == "beta":
					shape_M12_a = randomBeta(posterior['shape_M12_a'], nMultilocus)
					shape_M12_b = randomBeta(posterior['shape_M12_b'], nMultilocus)
					shape_M21_a = randomBeta(posterior['shape_M21_a'], nMultilocus)
					shape_M21_b = randomBeta(posterior['shape_M21_b'], nMultilocus)
				else:
					for post in range(nLines):
						posterior['nBarriersM12'][post] = int(posterior['nBarriersM12'][post])
						posterior['nBarriersM21'][post] = int(posterior['nBarriersM21'][post])
					nBarriersM12 = randomBeta(posterior['nBarriersM12'], nMultilocus)
					nBarriersM21 = randomBeta(posterior['nBarriersM21'], nMultilocus)
					for i in range(nMultilocus):
						nBarriersM12[i] = int(nBarriersM12[i])
						nBarriersM21[i] = int(nBarriersM21[i])
					nBarriersM12 = [ i if i <= nLoci else nLoci for i in nBarriersM12 ]
					nBarriersM21 = [ i if i <= nLoci else nLoci for i in nBarriersM21 ]


if sys.argv[1] == "SC_1M_1N":
	# secondary contact
	# ms tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs
	# nsamtot theta rho L nsamA nsamB M12 M21 N1 N2 Tsc Tsplit Tsplit Na
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_1M_2N":
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_2M_1N":
	# param monolocus: values that will be read by ms
	if modeBarrier == "beta":
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTsc\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	else:
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTsc\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
	for sim in range(nMultilocus):
		if modeBarrier == "beta":
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		else:
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10}\t{11:.5f}\t{12}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tsc[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])

		# vectors of size 'nLoci' containing parameters
		if modeBarrier == "beta":
			scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
			scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
			rescaleM12 = shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
			rescaleM21 = shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
			M12_vec = [ M12[sim] * i / rescaleM12 for i in scalar_M12 ]
			M21_vec = [ M21[sim] * i / rescaleM21 for i in scalar_M21 ]
		else:
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]

		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim]))
	
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SC_2M_2N":
	# param monolocus: values that will be read by ms
	if modeBarrier == "beta":
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	else:
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
	for sim in range(nMultilocus):
		if modeBarrier == "beta":
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		else:
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12}\t{13:.5f}\t{14}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])

		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim]) # to centerize the beta distribution around 1
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]

		# vectors of size 'nLoci' containing parameters
		if modeBarrier == "beta":
			scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
			scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
			rescaleM12 = shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
			rescaleM21 = shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
			M12_vec = [ M12[sim] * i / rescaleM12 for i in scalar_M12 ]
			M21_vec = [ M21[sim] * i / rescaleM21 for i in scalar_M21 ]
		else:
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]

	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_1N":
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_2N":
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_1N":
	# param monolocus: values that will be read by ms
	if modeBarrier == "beta":
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTam\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	else:
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTam\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
	for sim in range(nMultilocus):
		if modeBarrier == "beta":
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tam[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		else:
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10}\t{11:.5f}\t{12}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tam[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])

		# vectors of size 'nLoci' containing parameters
		if modeBarrier == "beta":
			scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
			scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
			rescaleM12 = shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
			rescaleM21 = shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
			M12_vec = [ M12[sim] * i / rescaleM12 for i in scalar_M12 ]
			M21_vec = [ M21[sim] * i / rescaleM21 for i in scalar_M21 ]
		else:
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]

		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_2N":
	# param monolocus: values that will be read by ms
	if modeBarrier == "beta":
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	else:
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
	for sim in range(nMultilocus):
		if modeBarrier == "beta":
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		else:
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12}\t{13:.5f}\t{14}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim]) # to centerize the beta distribution around 1
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]

		if modeBarrier == "beta":
			scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
			scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
			rescaleM12 = shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
			rescaleM21 = shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
			M12_vec = [ M12[sim] * i / rescaleM12 for i in scalar_M12 ]
			M21_vec = [ M21[sim] * i / rescaleM21 for i in scalar_M21 ]
		else:
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_1N":
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_2N":
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]
	
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_1N":
	# param monolocus: values that will be read by ms
	if modeBarrier == "beta":
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	else:
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
	for sim in range(nMultilocus):
		if modeBarrier == "beta":
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		else:
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9}\t{10:.5f}\t{11}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])

		# vectors of size 'nLoci' containing parameters
		if modeBarrier == "beta":
			scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
			scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
			rescaleM12 = shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
			rescaleM21 = shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
			M12_vec = [ M12[sim] * i / rescaleM12 for i in scalar_M12 ]
			M21_vec = [ M21[sim] * i / rescaleM21 for i in scalar_M21 ]
		else:
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_2N":
	# param monolocus: values that will be read by ms
	if modeBarrier == "beta":
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	else:
		priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
	for sim in range(nMultilocus):
		if modeBarrier == "beta":
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		else:
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11}\t{12:.5f}\t{13}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])
		
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim]) # to centerize the beta distribution around 1
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]
		
		# vectors of size 'nLoci' containing parameters
		if modeBarrier == "beta":
			scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
			scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
			rescaleM12 = shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
			rescaleM21 = shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
			M12_vec = [ M12[sim] * i / rescaleM12 for i in scalar_M12 ]
			M21_vec = [ M21[sim] * i / rescaleM21 for i in scalar_M21 ]
		else:
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SI_1N":
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "SI_2N":
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "PAN_1N":
	# PAN : msnsam tbs 10000 -t tbs -r tbs tbs -eN 0 tbs -eN tbs tbs
	# param monolocus: values that will be read by ms
	priorfile = "N1\tfounders1\tTdem1\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\n".format(N1[sim], founders1[sim], Tdem1[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], N1[sim], Tdem1[sim], N1[sim]*founders1[sim]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "PAN_2N":
	# param monolocus: values that will be read by ms
	priorfile = "N1\tfounders1\tTdem1\tshape_N_a\tshape_N_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\n".format(N1[sim], founders1[sim], Tdem1[sim], shape_N_a[sim], shape_N_b[sim])
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim]) # to centerize the beta distribution around 1
		N1_vec = [ i/rescale for i in scalar_N ]
		
		for locus in range(nLoci):
			# PAN : msnsam tbs 10000 -t tbs -r tbs tbs -eN 0 tbs -eN tbs tbs
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(nsam_tot[locus], float(theta[locus])*N1_vec[locus], rho[locus], L[locus], N1[sim], Tdem1[sim], N1[sim]*founders1[sim]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

