#!/usr/bin/python
# #!/home/roux/python/Python-2.7.14/python
# -*- coding: utf-8 -*-

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

import sys
from numpy import median
from numpy.random import randint
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from random import shuffle
from random import randint
from random import choice

#Tsc = provideTimes(Tsplit=posterior['Tsplit'], Tsmall=posterior['Tsc'], sizePosterior=cnt, nMultilocus=nMultilocus)
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


help = "\t\033[1;31;40mTakes one model specifier with heterogeneous migration, a number of monolocus simulations and a posterior file containing parameter estimates:\033[0m\n\t\t"
help += "\n\t\t".join(["SC_2M_1N", "SC_2M_2N", "IM_2M_1N", "IM_2M_2N"])
help += "\n\n"
help += "\t\033[1;32;40m#IM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;40m#SC\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;40m\nExample: ./priorgen_gof_2pop_test_monolocus.py [SC_2M_1N, SC_2M_2N, IM_2M_1N, IM_2M_2N] [number of multilocus simulations] [name of the posterior posterior_file] [barriers are modelled by a 'beta' or 'bimodal' distribution] ['migration' or 'isolation'] [population_growth = 'constant' or 'variable']\033[0m\n\n"
help += "\t\033[1;32;40mExample: ./priorgen_gof_2pop_test_monolocus.py SC_2M_2N 10000 posterior_file migration population_growth\033[0m\n"

if len(sys.argv) != 7:
	print(help)
	sys.exit()

# Configuration of the prior distribution
nMultilocus = int(sys.argv[2])
modeBarrier = sys.argv[4] # bimodal / beta
if sys.argv[5]=='isolation': # a locus can be "migration" (unlinked to barrier, with M>0) or "isolation" (linked to a barrier, with M=0)
	locus_category = 0
else:
	locus_category = 1

population_growth = sys.argv[6]

# read bpfile
infile = open("bpfile", "r")
tmp = infile.readline()
L = infile.readline().strip().split("\t")
nsamA = infile.readline().strip().split("\t")
nsamB = infile.readline().strip().split("\t")
theta = infile.readline().strip().split("\t")
rho = infile.readline().strip().split("\t")
infile.close()

L = median( [ float(i) for i in L ] )
nsamA = median( [ float(i) for i in nsamA ] )
nsamB = median( [ float(i) for i in nsamB ] )
theta = median( [ float(i) for i in theta ] )
rho = median( [ float(i) for i in rho ] )

# number of loci
nLoci = 1

# sum of nsamA + nsamB
L = [L]
nsam_tot = [nsamA + nsamB]
nsamA = [nsamA]
nsamB = [nsamB]
theta = [theta]
rho = [rho]

# rewrite the bpfile
tmp_bpfile = '# tmp\n'
tmp_bpfile += str(float(L[0])) + '\n'
tmp_bpfile += str(int(nsamA[0])) + '\n'
tmp_bpfile += str(int(nsamB[0])) + '\n'
tmp_bpfile += str(theta[0]) + '\n'
tmp_bpfile += str(rho[0])

outfile_bpfile = open("bpfile", "w")
outfile_bpfile.write(tmp_bpfile)
outfile_bpfile.close()

# get the posterior
infile = open(sys.argv[3], 'r')
posterior = {}
params_posterior = []
header = infile.readline().strip().split('\t')
nLines = 0
for i in header:
	params_posterior.append(i)
	if i not in posterior:
		posterior[i] = []
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
infile.close()

migration_models = ['SC_2M_1N', 'SC_2M_2N', 'IM_2M_1N', 'IM_2M_2N']

# get the lines of the posterior used for the simulations: vector of length nMultilocus
# used_posterior = randint(cnt, size=nMultilocus)

used_posterior = [ randint(0, nLines-1) for i in range(nMultilocus) ]
N1 = [ posterior['N1'][i] for i in used_posterior ]
N2 = [ posterior['N2'][i] for i in used_posterior ]
Na = [ posterior['Na'][i] for i in used_posterior ]
Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

if population_growth == 'variable':
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]
	Tdem1 = [ posterior['Tdem1'][i] if posterior['Tdem1'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] if posterior['Tdem2'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

if '2N' in sys.argv[1]:
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

if sys.argv[1] in migration_models:
	M12 = [ posterior['M12'][i]*locus_category for i in used_posterior ] # locus_category = 1 for unlinked loci to a barrier, locus_category = 0 if linked to a barrier
	M21 = [ posterior['M21'][i]*locus_category for i in used_posterior ]
	if 'SC' in sys.argv[1]:
		Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]
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


if population_growth == 'constant':
	# print the prior for simulations
	if sys.argv[1] == "SC_2M_1N":
		if modeBarrier == "beta":
			priorfile = "N1\tN2\tNa\tTsplit\tTsc\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
		else:
			priorfile = "N1\tN2\tNa\tTsplit\tTsc\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
		for sim in range(nMultilocus):
			if modeBarrier == "beta":
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
			else:
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7:.5f}\t{8}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tsc[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])
			# vectors of size 'nLoci' containing parameters
			if modeBarrier == "beta":
				scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
				scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
				rescaleM12 = shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
				rescaleM21 = shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
				M12_vec = [ M12[sim] * i / rescaleM12 for i in scalar_M12 ]
				M21_vec = [ M21[sim] * i / rescaleM21 for i in scalar_M21 ]
			else:
				M12_vec = [ M12[sim]*1 for i in range(nLoci) ]
				M21_vec = [ M21[sim]*1 for i in range(nLoci) ]
			for locus in range(nLoci):
				print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1[sim], N2[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim]))
		
		outfile = open("priorfile.txt", "w")
		outfile.write(priorfile)
		outfile.close()

	if sys.argv[1] == "SC_2M_2N":
		if modeBarrier == "beta":
			priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
		else:
			priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
		for sim in range(nMultilocus):
			if modeBarrier == "beta":
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
			else:
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8}\t{9:.5f}\t{10}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])
			# vectors of size 'nLoci' containing parameters
			scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
			rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
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
				M12_vec = [ M12[sim]*1 for i in range(nLoci) ]
				M21_vec = [ M21[sim]*1 for i in range(nLoci) ]
			for locus in range(nLoci):
				print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

		outfile = open("priorfile.txt", "w")
		outfile.write(priorfile)
		outfile.close()

	if sys.argv[1] == "IM_2M_1N":
		if modeBarrier == "beta":
			priorfile = "N1\tN2\tNa\tTsplit\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
		else:
			priorfile = "N1\tN2\tNa\tTsplit\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
		for sim in range(nMultilocus):
			if modeBarrier == "beta":
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
			else:
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5}\t{6:.5f}\t{7}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])
		
			# vectors of size 'nLoci' containing parameters
			if modeBarrier == "beta":
				scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
				scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
				rescaleM12 = shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
				rescaleM21 = shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
				M12_vec = [ M12[sim] * i / rescaleM12 for i in scalar_M12 ]
				M21_vec = [ M21[sim] * i / rescaleM21 for i in scalar_M21 ]
			else:
				M12_vec = [ M12[sim]*1 for i in range(nLoci) ]
				M21_vec = [ M21[sim]*1 for i in range(nLoci) ]

			for locus in range(nLoci):
				# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
				print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim]))

		outfile = open("priorfile.txt", "w")
		outfile.write(priorfile)
		outfile.close()


	if sys.argv[1] == "IM_2M_2N":
		if modeBarrier == "beta":
			priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
		else:
			priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
		for sim in range(nMultilocus):
			if modeBarrier == "beta":
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
			else:
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7}\t{8:.5f}\t{9}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])
			# vectors of size 'nLoci' containing parameters
			scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
			rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
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
				M12_vec = [ M12[sim]*1 for i in range(nLoci) ]
				M21_vec = [ M21[sim]*1 for i in range(nLoci) ]
			for locus in range(nLoci):
				# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
				print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
		outfile = open("priorfile.txt", "w")
		outfile.write(priorfile)
		outfile.close()
else:
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
				M12_vec = [ M12[sim]*1 for i in range(nLoci) ]
				M21_vec = [ M21[sim]*1 for i in range(nLoci) ]

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
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
			else:
				priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12}\t{12:.5f}\t{13}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])

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
				M12_vec = [ M12[sim]*1 for i in range(nLoci) ]
				M21_vec = [ M21[sim]*1 for i in range(nLoci) ]

		
			for locus in range(nLoci):
				print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

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
				M12_vec = [ M12[sim]*1 for i in range(nLoci) ]
				M21_vec = [ M21[sim]*1 for i in range(nLoci) ]
			
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
				M12_vec = [ M12[sim]*1 for i in range(nLoci) ]
				M21_vec = [ M21[sim]*1 for i in range(nLoci) ]
			
			for locus in range(nLoci):
				#print(nsam_tot[locus])
				print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

		outfile = open("priorfile.txt", "w")
		outfile.write(priorfile)
		outfile.close()

