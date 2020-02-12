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
from random import shuffle
from random import randint
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from random import shuffle

def produceBarriers(nLoci, nBarriers):
	# produces a vector of 0 (non barrier) or 1 (barrier), of size equal to the number of loci
	barriers = [0]*nBarriers + [1]*(nLoci-nBarriers)
	shuffle(barriers)
	return(barriers)

	
help = "\t\033[1;31;40mTakes one model specifier, a number of multilocus simulations and a config.yaml file containing prior boundaries as arguments:\033[0m\n\t\t"
help += "\n\t\t".join(["SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N", "AM_1M_1N", "AM_1M_2N", "AM_2M_1N", "AM_2M_2N", "IM_1M_1N", "IM_1M_2N", "IM_2M_1N", "IM_2M_2N", "SI_1N", "SI_2N"])
help += "\n\n"
help += "\t\033[1;32;40m#SI\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;40m#AM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;40m#PAN\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -eN 0 tbs\n"
help += "\t\033[1;32;40m#IM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs\n"
help += "\t\033[1;32;40m#SC\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs\n\n"
help += "\t\033[1;32;40mExample: ./priorgen.py SC_2M_2N 1000 config_yaml\033[0m\n" # last argument is "beta" or "bimodal"

if len(sys.argv) != 4:
	print(help)
	sys.exit()

# Configuration of the prior distribution
nMultilocus = int(sys.argv[2])

shape_bound = [0.01, 5]
N_bound = [0, 0] # number of diploid individuals in the population
T_bound = [0, 0] # number of generations
M_bound = [0, 0] # 4.N.m , so the number of diploid migrant copies is 2.N.m
config_yaml = open(sys.argv[3], 'r')
for i in config_yaml:
	i = i.strip().split(':')
	if(i[0] == 'N_min'):
		N_bound[0] = float(i[1])
	if(i[0] == 'N_max'):
		N_bound[1] = float(i[1])
	if(i[0] == 'Tsplit_min'):
		T_bound[0] = float(i[1])
	if(i[0] == 'Tsplit_max'):
		T_bound[1] = float(i[1])
	if(i[0] == 'M_min'):
		M_bound[0] = float(i[1])
	if(i[0] == 'M_max'):
		M_bound[1] = float(i[1])
	if(i[0] == 'modeBarrier'): # is equal to "beta" or "bimodal"
		modeBarrier = i[1].replace(" ", "")
config_yaml.close()

# convert parameter values in coalescent units
#Nref = (N_bound[1]+N_bound[0])/2.0
Nref = (0+N_bound[1])/2.0
N_bound[0] /= Nref
N_bound[1] /= Nref
T_bound[0] /= (4*Nref)
T_bound[1] /= (4*Nref)

min_Tsc = T_bound[0]
max_Tsc = 0.2
min_Tam = 0.5

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



if sys.argv[1] == "SC_1M_1N":
	# secondary contact
	# ms tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs
	# nsamtot theta rho L nsamA nsamB M12 M21 N1 N2 Tsc Tsplit Tsplit Na

	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#	Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = min_Tsc, high = max_Tsc * Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], N1[sim], N2[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SC_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#       Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = min_Tsc, high = max_Tsc * Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "SC_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#       Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	
	## factor of variation in M.
	if modeBarrier == "beta":
		shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	else:
		nBarriersM12 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]
		nBarriersM21 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = min_Tsc, high = max_Tsc * Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]
	
	# param monolocus: values that will be read by ms
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
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1[sim], N2[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim]))
	
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SC_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## factor of variation in M.
	if modeBarrier == "beta":
		shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	else:
		nBarriersM12 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]
		nBarriersM21 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = min_Tsc, high = max_Tsc * Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	
	
	# param monolocus: values that will be read by ms
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
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = min_Tam*Tsplit[i], high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = min_Tam*Tsplit[i], high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## number of neutral loci
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = min_Tam*Tsplit[i], high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## factor of variation in M.
	if modeBarrier == "beta":
		shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	else:
		nBarriersM12 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]
		nBarriersM21 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	if modeBarrier == "beta":
		priorfile = "N1\tN2\tNa\tTsplit\tTam\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	else:
		priorfile = "N1\tN2\tNa\tTsplit\tTam\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
	for sim in range(nMultilocus):
		if modeBarrier == "beta":
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tam[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		else:
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7:.5f}\t{8}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tam[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])
		
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
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## factor of variation in M.
	if modeBarrier == "beta":
		shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	else:
		nBarriersM12 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]
		nBarriersM21 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	if modeBarrier == "beta":
		priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	else:
		priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tnBarriersM12\tM21\tnBarriersM21\n"
	for sim in range(nMultilocus):
		if modeBarrier == "beta":
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		else:
			priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8}\t{9:.5f}\t{10}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], nBarriersM12[sim], M21[sim], nBarriersM21[sim])
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
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## bf = factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## factor of variation in M.
	if modeBarrier == "beta":
		shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	else:
		nBarriersM12 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]
		nBarriersM21 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
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
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]

	
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## bf = factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	## factor of variation in M.
	if modeBarrier == "beta":
		shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
		shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	else:
		nBarriersM12 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]
		nBarriersM21 = [ randint(1, int(0.9*nLoci-1)) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
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
			M12_vec = [ M12[sim]*i for i in produceBarriers(nLoci, nBarriersM12[sim]) ]
			M21_vec = [ M21[sim]*i for i in produceBarriers(nLoci, nBarriersM21[sim]) ]

		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SI_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SI_2N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## bf = factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ]
		N2_vec = [ N2[sim]*i/rescale for i in scalar_N ]
		Na_vec = [ Na[sim]*i/rescale for i in scalar_N ]

	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "PAN_1N":
	# PAN : msnsam tbs 10000 -t tbs -r tbs tbs -eN 0 tbs"
	# param multilocus: values that will be printed in priorfile.txt
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\n".format(N1[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], N1[sim]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "PAN_2N":
	# PAN : msnsam tbs 10000 -t tbs -r tbs tbs -eN 0 tbs"
	# param multilocus: values that will be printed in priorfile.txt
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)

	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	# param monolocus: values that will be read by ms
	priorfile = "N1\tshape_N_a\tshape_N_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\n".format(N1[sim], shape_N_a[sim], shape_N_b[sim])
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci) # beta distribution of scalars (size = nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim]) # rescaling factor in order to centered the beta distribution around 1
		N1_vec = [ N1[sim]*i/rescale for i in scalar_N ] # local theta for eahc of the nLoci locus
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], N1_vec[locus]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

