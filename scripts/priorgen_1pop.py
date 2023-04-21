#!/shared/software/miniconda/envs/python-2.7/bin/python2.7
# #!/home/roux/python/Python-2.7.14/python
# #!/usr/bin/python
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
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from numpy import log
from random import shuffle
help = "\t\033[1;31;40mTakes one model specifier, a number of multilocus simulations and a config.yaml file containing prior boundaries as arguments:\033[0m\n\t\t"
#help += "\n\t\t".join(["Constant_1N", "Constant_2N", "Discrete_1N", "Discrete_2N", "Expo_1N", "Expo_2N"])
help += "\n\t\t".join(["Constant_1N", "Constant_2N", "Expansion_1N", "Expansion_2N", "Contraction_1N", "Contraction_2N"])
help += "\n\n"
help += "\t\033[1;32;40m#Constant\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs\n"
help += "\t\033[1;32;40m#Discrete\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs\n"
help += "\t\033[1;32;40m#Expo\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -G tbs -eG tbs 0.0 -eN tbs tbs\n"

help += "\t\033[1;32;40mExample: ./priorgen_1pop.py Constant_1N 1000 config.yaml\033[0m\n"

if len(sys.argv) != 4:
	print(help)
	sys.exit()

# Configuration of the prior distribution
nMultilocus = int(sys.argv[2])

shape_bound = [1, 5]
N_bound = [0, 0] # number of diploid individuals in the population
T_bound = [0, 10]
config_yaml = open(sys.argv[3], 'r')
for i in config_yaml:
	i = i.strip().split(':')
	if(i[0] == 'N_min'):
		N_bound[0] = float(i[1])
	if(i[0] == 'N_max'):
		N_bound[1] = float(i[1])
	if(i[0] == 'Tchanges_min'):
		T_bound[0] = float(i[1])
	if(i[0] == 'Tchanges_max'):
		T_bound[1] = float(i[1])
config_yaml.close()

# read bpfile
infile = open("bpfile", "r")
tmp = infile.readline()
L = infile.readline().strip().split("\t") # in number of nucleotides
nsamA = infile.readline().strip().split("\t") # number of individuals within the population 
mu = infile.readline().strip().split("\t") # mutation rate per generation and per base pair
rec = infile.readline().strip().split("\t") # recombiation rate per generation and per base pair
infile.close()

# converts in integer/floating values
L = [ float(i) for i in L ]
nsamA = [ int(i) for i in nsamA ]
mu = [ float(i) for i in mu ]
rec = [ float(i) for i in rec ]

# number of loci
nLoci = len(L)


if sys.argv[1] == "Constant_1N":
	# msnsam tbs 10000 -t tbs
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N

	# param monolocus: values that will be read by ms
	priorfile = "N\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\n".format(N[sim])
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]
			rho_locus = rho[locus]*N[sim]
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}".format(nsamA[locus], theta_locus, rho_locus, L[locus]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Constant_2N":
	# msnsam tbs 10000 -t tbs
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N

	## bf = factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N\tshape_N_a\tshape_N_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\n".format(N[sim], shape_N_a[sim], shape_N_b[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
	
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]*scalar_N[locus]/rescale
			rho_locus = rho[locus]*N[sim]*scalar_N[locus]/rescale
			
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}".format(nsamA[locus], theta_locus, rho_locus, L[locus]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Expansion_1N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs
	# nsamA theta 
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N
	
	T = [ uniform(low = 0.1*i, high = 4*i, size = 1)[0] for i in N ] # uncommented for the distributed version
#	T = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus) # explore the whole space of parameters for the paper, to show the limits
	
	Nanc = [ uniform(low = 0, high = 0.5*i, size = 1)[0] for i in N ]# uncommented for the distributed version
#	Nanc = [ uniform(low = 0, high = 1*i, size = 1)[0] for i in N ]#  for the paper
	
	# param monolocus: values that will be read by ms
	priorfile = "N\tNpast\tTdem\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\n".format(N[sim], Nanc[sim], T[sim])
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]
			rho_locus = rho[locus]*N[sim]
			
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(nsamA[locus], theta_locus, rho_locus, L[locus], T[sim]/(4.0*N[sim]), Nanc[sim]/(1.0*N[sim])) )
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Expansion_2N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs
	# nsamA theta 
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N
	
	## bf = factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	
	T = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)# uncommented for the distributed version
#	T = [ uniform(low = 0.1*i, high = 4*i, size = 1)[0] for i in N ]# explore the whole space of parameters for the paper, to show the limits
	Nanc = [ uniform(low = N_bound[0], high = 1*i, size = 1)[0] for i in N ]# uncommented for the distributed version
#	Nanc = [ uniform(low = 0, high = 0.5*i, size = 1)[0] for i in N ]# explore the whole space of parameters for the paper, to show the limits
	
	# param monolocus: values that will be read by ms
	priorfile = "N\tNpast\tshape_N_a\tshape_N_b\tTdem\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\n".format(N[sim], Nanc[sim], shape_N_a[sim], shape_N_b[sim], T[sim])
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]*scalar_N[locus]/rescale
			rho_locus = rho[locus]*N[sim]*scalar_N[locus]/rescale
			
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(nsamA[locus], theta_locus, rho_locus, L[locus], T[sim]/(4.0*N[sim]), Nanc[sim]/(1.0*N[sim])))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Contraction_1N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs
	# nsamA theta 
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N
	
	Nanc = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	
#	N = [ uniform(low = N_bound[0], high = i, size = 1)[0] for i in Nanc ] # explore the whole space of parameters for the paper, to show the limits
	N = [ uniform(low = N_bound[0], high = 0.25*i, size = 1)[0] for i in Nanc ]# uncommented for the distributed version

#	T = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)# explore the whole space of parameters for the paper, to show the limits
	T = [ uniform(low = 0.1*i, high = 4*i, size = 1)[0] for i in N ]# uncommented for the distributed version
	
	
	# param monolocus: values that will be read by ms
	priorfile = "N\tNpast\tTdem\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\n".format(N[sim], Nanc[sim], T[sim])
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]
			rho_locus = rho[locus]*N[sim]
			
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(nsamA[locus], theta_locus, rho_locus, L[locus], T[sim]/(4.0*N[sim]), Nanc[sim]/(1.0*N[sim])))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Contraction_2N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs
	# nsamA theta 
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N
	
	Nanc = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	
	
#	N = [ uniform(low = N_bound[0], high = i, size = 1)[0] for i in Nanc ] # explore the whole space of parameters for the paper, to show the limits
	N = [ uniform(low = N_bound[0], high = 0.25*i, size = 1)[0] for i in Nanc ]# uncommented for the distributed version

#	T = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)# explore the whole space of parameters for the paper, to show the limits
	T = [ uniform(low = 0.1*i, high = 4*i, size = 1)[0] for i in N ]# uncommented for the distributed version
	
	# param monolocus: values that will be read by ms
	priorfile = "N\tNpast\tshape_N_a\tshape_N_b\tTdem\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\n".format(N[sim], Nanc[sim], shape_N_a[sim], shape_N_b[sim], T[sim])
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]*scalar_N[locus]/rescale
			rho_locus = rho[locus]*N[sim]*scalar_N[locus]/rescale
			
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(nsamA[locus], theta_locus, rho_locus, L[locus], T[sim]/(4.0*N[sim]), Nanc[sim]/(1.0*N[sim])))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

