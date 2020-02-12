#!/usr/bin/python

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
import os
import time

# Example to submit jobs using slurm:
#    for model in AM IM SC; do for M in 1M 2M; do for N in 1N 2N; do ./submit.py 100000 10 ${model}_${M}_${N}; done; done; done
#    for model in SI; do for N in 1N 2N; do ./submit.py 100000 10 ${model}_${N}; done; done


if len(sys.argv) != 13:
	print("\n\tsubmit_simulations_2pop_test_monolocus.py [outgroup] [nmultilocus] [iteration] [model: SI_x AM_x IM_x SC_x PSC_x PAM_x] [nameA] [nameB] [sub_dir_sim] [posterior_file] [beta or bimodal] [ joint; disjoint; randomBeta ] [bin folder] [binpath]")
	print("\n\tex: submit_simulations_2pop_test_monolocus.py 1 1000 2 SI_1N flo mal sim_SI_1N posterior_IM_1M_2N.txt Ng4PymB1dy [beta]\n\tto simulate 1000 multilocus simulations at the second iteration, in the folder sim_SI_1N variable") 
	sys.exit(0)
#submit_simulations_2pop_test_monolocus.py 0 10000 1 ${best_model} txn ama locus_modelComp ${PWD}/variableBeta3roundsTxnAma/best_model_5/posterior_bestModel.txt variableBeta3roundsTxnAma beta /shared/mfs/data/home/croux/softwares/ABConline/bin variable
outgroup = int(sys.argv[1])
nmultilocus = int(sys.argv[2]) # 10000
iteration = int(sys.argv[3]) # 2
model = sys.argv[4]
nameA = sys.argv[5] # name of the species A
nameB = sys.argv[6] # name of the species B
sub_dir_sim = sys.argv[7] # name of the subdir where the simulations will be run
posterior_file = sys.argv[8]
timeStamp = sys.argv[9]
modeBarrier = sys.argv[10] # beta / bimodal
binpath = sys.argv[11] 
population_growth = sys.argv[12] # constant / variable

path = os.getcwd() + '/{0}'.format(timeStamp)

test_bpfile = os.path.isfile('{0}/bpfile'.format(path))
if test_bpfile == False:
	sys.exit('\n\tERROR in submit_simulations_2pop.py : the file {0}/bpfile is not found\n'.format(path))
else:
	infile = open('{0}/bpfile'.format(path), 'r')
	tmp = infile.readline()
	tmp = infile.readline().strip().split('\t')
	nlocus = len(tmp)
	infile.close()
mscommand = ""

if population_growth == 'constant':
	if "SC" in model:
		mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs"
	if "IM" in model:
		mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs"
else:
	if "SC" in model:
		mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs"
	if "IM" in model:
		mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs"

#if mscommand == "":
#	print("You specified a wrong model: SI_x, AM_x, AM_x or SC_x\n")
#	sys.exit()

if mscommand == "":
	# if model of isolation
	os.system("touch {0}/locus_modelComp/isolation/ABCstat.txt".format(path))
	os.system("touch {0}/locus_modelComp/migration/ABCstat.txt".format(path))
	sys.exit()
else:
	# migration
	sub_dir_model='migration'
	tmp_migration = "cp {0}/bpfile {0}/locus_modelComp/{1}/; ".format(path, sub_dir_model)
	tmp_migration += "cd {0}/locus_modelComp/{1}/; ".format(path, sub_dir_model)
	tmp_migration += "python {0}/priorgen_gof_2pop_test_monolocus.py {1} {2} {3} {4} migration {5} | {0}/msnsam tbs {6} {7} >tmp.ms ;".format(binpath, model, nmultilocus, posterior_file, modeBarrier, population_growth, nmultilocus*nlocus, mscommand, outgroup)
	tmp_migration += "cat tmp.ms | pypy {0}/mscalc_2pop_SFS.py {1}; ".format(binpath, outgroup)
	print(tmp_migration)
	os.system(tmp_migration)

	# isolation 
	sub_dir_model='isolation'
	tmp_isolation = "cp {0}/bpfile {0}/locus_modelComp/{1}/; ".format(path, sub_dir_model)
	tmp_isolation += "cd {0}/locus_modelComp/{1}/; ".format(path, sub_dir_model)
	tmp_isolation += "python {0}/priorgen_gof_2pop_test_monolocus.py {1} {2} {3} {4} isolation {5} | {0}/msnsam tbs {6} {7} >tmp.ms ;".format(binpath, model, nmultilocus, posterior_file, modeBarrier, population_growth, nmultilocus*nlocus, mscommand, outgroup)
	tmp_isolation += "cat tmp.ms | pypy {0}/mscalc_2pop_SFS.py {1}; ".format(binpath, outgroup)
	print(tmp_isolation)
	os.system(tmp_isolation)

