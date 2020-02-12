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


if len(sys.argv) != 14:
	print("\n\t submit_simulations_gof_2pop_popGrowth.py [nmultilocus] [iteration] [model: SI_x AM_x IM_x SC_x PSC_x PAM_x] [nameA] [nameB] [sub_dir_sim] [sub_dir_model] [posterior_file] [beta or bimodal] [binpath]")
	print("\n\tex: submit_simulations_gof_2pop_popGrowth.py 1000 2 SI_1N flo mal sim_SI_1N SI_1N posterior_IM_1M_2N.txt beta\n\tto simulate 1000 multilocus simulations at the second iteration, in the folder sim_SI_1N") 
	sys.exit(0)

outgroup = int(sys.argv[1])
nmultilocus = int(sys.argv[2]) # 10000
iteration = int(sys.argv[3]) # 2
model = sys.argv[4]
nameA = sys.argv[5] # name of the species A
nameB = sys.argv[6] # name of the species B
sub_dir_sim = sys.argv[7] # name of the subdir where the simulations will be run
sub_dir_model = sys.argv[8] # name of the sub_sub_dir containing ABCstat.txt
posterior_file = sys.argv[9]
timeStamp = sys.argv[10]
modeBarrier = sys.argv[11] # beta / bimodal
modePrior = sys.argv[12] # joint / disjoint / randomBeta
binpath = sys.argv[13]

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
if "SI" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs"
if "AM" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs"
if "SC" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs"
if "IM" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs"
if "PAN" in model:
	mscommand = "-t tbs -r tbs tbs -eN 0 tbs -eN tbs tbs"

if mscommand == "":
	print("You specified a wrong model: SI_x, AM_x, AM_x, SC_x, PAN_x\n")
	sys.exit()

#tmp = "mkdir {0}/{1}; ".format(path, sub_dir_sim)
#tmp += "mkdir {0}/{1}/{2}_{3}; ".format(path, sub_dir_sim, sub_dir_model, iteration)
tmp = "cp {0}/bpfile {0}/{1}/{2}_{3}; ".format(path, sub_dir_sim, sub_dir_model, iteration)
tmp += "cd {0}/{1}/{2}_{3}; ".format(path, sub_dir_sim, sub_dir_model, iteration)

#tmp += "{8}/priorgen_gof_2pop_popGrowth.py {0} {1} {2} {6} {7} | {8}/msnsam tbs {3} {4} | {8}/mscalc_2pop_SFS.py {5}".format(model, nmultilocus, posterior_file, nmultilocus*nlocus, mscommand, outgroup, modeBarrier, modePrior, binpath)
tmp += "python {0}/priorgen_gof_2pop_popGrowth.py {1} {2} {3} {4} {5} | {0}/msnsam tbs {6} {7} | pypy {0}/mscalc_2pop_SFS.py {8}".format(binpath, model, nmultilocus, posterior_file, modeBarrier, modePrior, nmultilocus*nlocus, mscommand, outgroup)
#tmp2 = 'sbatch --nodes=1 --ntasks-per-node=1 --time=02:00:00 -J {0}_{1} --wrap="{2}"\n'.format(model, iteration, tmp)

print(tmp)
os.system(tmp)

