#!/shared/mfs/data/software/miniconda/envs/pypy-2.7-5.10.0/bin/pypy
# #!/usr/bin/pypy
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
from math import ceil
import random

# check the arguments
if len(sys.argv) != 12:
	print("\n\tfasta2ABC_1pops.py produces: bpfile (for simulations) and summary statistics (for inferences)")
	print("\n\033[1;33m\tExample: ./fasta2ABC_1pops.py all_loci.fasta flo coding 30 0.1 10 100000 0.00000002 1\033[0m\n")
	print("\t\targ1 =\tname of the fasta file containing all of the sequences")
	print("\t\targ2 =\tname of the project's directory (in practice: a random timestamp")
	print("\t\targ3 =\tID of species A (example: flo)")
	print("\t\targ4 =\tID of the outgroup (example: num)")
	print("\t\targ5 =\t'coding' or 'noncoding', to study only synonymous polymorphisms (if coding) or all SNPs (if noncoding)")
	print("\t\targ6 =\tminimum length of a locus to be considered, i.e, number of total positions minus the positions containing a N")
	print("\t\targ7 =\tvalue in [0-1]. Corresponds to a threshold of %N above which a sequence is rejected")
	print("\t\targ8 =\tinteger, corresponding to the minimum number of retained sequences (after rejection).\n\t\t\tif not enough sequences are retained, the loci is excluded from the analysis")
	print("\t\targ9 =\tmutation rate by bp and by generation. example: 0.00000002")
	print("\t\targ10 =\tratio of the recombination rate over mutation. example: 1")
	print("\t\targ11 =\tbinpath. example: /tools/ABConline/bin")
	if(len(sys.argv)<12):
		sys.exit("\n\033[1;31m ERROR in fasta2ABC_1pops.py: 12 arguments are required: {0} missing\033[0m\n".format(12-len(sys.argv)))
	if(len(sys.argv)>12):
		sys.exit("\n\033[1;31m ERROR in fasta2ABC_1pops.py: 12 arguments are required: {0} too much\033[0m\n".format(len(sys.argv)-12))

fileName = sys.argv[1] # example: all_loci.fasta
timeStamp = sys.argv[2] # example: timeStamp used as directory for the project
nameA = sys.argv[3] # name of species A. example: flo
nameOut = sys.argv[4] # name of the outgroup. example: num
region = sys.argv[5] # if == coding: will only deal with synonymous codons; if == noncoding: will deal with all positions
Lmin = int(sys.argv[6]) # minimum length for a locus to be retained. example: 30
max_N_tolerated = float(sys.argv[7]) # if an allele has %N > threshold_N --> sequence is rejected
nMin = int(sys.argv[8]) # minimum number of individuals within a species. example: 10
mu = float(sys.argv[9]) # mutation rate by bp and by generation. example: 0.00000002
rho_over_theta = float(sys.argv[10]) # ratio of the recombination rate over mutation. example: 1
binpath = sys.argv[11] # path to the bin directory of DILS containing all executables 

test = os.path.isfile(fileName)
if test == False:
	sys.exit("\n\t\033[1;31m ERROR in fasta2ABC_1pops.py: alignement '{0}' is not found\033[0m\n".format(fileName))

if region not in ['coding', 'noncoding']:
	sys.exit("\n\t\033[1;31m ERROR in fasta2ABC_1pops.py: '{0}' is not an expected argument. 'coding' or 'noncoding' are expected here\033[0m\n".format(region))

if os.path.isdir('{0}'.format(timeStamp)) == True:
	commande = 'rm -rf {0}'.format(timeStamp)
	os.system(commande)
commande = 'mkdir {0}'.format(timeStamp)
os.system(commande)

def coloredSeq(seq):
	# print sequences with the standard color code
	seq = seq.replace("A", '\x1b[5;30;41m' + 'A' + '\x1b[0m')
	seq = seq.replace("T", '\x1b[5;30;44m' + 'T' + '\x1b[0m')
	seq = seq.replace("G", '\x1b[6;30;43m' + 'G' + '\x1b[0m')
	seq = seq.replace("C", '\x1b[5;30;42m' + 'C' + '\x1b[0m')
	return(seq)


def trunc2triplets(size):
	# trunc a value to its closest and smaller multiple of 3
	size = int(size)
	for i in range(2):
		if size%3 != 0:
			size -= 1
	return(size)


# nN = number of non-synonymous sites in the codon i: example for CGG -> nN = 2/3 + 3/3 + 0/3
# nS = number of synonymous sites in the codon i: example for CGG -> n> = 1/3 + 0/3 + 3/3
codonTable = {'AAA': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAC': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAG': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAT': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ACA': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACC': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACG': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACT': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'AGA': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGC': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AGG': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGT': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ATA': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATC': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATG': {'aa': 'M', 'nN': 3.0, 'nS': 0.0}, 'ATT': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'CAA': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAC': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAG': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAT': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CCA': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCC': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCG': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCT': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CGA': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGC': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CGG': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGT': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CTA': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTC': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'CTG': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTT': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'GAA': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAC': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAG': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAT': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GCA': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCC': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCG': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCT': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GGA': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGC': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGG': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGT': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GTA': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTC': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTG': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTT': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'TAC': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TAT': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TCA': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCC': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCG': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCT': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TGC': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TGG': {'aa': 'W', 'nN': 3.0, 'nS': 0.0}, 'TGT': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTA': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTC': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTG': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTT': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}}


#def fasta2dic(fastaFile):
#	fasta = open(fastaFile).readlines()
#	seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
#	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
#	res = {}
#	for i in range(len(seq)):
#		res[seqName[i]] = seq[i]
#	return (res)

#alignA = fasta2dic(seqA)


def fasta2list(fastaFile, nameA, nameOut, nMin, max_N_tolerated):
	L = {}
	res = {} # res[species][locus]['seq', 'id']
	res[nameA] = {}
	if nameOut != 'NA':
		res[nameOut] = {}
	fasta = open(fastaFile).readlines()
	seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	
	nsam = {} # number of individuals in both species
	
	if nameOut == 'NA':	
		list_species = [nameA]
	else:
		list_species = [nameA, nameOut]
	for i in range(len(seqName)): # loop over loci
		tmp = seqName[i].split('|') # split a string similar to : Hmel219002_6|flo|flo.CS12|allele1
		locus = tmp[0]
		species = tmp[1]
		if locus not in nsam:
			nsam[locus] = {}
			for j in list_species:
				nsam[locus][j] = 0
		
		if species in list_species:
			if locus not in res[species]:
				res[species][locus] = {}
				res[species][locus]['seq'] = []
				res[species][locus]['id'] = []

			# remove sequences with too many N
			propN = seq[i].count("N")/(1.0 * len(seq[i]))
			if propN <= max_N_tolerated:
				res[species][locus]['seq'].append(seq[i])
				res[species][locus]['id'].append(seqName[i])
				nsam[locus][species] += 1
		
	# remove loci not found in sufficient individuals in both species
	for i in nsam: # loop along loci
		if nameOut == 'NA':
			if nsam[i][nameA] < nMin:
				for j in res.keys(): # loop over species
					if i in res[j]:
						del res[j][i] # delete locus i
			else:
				L[i] = len(res[nameA][i]['seq'][0]) - 3 # remove the last 3 bases to excluse final stop codon
				L[i] = trunc2triplets(L[i]) # convert the remaining length into a multiple of 3
				# sub sample nMin sequences from species A
				sub_sample_A = random.sample(range(nsam[i][nameA]), nMin)
				sub_seq_A = [ res[nameA][i]['seq'][seq] for seq in sub_sample_A ]
				sub_id_A = [ res[nameA][i]['id'][seq] for seq in sub_sample_A ]
				res[nameA][i]['seq'] = sub_seq_A
				res[nameA][i]['id'] = sub_id_A
				nsam[i][nameA] = nMin
		else:
			if nsam[i][nameA] < nMin or nsam[i][nameOut]==0:
				for j in res.keys(): # loop over species
					if i in res[j]:
						del res[j][i] # delete locus i
			else:
				L[i] = len(res[nameA][i]['seq'][0]) - 3 # remove the last 3 bases to excluse final stop codon
				L[i] = trunc2triplets(L[i]) # convert the remaining length into a multiple of 3
				# sub sample nMin sequences from species A
				sub_sample_A = random.sample(range(nsam[i][nameA]), nMin)
				sub_seq_A = [ res[nameA][i]['seq'][seq] for seq in sub_sample_A ]
				sub_id_A = [ res[nameA][i]['id'][seq] for seq in sub_sample_A ]
				res[nameA][i]['seq'] = sub_seq_A
				res[nameA][i]['id'] = sub_id_A
				nsam[i][nameA] = nMin
	return ({'align': res, 'L': L})


def getConsensus(align, L):
	# align[locus]['seq', 'id']
	# L[locus] = n nucleotides
	consensus = {}
	for locus in align: # loop over loci
		L_locus = L[locus]
		consensus[locus]='' # consensus[locus] = sequence
		for pos in range(L_locus):
			position = []
			for sequence in align[locus]['seq']:
				base = sequence[pos]
				if base != 'N':
					position.append(base)
			if len(position) == 0:
				consensus[locus] += 'N'
			else:
				consensus[locus] += random.sample(position, 1)[0]
	return(consensus)


def getScalar(align, L, consensus, region, nameA):
	# align[species][locus]['seq', 'id']
	# L[locus]
	# consensus[locus]
	scalar = {} # scalar[locus] = divergence_of_the_locus / mean(divergence)
	divergence = [] # array containing all of the divergence_i of all loci. Used to compute de the mean(divergence)
	nLocus = 0
	for locus in L: # loop over loci
		if locus in align[nameA]:
			nLocus += 1
			divergence_locus = 0.0
			if region == 'coding':
				# use k3, the divergence at third coding position
				list_positions = range(2, L[locus], 3)
			else:
				# use all positions
				list_positions = range(L[locus])
			nComb = 0
			for pos in list_positions: # loop over third coding positions
				out = consensus[locus][pos]
				for seqA in align[nameA][locus]['seq']: # loop over indA
					A = seqA[pos]
					if A!='N' and out!='N':
						nComb += 1
						if A!=out:
							divergence_locus += 1.0
			# end of loop over locus
			if nComb > 0:
				divergence_locus /= nComb # hypothesis : if no possible comparison between the ingroup and the outgroup, then divergence = 0
			else:
				divergence_locus = 0.0
			divergence.append(divergence_locus)
			scalar[locus] = divergence_locus
			
	# compute the average divergence
	mean = 0.0
	if nLocus > 0:
		for i in range(nLocus):
			mean += divergence[i]
		mean /= nLocus
	else:
		mean = 1
	
	# attribute a scalar of mutation rate for each locus
	for locus in scalar:
			scalar[locus] /= mean 
			# to constrain the range of mutation rates in [0.1 - 10] x mean(divergence)
			if scalar[locus]<0.1:
				scalar[locus] = 0.1
			if scalar[locus]>10:
				scalar[locus] = 10

	return(scalar)


# read the input file
align = fasta2list(fileName, nameA, nameOut, nMin, max_N_tolerated)  # align[species][locus]['id', 'seq']


if len(align['align'][nameA]) == 0:
	sys.exit('\n\tERROR in fasta2ABC_1pops.py: no locus found in file {0} corresponding to a correct alignement within {1}\n'.format(fileName, nameA))
	if nameOut != 'NA':
		if len(align['align'][nameOut]) == 0:
			sys.exit('\n\tERROR in fasta1ABC_2pops.py: no locus found in file {0} corresponding {1}\n'.format(fileName, nameOut))


# if there is an outgroup -> make a concensus
if nameOut != 'NA':
	consensus = getConsensus(align['align'][nameOut], align['L']) # consensus[locus] = sequence
	scalar = getScalar(align['align'], align['L'], consensus, region, nameA) # with outgroup: locus specific mutation rate mu_i = mean(mu) * div_i / mean(div)
else:
	scalar = {}
	for i in align['align'][nameA]:
		scalar[i] = 1 # without outgroup: all loci share the same per nucleotide mutation rate


# treat the input file
bpfile_L1 = '# spA={0} mu={1}'.format(nameA, mu)
bpfile_L2 = [] # L
bpfile_L3 = [] # nA
bpfile_L4 = [] # theta
bpfile_L5 = [] # rho

output_ms = "./msnsam tbs 20 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 1 2 -eN tbs tbs\n3579 27011 59243\n\n"
outfile_ms = open('{0}/{1}.ms'.format(timeStamp, nameA), 'w')
outfile_ms.write(output_ms)

if nameOut == 'NA':
	# For coding loci
	if region == 'coding':
		output_info = "locusName\tL_including_N\tLsyno\tnSynSegSite\tnsamA\n"
		outfile_info = open('{0}/{1}_infos.txt'.format(timeStamp, nameA), 'w')
		outfile_info.write(output_info)
		for locus_i in align['L'].keys():
			geneName = locus_i
			
			nA = len(align['align'][nameA][locus_i]['id'])
			
			L = align['L'][locus_i] 
			interspe = [] # contains the interspecific alignment
			interspeName = [] # contains the id of sequences in the interspecific alignment
			for i in range(nA):
				interspe.append(align['align'][nameA][locus_i]['seq'][i])
				interspeName.append(align['align'][nameA][locus_i]['id'][i])

			nSites = 0 # total number of synonymous sites within the sequence, computed using codonTable
			nSynSegSite = 0 # number of synonymous segregating sites among the nSites
			positions = [] # list of synonymous polymorphic positions: doesn't correspond to the SNP position, but to the first codon position
			msStyle = [] # contains the msStyle format
			for ind in range(nA):
				msStyle.append([])

			# loop over codons:
			for pos in range(L)[::3]:
				alignmentOfCodons = [] # set of codons in the alignment, starting at the position 'pos1'
				# loop over individuals:
				# get all codons in the alignment
				for ind in range(nA):
					pos1 = interspe[ind][pos]
					pos2 = interspe[ind][pos + 1]
					pos3 = interspe[ind][pos + 2]
					base = pos1 + pos2 + pos3 
					alignmentOfCodons.append(base)
				
				polyMcodons = list(set(alignmentOfCodons)) # list of codons found in the alignment
				nCodons = 0
				nCodons = len(polyMcodons)
				testN = False # False if no codon with 'N'; True if a 'N' is found in at least one codon for one individual
				testStopCodon = False # False if no stop codon was found; True if a stop codon was found
				for i in polyMcodons: # loop to test for some 'N'
					if 'N' in i:
						testN = True
					if i not in codonTable:
						testStopCodon = True
				
				# if: 1) a maximum of 2 polymorphic codons, and, 2) no codon with 'N', and, 3) all codons effectively code for an amino acid
				if nCodons <= 2 and testN==False and testStopCodon==False: 
					nSites_pos = 0.0
					for i in alignmentOfCodons:
						nSites_pos += codonTable[i]['nS']
					nSites += nSites_pos/len(alignmentOfCodons)
					
					# if two codons --> there is a polymorphism
					if nCodons == 2:
						alignmentOfAminoAcids = []
						for i in alignmentOfCodons:
							alignmentOfAminoAcids.append(codonTable[i]['aa'])
						setOfAminoAcids = list(set(alignmentOfAminoAcids))
						
						# if two codons but one amino acids --> synonymous polymorphism
						if len(setOfAminoAcids) == 1:
							nSynSegSite += 1
							positions.append(pos) # positions: list of first codon position of polymorphic synonymous codons
							ancestralAllele = polyMcodons[0] # in absence of outgroup --> the ancestral allele is the first in the alignement
							derivedAllele = polyMcodons[1] # without outgroup --> the derived allele is the one who is not the first...
							for i in range(nA):
								if alignmentOfCodons[i] == ancestralAllele:
									msStyle[i].append('0')
								if alignmentOfCodons[i] == derivedAllele:
									msStyle[i].append('1')

			if nSites >= Lmin: # if the locus is big enough to be considered
				# ms_like output files
				locus_ms = ''
				locus_ms = locus_ms + "//{0}\n".format(geneName)
				locus_ms = locus_ms + "segsites: {0}\n".format(int(nSynSegSite))
				if nSynSegSite != 0:
					locus_ms += "positions: {0}\n".format( " ".join([ str(round((1.0*i)/L, 4)) for i in positions ]))
					for i in msStyle:
						locus_ms = locus_ms + "".join( [ str(j) for j in i ] ) + "\n"
				
				outfile_ms.write(locus_ms + '\n')
				
				bpfile_L2.append(int(ceil(nSites)))
				bpfile_L3.append(nA)
				bpfile_L4.append(mu) #
				bpfile_L5.append(mu)
					
				# informations about locus
				res = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(geneName, L, nSites, nSynSegSite, nA)
				outfile_info.write(res)
				
			#	res = ""
			#	for i in range(len(interspe)):
			#		res += ">{0}\n{1}\n".format(interspeName[i], interspe[i])
			#	outfile = open('{0}.fas'.format(geneName), "w")
			#	outfile.write(res)
			#	outfile.close()
		outfile_ms.close()
		outfile_info.close()

		bpfile = bpfile_L1 + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L2 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L3 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L4 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L5 ]) + '\n'

		outfile = open('{0}/nLoci.txt'.format(timeStamp), 'w')
		outfile.write('{0}'.format(len(bpfile_L2)))
		outfile.close()

		outfile = open('{0}/bpfile'.format(timeStamp), 'w')
		outfile.write(bpfile)
		outfile.close()
		
	# For non coding loci
	if region == 'noncoding':
		output_info = "locusName\tL_including_N\tL\tnSegSite\tnsamA\n"
		outfile_info = open('{0}/{1}_{2}_infos.txt'.format(timeStamp, nameA), 'w')
		outfile_info.write(output_info)
		for locus_i in align['L'].keys():
			geneName = locus_i
			
			nA = len(align['align'][nameA][locus_i]['id'])
			
			L = align['L'][locus_i] 
			interspe = [] # contains the interspecific alignment
			interspeName = [] # contains the id of sequences in the interspecific alignment
			for i in range(nA):
				interspe.append(align['align'][nameA][locus_i]['seq'][i])
				interspeName.append(align['align'][nameA][locus_i]['id'][i])

			nSites = 0 # total number of sites within the sequence
			nSegSite = 0 # number of segregating sites among the nSites
			positions = [] # list of polymorphic positions: correspond to the SNP position
			msStyle = [] # contains the msStyle format
			for ind in range(nA):
				msStyle.append([])

			# loop over pos:
			for pos in range(L):
				alignmentOfPos = [] # set of pos in the alignment, starting at the position 'pos1'
				# loop over individuals:
				# get all pos in the alignment
				for ind in range(nA):
					pos1 = interspe[ind][pos]
					base = pos1 
					alignmentOfPos.append(base)
				
				polyMpos = list(set(alignmentOfPos)) # list of pos found in the alignment
				nPos = 0
				nPos = len(polyMpos)
				testN = False # False if no codon with 'N'; True if a 'N' is found in at least one codon for one individual
				for i in polyMpos: # loop to test for some 'N'
					if 'N' in i:
						testN = True
				
				# if: 1) a maximum of 2 polymorphic pos, and, 2) no codon with 'N'
				if nPos <= 2 and testN==False: 
					nSites += 1
					
					# if two pos --> there is a polymorphism
					if nPos == 2:
						nSegSite += 1
						positions.append(pos) # positions: list of first codon position of polymorphic synonymous pos
						ancestralAllele = polyMpos[0] # in absence of outgroup --> the ancestral allele is the first in the alignement
						derivedAllele = polyMpos[1] # without outgroup --> the derived allele is the one who is not the first...
						for i in range(nA):
							if alignmentOfPos[i] == ancestralAllele:
								msStyle[i].append('0')
							if alignmentOfPos[i] == derivedAllele:
								msStyle[i].append('1')

			if nSites >= Lmin: # if the locus is big enough to be considered
				# ms_like output files
				locus_ms = ''
				locus_ms = locus_ms + "//{0}\n".format(geneName)
				locus_ms = locus_ms + "segsites: {0}\n".format(int(nSegSite))
				if nSegSite != 0:
					locus_ms += "positions: {0}\n".format( " ".join([ str(round((1.0*i)/L, 4)) for i in positions ]))
					for i in msStyle:
						locus_ms = locus_ms + "".join( [ str(j) for j in i ] ) + "\n"
				
				outfile_ms.write(locus_ms + '\n')
				
				bpfile_L2.append(int(ceil(nSites)))
				bpfile_L3.append(nA)
				bpfile_L4.append(mu) # 
				bpfile_L5.append(mu*rho_over_theta)
					
				# informations about locus
				res = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(geneName, L, nSites, nSegSite, nA)
				outfile_info.write(res)
				
			#	res = ""
			#	for i in range(len(interspe)):
			#		res += ">{0}\n{1}\n".format(interspeName[i], interspe[i])
			#	outfile = open('{0}.fas'.format(geneName), "w")
			#	outfile.write(res)
			#	outfile.close()
		outfile_ms.close()
		outfile_info.close()

		bpfile = bpfile_L1 + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L2 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L3 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L4 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L5 ]) + '\n'

		outfile = open('{0}/nLoci.txt'.format(timeStamp), 'w')
		outfile.write('{0}'.format(len(bpfile_L2)))
		outfile.close()
		
		outfile = open('{0}/bpfile'.format(timeStamp), 'w')
		outfile.write(bpfile)
		outfile.close()
else: # if there is an outgroup
	# For coding loci
	if region == 'coding':
		output_info = "locusName\tL_including_N\tLsyno\tnSynSegSite\tnsamA\tmutation_scalar\n"
		outfile_info = open('{0}/{1}_infos.txt'.format(timeStamp, nameA), 'w')
		outfile_info.write(output_info)
		for locus_i in align['L'].keys():
			if locus_i in consensus:
				geneName = locus_i
				
				nA = len(align['align'][nameA][locus_i]['id'])
				L = align['L'][locus_i] 
				interspe = [] # contains the interspecific alignment
				interspeName = [] # contains the id of sequences in the interspecific alignment
				for i in range(nA):
					interspe.append(align['align'][nameA][locus_i]['seq'][i])
					interspeName.append(align['align'][nameA][locus_i]['id'][i])

				nSites = 0 # total number of synonymous sites within the sequence, computed using codonTable
				nSynSegSite = 0 # number of synonymous segregating sites among the nSites
				positions = [] # list of synonymous polymorphic positions: doesn't correspond to the SNP position, but to the first codon position
				msStyle = [] # contains the msStyle format
				for ind in range(nA):
					msStyle.append([])

				# loop over codons:
				for pos in range(L)[::3]:
					alignmentOfCodons = [] # set of codons in the alignment, starting at the position 'pos1'
					codon_outgroup = consensus[locus_i][pos:(pos+3)]
					# loop over individuals:
					# get all codons in the alignment
					for ind in range(nA):
						pos1 = interspe[ind][pos]
						pos2 = interspe[ind][pos + 1]
						pos3 = interspe[ind][pos + 2]
						base = pos1 + pos2 + pos3 
						alignmentOfCodons.append(base)
					
					polyMcodons = list(set(alignmentOfCodons)) # list of codons found in the alignment
					
					nCodons = 0
					nCodons = len(polyMcodons)
					testN = False # False if no codon with 'N'; True if a 'N' is found in at least one codon for one individual
					testStopCodon = False # False if no stop codon was found; True if a stop codon was found
					for i in polyMcodons: # loop to test for some 'N'
						if 'N' in i:
							testN = True
						if i not in codonTable:
							testStopCodon = True
					
					# if: 1) a maximum of 2 polymorphic codons, and, 2) no codon with 'N', and, 3) all codons effectively code for an amino acid
					if nCodons <= 2 and codon_outgroup in polyMcodons and testN==False and testStopCodon==False: 
						nSites_pos = 0.0
						for i in alignmentOfCodons:
							nSites_pos += codonTable[i]['nS']
						nSites += nSites_pos/len(alignmentOfCodons)
						
						# if two codons --> there is a polymorphism
						if nCodons == 2:
							alignmentOfAminoAcids = []
							for i in alignmentOfCodons:
								alignmentOfAminoAcids.append(codonTable[i]['aa'])
							setOfAminoAcids = list(set(alignmentOfAminoAcids))
							
							# if two codons but one amino acids --> synonymous polymorphism
							if len(setOfAminoAcids) == 1:
								nSynSegSite += 1
								positions.append(pos) # positions: list of first codon position of polymorphic synonymous codons
								#ancestralAllele = polyMcodons[0] # in absence of outgroup --> the ancestral allele is the first in the alignement
								#derivedAllele = polyMcodons[1] # without outgroup --> the derived allele is the one who is not the first...
								ancestralAllele = codon_outgroup # in absence of outgroup --> the ancestral allele is the first in the alignement
								derivedAllele = polyMcodons[abs(1-polyMcodons.index(codon_outgroup))] # without outgroup --> the derived allele is the one who is not the first...
								for i in range(nA):
									if alignmentOfCodons[i] == ancestralAllele:
										msStyle[i].append('0')
									if alignmentOfCodons[i] == derivedAllele:
										msStyle[i].append('1')

				if nSites >= Lmin: # if the locus is big enough to be considered
					# ms_like output files
					locus_ms = ''
					locus_ms = locus_ms + "//{0}\n".format(geneName)
					locus_ms = locus_ms + "segsites: {0}\n".format(int(nSynSegSite))
					if nSynSegSite != 0:
						locus_ms += "positions: {0}\n".format( " ".join([ str(round((1.0*i)/L, 4)) for i in positions ]))
						for i in msStyle:
							locus_ms = locus_ms + "".join( [ str(j) for j in i ] ) + "\n"
					
					outfile_ms.write(locus_ms + '\n')
					bpfile_L2.append(int(ceil(nSites)))
					bpfile_L3.append(nA)
					bpfile_L4.append(mu*scalar[locus_i]) #
					bpfile_L5.append(mu*rho_over_theta*scalar[locus_i])
						
					# informations about locus
					res = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(geneName, L, nSites, nSynSegSite, nA, scalar[locus_i])
					outfile_info.write(res)
					
				#	res = ""
				#	for i in range(len(interspe)):
				#		res += ">{0}\n{1}\n".format(interspeName[i], interspe[i])
				#	outfile = open('{0}.fas'.format(geneName), "w")
				#	outfile.write(res)
				#	outfile.close()
		outfile_ms.close()
		outfile_info.close()

		bpfile = bpfile_L1 + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L2 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L3 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L4 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L5 ]) + '\n'

		outfile = open('{0}/nLoci.txt'.format(timeStamp), 'w')
		outfile.write('{0}'.format(len(bpfile_L2)))
		outfile.close()
		
		outfile = open('{0}/bpfile'.format(timeStamp), 'w')
		outfile.write(bpfile)
		outfile.close()
	
	# For coding loci
	if region == 'noncoding':
		output_info = "locusName\tL_including_N\tL\tnSegSite\tnsamA\tmutation_scalar\n"
		outfile_info = open('{0}/{1}_infos.txt'.format(timeStamp, nameA), 'w')
		outfile_info.write(output_info)
		for locus_i in align['L'].keys():
			if locus_i in consensus:
				geneName = locus_i
				
				nA = len(align['align'][nameA][locus_i]['id'])
				
				L = align['L'][locus_i] 
				interspe = [] # contains the interspecific alignment
				interspeName = [] # contains the id of sequences in the interspecific alignment
				for i in range(nA):
					interspe.append(align['align'][nameA][locus_i]['seq'][i])
					interspeName.append(align['align'][nameA][locus_i]['id'][i])

				nSites = 0 # total number of synonymous sites within the sequence, computed using codonTable
				nSynSegSite = 0 # number of synonymous segregating sites among the nSites
				positions = [] # list of synonymous polymorphic positions: doesn't correspond to the SNP position, but to the first codon position
				msStyle = [] # contains the msStyle format
				for ind in range(nA):
					msStyle.append([])

				# loop over codons:
				for pos in range(L):
					alignmentOfCodons = [] # set of codons in the alignment, starting at the position 'pos1'
					codon_outgroup = consensus[locus_i][pos:(pos+3)]
					# loop over individuals:
					# get all codons in the alignment
					for ind in range(nA):
						pos1 = interspe[ind][pos]
						pos2 = interspe[ind][pos + 1]
						pos3 = interspe[ind][pos + 2]
						base = pos1 + pos2 + pos3 
						alignmentOfCodons.append(base)
					
					polyMcodons = list(set(alignmentOfCodons)) # list of codons found in the alignment
					
					nCodons = 0
					nCodons = len(polyMcodons)
					testN = False # False if no codon with 'N'; True if a 'N' is found in at least one codon for one individual
					testStopCodon = False # False if no stop codon was found; True if a stop codon was found
					for i in polyMcodons: # loop to test for some 'N'
						if 'N' in i:
							testN = True
						if i not in codonTable:
							testStopCodon = True
					
					# if: 1) a maximum of 2 polymorphic codons, and, 2) no codon with 'N', and, 3) all codons effectively code for an amino acid
					if nCodons <= 2 and codon_outgroup in polyMcodons and testN==False and testStopCodon==False: 
						nSites_pos = 0.0
						for i in alignmentOfCodons:
							nSites_pos += codonTable[i]['nS']
						nSites += nSites_pos/len(alignmentOfCodons)
						
						# if two codons --> there is a polymorphism
						if nCodons == 2:
							alignmentOfAminoAcids = []
							for i in alignmentOfCodons:
								alignmentOfAminoAcids.append(codonTable[i]['aa'])
							setOfAminoAcids = list(set(alignmentOfAminoAcids))
							
							# if two codons but one amino acids --> synonymous polymorphism
							if len(setOfAminoAcids) == 1:
								nSynSegSite += 1
								positions.append(pos) # positions: list of first codon position of polymorphic synonymous codons
								#ancestralAllele = polyMcodons[0] # in absence of outgroup --> the ancestral allele is the first in the alignement
								#derivedAllele = polyMcodons[1] # without outgroup --> the derived allele is the one who is not the first...
								ancestralAllele = codon_outgroup # in absence of outgroup --> the ancestral allele is the first in the alignement
								derivedAllele = polyMcodons[abs(1-polyMcodons.index(codon_outgroup))] # without outgroup --> the derived allele is the one who is not the first...
								for i in range(nA):
									if alignmentOfCodons[i] == ancestralAllele:
										msStyle[i].append('0')
									if alignmentOfCodons[i] == derivedAllele:
										msStyle[i].append('1')

				if nSites >= Lmin: # if the locus is big enough to be considered
					# ms_like output files
					locus_ms = ''
					locus_ms = locus_ms + "//{0}\n".format(geneName)
					locus_ms = locus_ms + "segsites: {0}\n".format(int(nSynSegSite))
					if nSynSegSite != 0:
						locus_ms += "positions: {0}\n".format( " ".join([ str(round((1.0*i)/L, 4)) for i in positions ]))
						for i in msStyle:
							locus_ms = locus_ms + "".join( [ str(j) for j in i ] ) + "\n"
					
					outfile_ms.write(locus_ms + '\n')
					
					bpfile_L2.append(int(ceil(nSites)))
					bpfile_L3.append(nA)
					bpfile_L4.append(mu*scalar[locus_i]) #
					bpfile_L5.append(mu*rho_over_theta*scalar[locus_i])
						
					# informations about locus
					res = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(geneName, L, nSites, nSegSite, nA, scalar[locus_i])
					outfile_info.write(res)
					
				#	res = ""
				#	for i in range(len(interspe)):
				#		res += ">{0}\n{1}\n".format(interspeName[i], interspe[i])
				#	outfile = open('{0}.fas'.format(geneName), "w")
				#	outfile.write(res)
				#	outfile.close()
		outfile_ms.close()
		outfile_info.close()

		bpfile = bpfile_L1 + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L2 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L3 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L4 ]) + '\n'
		bpfile += '\t'.join([ str(i) for i in bpfile_L5 ]) + '\n'

		outfile = open('{0}/nLoci.txt'.format(timeStamp), 'w')
		outfile.write('{0}'.format(len(bpfile_L2)))
		outfile.close()
		
		outfile = open('{0}/bpfile'.format(timeStamp), 'w')
		outfile.write(bpfile)
		outfile.close()


# compute the summary statistics for ABC
if nameOut == 'NA':
	outgroup_present = 0
else:
	outgroup_present = 1
#commande = 'cat {0}/{1}.ms | mscalc_1pop_observedDataset.py {0} {2}'.format(timeStamp, nameA, outgroup_present)
commande = 'cat {0}/{1}.ms | pypy {3}/mscalc_1pop_observedDataset_SFS.py {0} {2}'.format(timeStamp, nameA, outgroup_present, binpath)
#print(commande)
os.system(commande)

# remove the useless ms file
commande = 'rm {0}/{1}.ms'.format(timeStamp, nameA)
#print(commande)
#os.system(commande)

