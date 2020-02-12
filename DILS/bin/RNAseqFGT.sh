#!/usr/bin/bash
params_nameA=$1
params_nameB=$2
input_infos=$3
input_fasta=$4
input_nLoci=$5
output_tmp=$6
output_results=$7
binpath=$8


# RNAseqFGT.sh {params.nameA} {params.nameB} {input.infos} {input.fasta} {input.nLoci} {output.tmp} {output.results} {binpath}

x=$(cat ${input_infos}  | grep -v locusName | awk '{print $2}' | sort -n | tail -n1)
if (($x<99990)); then
	${binpath}/RNAseqFGT ${input_fasta} ${output_tmp}
	head -n1 ${output_tmp} > ${output_results}
	cat ${output_tmp} | grep ${params_nameA} >> ${output_results}
	cat ${output_tmp} | grep ${params_nameB} >> ${output_results}
	
else

	echo "${input_fasta} not in correct FASTA format (line length > 100000)" > ${output_results}
	echo "${input_fasta} not in correct FASTA format (line length > 100000)" > ${output_tmp}
fi

return 0

