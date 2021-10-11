#!/bin/bash

# run STRsimulation
REFERENCE_FA=$1
OUTPUT_FA=$2
INPUT_BED=$3
OUTPUT_BED=$4
STR_MOD=$5
HAP_DI=$6
SUBSTITUTIONS=$7
INDELS=$8
HOMOZYGOUS=$9

python ./STRsimulator_v7.py 	"${REFERENCE_FA}" "${OUTPUT_FA}" "${INPUT_BED}" "${OUTPUT_BED}" \
			${STR_MOD} "${HAP_DI}" ${SUBSTITUTIONS}  ${INDELS} ${HOMOZYGOUS} 

# No parameters provided:
#   1: Reference fasta file
#   2: output fasta file (simulated genome/chromosome)
#   3: Input Bedfile (sorted!!!) 
#   4: Output Bedfile 
#   5: [0.00-1.00] relative amount of STR Regions that should be modified 
#   6: ["d"|"h"] h = haploid or d = diploid output genome
#   7: [0.00-1.00] chance of bases in the STR regions SNV-mutated (substituted)  
#   8: [0-100 or more] times more unlikely Indel than substitution 
#   9: [0.00-1.00] percentage of homozygous regions
