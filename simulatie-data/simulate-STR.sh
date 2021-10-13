#!/bin/bash

# depends on the next software being in PATH:
# SURVIVOR         https://github.com/fritzsedlazeck/SURVIVOR
# STRsimulator.py  https://github.com/DamarisLa/STRsimulator
# samtools         https://github.com/samtools/samtools
# minimap2         https://github.com/lh3/minimap2
# longshot         https://github.com/pjedge/longshot


echo "STRsimulator"
REFERENCE_FA=$1
STRSIM_FA="strsim-genome.fasta"
INPUT_BED=$2
STRSIM_BED="strsim.bed"

[[ "${REFERENCE_FA}" == "" ]] && echo "'REFERENCE_FA' is missing" && exit
[[ "${INPUT_BED}" == "" ]] && echo "'INPUT_BED' is missing" && exit

# STROsimulator parameters
#   1: Reference fasta file
#   2: output fasta file (simulated genome/chromosome)
#   3: Input Bedfile (sorted!!!)
#   4: Output Bedfile
#   5: [0.00-1.00] relative amount of STR Regions that should be modified
#   6: ["d"|"h"] h = haploid or d = diploid output genome
#   7: [0.00-1.00] chance of bases in the STR regions SNV-mutated (substituted)
#   8: [0-100 or more] times more unlikely Indel than substitution
#   9: [0.00-1.00] percentage of homozygous regions
STR_MOD="0.50"
HAP_DI="d"
SUBSTITUTIONS="0"
INDELS="0"
HOMOZYGOUS="0.35"

python3 STRsimulator.py \
  "${REFERENCE_FA}" \
  "${STRSIM_FA}" \
  "${INPUT_BED}" \
  "${STRSIM_BED}" \
  ${STR_MOD} \
  "${HAP_DI}" \
  ${SUBSTITUTIONS} \
  ${INDELS} \
  ${HOMOZYGOUS}

# add SNVs
SURV_PARAM="survivor-parameters.txt"
STR_WITH_SNV_FAS="strsim-snv-genome.fasta"
SURVIVOR simSV  ${STRSIM_FA}  ${SURV_PARAM}  0.001  0  ${STR_WITH_SNV_FAS}

echo "simulate reads with SURVIVOR"
ONT_STATS="ont-survivor-error-dist.tsv"
COVERAGE=30
READS="reads-survivor.fasta"
SURVIVOR simreads ${STR_WITH_SNV_FAS}.fasta ${ONT_STATS} ${COVERAGE} ${READS}


echo "map with minimap2"
SORTED_BAM_OUT="alignment-minimap2-sorted.bam"
THREADS="4"

minimap2  -ax map-ont \
    "${REFERENCE_FA}" \
    "${READS}" \
    "--MD" -t "${THREADS}" \
    -Y | \
samtools sort -@ ${THREADS} - > "${SORTED_BAM_OUT}"
samtools index ${SORTED_BAM_OUT}


echo "phase reads with lognshot"
PHASED_VCF="phaser-alignment-longshot-vars.vcf"
PHASED_BAM="phaser-alignment-longshot.bam"

longshot --bam ${SORTED_BAM_OUT} --ref ${REFERENCE_FA} --out ${PHASED_VCF} --out_bam ${PHASED_BAM}
