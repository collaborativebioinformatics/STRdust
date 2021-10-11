#!/bin/bash

# -a:   out is SAM
# -x:   preset of data source (ont, pacbio)
# -t:   number of threads to use (def = 3)
# -Y:   use soft clipping for supp. alignments
# --MD: add MD tag

REFERENCE=$1
READS=$2
SORTED_BAM_OUT=$3
THREADS=$4

minimap2  -ax map-ont \
    "${REFERENCE}" \
    "${READS}" \
    "--MD" -t "${THREADS}" \
    -Y | \
samtools sort -@ ${THREADS} - > "${SORTED_BAM_OUT}"
