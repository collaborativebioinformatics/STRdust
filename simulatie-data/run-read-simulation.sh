#!/bin/bash

# read simulation
REFERENCE=$1
SIM_STATS=$2
OUT_READS=$3
COVERAGE=$4

survivor simreads ${REFERENCE} ${SIM_STATS} ${COVERAGE} ${OUT_READS}

# No parameters provided:
#   1: Reference fasta file
#   2: error profile file (see scan)
#   3: Coverage
#   4: output prefix
