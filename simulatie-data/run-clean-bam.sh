#!/bin/bash

# Clean based onf quality and flags
# q: align quality >= 30
# F: ignore reads with flag bits (1796)
#        read unmapped (4)
#        not primary align (256)
#        read fails vendor qual check (512)
#        read is PCR dup (1024)
INBAM=$1
OUTBAM=$2
echo "Only reads pass filter"
samtools view -q 30 -F 1796 -b -h ${INBAM} > ${OUTBAM}
