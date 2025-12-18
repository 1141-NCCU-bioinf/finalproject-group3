#!/bin/bash

echo "merging S2R_plus_R1..."
cat SRR5579177_1.fastq SRR5579178_1.fastq | pigz -c -p 8 >S2R_plus_R1.fastq.gz
echo "S2R_plus_R1 done!"

echo "merging S2R_plus_R1..."
cat SRR5579177_2.fastq SRR5579178_2.fastq | pigz -c -p 8 >S2R_plus_R2.fastq.gz
echo "S2R_plus_R2 done!"
