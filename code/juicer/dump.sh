#!/bin/bash

mkdir -p fastq

cat srr/SRR_Acc_List.txt | while read SRR; do
  echo "Processing $SRR"
  fasterq-dump "srr/$SRR" \
    --outdir fastq \
    --split-files \
    --threads 8
done
