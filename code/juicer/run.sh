#!/bin/bash

rm -rf aligned/
# rm -rf splits/

./scripts/juicer.sh \
  -D "$(pwd)" \
  -z "$(pwd)/references/dm3.fa" \
  -p "$(pwd)/references/dm3.chrom.sizes" \
  -y "$(pwd)/restriction_sites/dm3_DpnII.txt" \
  -s DpnII \
  -g dm3 \
  -t 12
  # -S chimeric
