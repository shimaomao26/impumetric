#!/bin/bash

export PATH=/home/shimm/software:$PATH

echo -e "JOB_ID: ${JOB_ID}"

set -eu

inpath=$1
outpath=$2
variants=overlapping_155297.idv
sample_file=jewel.n993.iid

# print
bcftools view \
  --include ID==@$variants \
  --force-samples \
  -S $sample_file \
  -Ou $inpath | bcftools query \
  -f'%CHROM\t%POS\t%REF\t%ALT[\t%HDS]\n' > $outpath
  