#!/bin/bash

export PATH=/home/shimm/software:$PATH

echo -e "JOB_ID: ${JOB_ID}"

set -eu

inpath=$1
outpath=$2
variants=overlapping_155297.idv
sample_file=jewel.n993.iid

# id
bcftools view \
  -r $variants \
  --force-samples \
  -S $sample_file \
  -Ou $inpath | bcftools view \
  -h | tail -n 1 | tr "\t" "\n" | tail -n+10 > $outpath.ids