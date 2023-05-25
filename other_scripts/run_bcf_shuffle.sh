#!/bin/bash

trial=$1
run=$2
# sample id file, like EAS.EUR

bcftools \
  view  --samples-file ./sample_shuffle/sample.shuffle1.${trial}.txt.${run} \
  --threads 1 \
  -Oz -o ./vcf_shuffle/jewel3k.chr19_test.shuffle1.${trial}.run${run}.vcf.gz 1000Gp3v5_JEWEL_3k_chr19_noLCR_noUnMatchInsDel.vcf.gz

