#!/bin/bash

trial=$1

/home/shimm/software/minimac4-1.0.2-Linux/bin/minimac4 \
  --refHaps ./1000Gp3v5_JEWEL_3k_chr19_noLCR_noUnMatchInsDel.${trial}.m3vcf \
  --haps ./eagle.chr19.993.vcf.gz \
  --prefix ./chr19_test.${trial} \
  --format GT,DS,HDS,GP,SD \
  --meta \
  --log \
  --cpus 10 \
  --params
