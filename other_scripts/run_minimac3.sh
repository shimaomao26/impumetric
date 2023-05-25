#!/bin/bash

trial=$1
# sample id file, like EAS.EUR
run=$2

/home/shimm/software/Minimac3Executable/bin/Minimac3-omp \
  --refHaps ./jewel3k.chr19_test.${trial}.vcf.gz \
  --processReference \
  --prefix ./jewel3k.chr19_test.${trial}.run${run} \
  --log

