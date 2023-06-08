#!/home/shi/anaconda3/bin/python
# -*- coding: utf-8 -*-

#$ -S /home/shimm/anaconda3/bin/python
#$ -q '!mjobs_rerun.q'
#$ -cwd
#$ -pe def_slot 1
#$ -l s_vmem=1G
#$ -V

from sys import argv

import random

script, in_filename, out_filename, downsample = argv

with open(in_filename,'r') as f:
    samples = []
    for line in f:
        samples.append(line.rstrip("\n"))
    
    
    random.shuffle(samples)
    
    
    textfile = open(out_filename, "w")
    
    for element in samples:
        textfile.write(str(element) + "\n")
    textfile.close()
    
    if int(downsample) < len(samples):
        textfile = open(out_filename + ".downsample.n" + downsample, "w")
        for element in samples[:int(downsample)]:
            textfile.write(str(element) + "\n")
        textfile.close()