#!/home/shi/anaconda3/bin/python
# -*- coding: utf-8 -*-


#$ -S /home/shimm/anaconda3/bin/python
#$ -q '!mjobs_rerun.q'
#$ -cwd
#$ -pe def_slot 1
#$ -l s_vmem=20G
#$ -V


import numpy as np
import pandas as pd

from sys import argv

script, filename = argv

def count_error(errfile):
    df = pd.read_table(errfile)
    
    sum2 = df["ErrorRate"].sum()/len(df)
    
    return sum2

sum2 = count_error(filename)

print(sum2)
