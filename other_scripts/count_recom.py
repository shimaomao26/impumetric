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

def count_spike_recom(recfile):
    df = pd.read_table(recfile)
    bg_recom = df.iloc[0, 1]
    
    # how many variants with background recom, and the sum
    sum0 = df[df["SwitchRate"] == bg_recom]["SwitchRate"].sum()
    count0 = df[df["SwitchRate"] == bg_recom]["SwitchRate"].count()
    
    # how many variants NOT with background recom, and the sum
    sum1 = df[df["SwitchRate"] != bg_recom]["SwitchRate"].sum()
    count1 = df[df["SwitchRate"] != bg_recom]["SwitchRate"].count()
    
    #total recom
    sum2 = sum0 + sum1
    
    return bg_recom, count0, sum0, count1, sum1, sum2

recom, a, b, c, d, e = count_spike_recom(filename)

print("\t".join([str(i) for i in [recom, a, b, c, d, e] ])  )
