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

def count_spike_pos(recfile):
    df = pd.read_table(recfile)
    bg_recom = df.iloc[0, 1]

    # how many variants NOT with background recom, and the sum
    #sum1 = df[df["SwitchRate"] != bg_recom]["SwitchRate"].sum()
    #count1 = df[df["SwitchRate"] != bg_recom]["SwitchRate"].count()
    df2 = df[df["SwitchRate"] != bg_recom]
    
    return df2

spike = count_spike_pos(filename)

spike.to_csv(filename + ".spike", sep="\t")
    
    