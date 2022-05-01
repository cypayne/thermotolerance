#!/usr/bin/env python

'''
usage: ./replace-chrom-name.py chrom-name-number.txt infile.tsv
'''

import sys
import csv
import pandas as pd

name2num = pd.read_csv(sys.argv[1], delimiter='\t')
#print(name2num.columns.tolist())
inf = pd.read_csv(sys.argv[2], delimiter='\t')
name2num_dict = dict(zip(list(name2num['chr']), list(name2num['group'])))
inf['chr'] = inf['chr'].map(name2num_dict)

inf.to_csv(sys.argv[2]+'_chr-renamed.csv',index=False)
