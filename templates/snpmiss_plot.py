#!/usr/bin/env python3
#Load SNP frequency file and generate histogram

import pandas as   pd
import numpy  as np
import sys
from matplotlib import use
use('Agg')
import argparse
import matplotlib
import matplotlib.pyplot as plt
import math

def parseArguments():
    if len(sys.argv)<=1:
        sys.argv="snpmissPlot.py $input $label $snpmiss_plot $cut_geno".split()
    parser=argparse.ArgumentParser()
    parser.add_argument('input', type=str, metavar='input'),
    parser.add_argument('label', type=str, metavar='label'),
    parser.add_argument('snpmiss_plot', type=str, metavar='snpmiss_plot'),
    parser.add_argument('cut_geno', type=str, metavar='cut_geno'),
    args = parser.parse_args()
    return args

args = parseArguments()

data = pd.read_csv(args.input,delim_whitespace=True)

fig = plt.figure(figsize=(17,14))
fig,ax = plt.subplots()
matplotlib.rcParams['ytick.labelsize']=13
matplotlib.rcParams['xtick.labelsize']=13
miss = data["F_MISS"]
interesting = interesting = -np.log10(miss[miss>0].sort_values().values)
interesting = np.sort(interesting)
n = np.arange(1,len(interesting)+1) / np.float(len(miss))
ax.step(interesting,n)
plt.axvline(x=-math.log10(float(args.cut_geno)),linewidth=0.5,color='r',linestyle='dashed')
ax.set_xlabel("-log10(Missingness)",fontsize=14)
ax.set_ylabel("Proportion of %s"%args.label,fontsize=14)
ax.set_title("Cumulative prop. of  %s with given missingness"%args.label,fontsize=16)
fig.tight_layout()
plt.savefig(args.snpmiss_plot,format='png')
