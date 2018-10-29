#! /usr/bin/env python3


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np
import argparse
import sys

EOL=chr(10)
gap = EOL*3

def parseArguments():
   if len(sys.argv)<=1:
      sys.argv=\
      "drawPCA.py $base $eigvals $eigvecs $output $pop_file".split()
   parser=argparse.ArgumentParser()
   parser.add_argument('input', type=str, metavar='input'),
   parser.add_argument('eigvals', type=str, metavar='label'),
   parser.add_argument('eigvecs', type=str, metavar='output'),
   parser.add_argument('output', type=str, metavar='output'),
   parser.add_argument('pop_file',type=str, metavar='pop_file')
   args = parser.parse_args()
   return args
args = parseArguments()

evals = np.loadtxt(args.eigvals)
pc1 = 100*evals[0]/evals[:10].sum()
pc2 = 100*evals[1]/evals[:10].sum()

df = pd.read_csv(args.eigvecs,delim_whitespace=True,names=list(range(22)))
top_2 = df[[0,1,2,3]]

pop_df = pd.read_csv(args.pop_file,delim_whitespace=True)
merge_df = top_2.merge(pop_df,left_on=1,right_on='id',how='outer')

merge_df = merge_df.fillna(value='sample')
sample_df = merge_df[merge_df['pop'] == 'sample']
chs_df = merge_df[merge_df['pop'] == 'CHS']
chb_df = merge_df[merge_df['pop'] == 'CHB']

fig = plt.figure(figsize=(10,8))
fig,ax = plt.subplots()
matplotlib.rcParams['ytick.labelsize']=10
matplotlib.rcParams['xtick.labelsize']=10
ax.scatter(chb_df[2],chb_df[3],s=20,label='CHB',facecolors='none',edgecolors='#00b8a9',linewidths=0.6)
ax.scatter(chs_df[2],chs_df[3],s=20,label='CHS',facecolors='none',edgecolors='#ffde7d',linewidths=0.6)
ax.scatter(sample_df[2],sample_df[3],s=20,label='Sample',facecolors='none',edgecolors='#f6416c',linewidths=0.6)
ax.set_xlabel("PC1 (variation %4.1f %%)"%pc1,fontsize=14)
ax.set_ylabel("PC2 (variation %4.1f %%)"%pc2,fontsize=14)
plt.xlim(-0.1, 0.6)
plt.ylim(-0.6, 0.1)
plt.xticks([0.0,0.1,0.2,0.3,0.4,0.5])
plt.yticks([-0.5,-0.4,-0.3,-0.2,-0.1,0.0])
#plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
fig.tight_layout()
fig.savefig(fname=args.output,format='png')