#! /usr/bin/env python3


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# from matplotlib.ticker import MaxNLocator
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

pop_df = pd.read_csv(args.pop_file,delim_whitespace=True,names=[1,2,3,4,5,6])
merge_df = top_2.merge(pop_df,left_on=1,right_on=2,how='outer')
merge_df = merge_df.fillna(value='CHIMGEN')

sample_df = merge_df[merge_df[6] == 'CHIMGEN']
ceu_df = merge_df[merge_df[6] == 'CEU']
chb_df = merge_df[merge_df[6] == 'CHB']
jpt_df = merge_df[merge_df[6] == 'JPT']
yri_df = merge_df[merge_df[6] == 'YRI']

fig = plt.figure(figsize=(10,8))
fig,ax = plt.subplots()
matplotlib.rcParams['ytick.labelsize']=10
matplotlib.rcParams['xtick.labelsize']=10
markers = ['o','^','s','p','*']
flatui = ["#00b8a9", "#D81159","#474f85" , "#84B1ED", "#f38181","#95e1d3","#aa96da"]

for i, item in enumerate(merge_df[6].value_counts().index):
    ax.scatter(merge_df[merge_df[6]==item]['2_x'],merge_df[merge_df[6]==item]['3_x'],s=30,label=item,linewidths=1,marker=markers[i],alpha=0.3,c=flatui[i])

ax.set_xlabel("PC1 (variation %4.1f %%)"%pc1,fontsize=14)
ax.set_ylabel("PC2 (variation %4.1f %%)"%pc2,fontsize=14)
# plt.xlim(-0.1, 0.6)
# plt.ylim(-0.6, 0.1)
# plt.xticks([0.0,0.1,0.2,0.3,0.4,0.5])
# plt.yticks([-0.5,-0.4,-0.3,-0.2,-0.1,0.0])
#plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
fig.tight_layout()
fig.savefig(fname=args.output,format='png')