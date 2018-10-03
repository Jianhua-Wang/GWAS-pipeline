#!/usr/bin/env python3

import pandas as pd
import numpy  as np
import sys
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import matplotlib
EOL=chr(10)


if len(sys.argv)<=1:
    sys.argv="showmaf.py $hwe  $base".split()

# We don't draw a Manhatten plot at the moment -- would be interesting to see if there are regions
# that are problematic -- problem is that the .hwe file doesn't have base position so we'd need to change
# the pipeline but it would be good to do
# def drawManhatten(pheno, result,outf):
#     fig, (ax1, ax2) =  plt.subplots(2, 1, sharey=True)
#     ax=ax1
#     delta=0
#     colours= ['crimson','blue','green']
#     xtick_pos=[]
#     xtick_label = []
#     chroms = result.groupby("CHR")
#     for chrom_num, chrom_res in chroms:
#         this_chrom = result['CHR']==chrom_num
#         result.loc[this_chrom,'BP']+=delta
#         old_delta=delta
#         delta = delta + int(chrom_res.tail(1)['BP'])
#         xtick_pos.append((delta+old_delta)/2)
#         xtick_label.append(str(chrom_num))
#         under_thresh = result['P']<0.005
#         ax.scatter(result.loc[this_chrom & under_thresh, 'BP'],\
#                    -np.log10(result.loc[this_chrom  & under_thresh,'P']),c=colours[chrom_num%3])
#         if chrom_num == 9:
#            ax.set_xticklabels(xtick_label)
#            ax.set_xticks(xtick_pos)
#            xtick_pos=[]
#            xtick_label=[]
#            ax=ax2
#            delta=0
#     ax.set_xticklabels(xtick_label)
#     ax.set_xticks(xtick_pos)
#     plt.savefig(outf)

def drawQQ(result,outf):
    plt.figure()
    fig, ax = plt.subplots()
    sort_p = -np.log10(result['P'].sort_values())
    n=len(sort_p)
    expected = -np.log10(np.linspace(1/n,1,n))
    ax.set_xlabel("HWE expected (-log p)",fontsize=14)
    ax.set_ylabel("HWE observed (-log p)",fontsize=14)
    plt.plot(expected,sort_p)
    plt.plot(expected,expected)
    plt.savefig(outf,format='png')

def getPic(frm,test,pdfout):
   fig = plt.figure(figsize=(17,14))
   fig,ax = plt.subplots()
   matplotlib.rcParams['ytick.labelsize']=13
   matplotlib.rcParams['xtick.labelsize']=13
   hwe = frm[frm["TEST"]==test]["P"]
   big = min(hwe.mean()+2*hwe.std(),hwe.nlargest(4).iloc[3])
   if big > max(0.95*len(hwe),100):
       hwe = hwe[hwe<big]
   hwe = np.sort(hwe)
   n = np.arange(1,len(hwe)+1) / np.float(len(hwe))
   ax.step(hwe,n)
   ax.set_xlabel("HWE p score",fontsize=14)
   ax.set_ylabel("Proportion of SNPs with HWE p-value or less",fontsize=14)
   ax.set_title("Cumulative prop. of SNPs with HWE or less",fontsize=16)
   fig.tight_layout()
   plt.savefig(pdfout,format='png')


f = open(sys.argv[1])
header=f.readline()
test=False
for i in range(5):
    line=f.readline()
    if "ALL(QT)" in line:
        test="ALL(QT)"
    elif "ALL(NP)" in line:
        test="ALL(NP)"
    elif "ALL" in line:
        test="ALL"
if not test:
    print((EOL*5)+"The Hardy-Weinberg file is malformed, can't find ALL test <%s>"%sys.argv[1])
    sys.exit(12)

frm = pd.read_csv(sys.argv[1],delim_whitespace=True)
frm = frm[frm["TEST"]==test]
base = sys.argv[2]
pdfout = "%s.png"%base
qqpdf  = "%s-qq.png"%base
getPic(frm,test,pdfout)
drawQQ(frm, qqpdf)