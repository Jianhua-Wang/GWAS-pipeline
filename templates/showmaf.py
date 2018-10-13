#!/usr/bin/env python3

import pandas as pd
import numpy  as np
import sys
import matplotlib
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt

EOL=chr(10)
TABLE=chr(9)

if len(sys.argv)<=1:
    sys.argv="showmaf.py $freq  $base $cut_maf".split()

def rename(x):
    if x.left < 0:
       return "0"
    else:
       return str(x)  


def getTable(frm,txtout):
   frm["MAF bin"]=pd.cut(frm['MAF'],[-0.001,0,0.005,0.01,0.02,0.03,0.04,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50],right=True)
   g1 = frm.groupby('MAF bin').size()
   g1 = g1.rename(rename)
   g1.name = "Num of SNPs"
   return g1.to_csv(txtout,index=True,index_label=True,header=True,sep=TABLE)


def getPic(frm,fname):
    mafs = np.sort(frm['MAF'])
    n = np.arange(1,len(mafs)+1) / np.float(len(mafs))
    fig,ax = plt.subplots()
    ax.step(mafs,n)
    plt.axvline(x=float(sys.argv[3]),linewidth=0.5,color='r',linestyle='dashed')
    matplotlib.rcParams['xtick.labelsize']=13
    matplotlib.rcParams['ytick.labelsize']=13
    ax.set_xlabel("Minor allele frequency",fontsize=14)
    ax.set_ylabel("Proportion of SNPS",fontsize=14)
    plt.title('Cumulative MAF-spectrum on QC-ed data',fontsize=16)
    plt.savefig(fname,format='png')



frm = pd.read_csv(sys.argv[1],delim_whitespace=True)
base = sys.argv[2]
pdfout = "%s.png"%base
txtout = "%s.txt"%base

getPic(frm,pdfout)
curr = getTable(frm,txtout)