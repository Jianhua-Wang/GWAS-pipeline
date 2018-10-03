#!/usr/bin/env python3

import pandas as pd
import numpy  as np
import sys
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt

EOL=chr(10)
TABLE=chr(9)

if len(sys.argv)<=1:
    sys.argv="showmaf.py $freq  $base".split()

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
    r = [-0.00001]+list(map(lambda x: x/100,range(0,51,2)))
    xs = list(map(lambda x:x/100,range(0,51,2)))
    frm['bin']=pd.cut(frm['MAF'],r,right=True)
    g2  = frm.groupby('bin').size()
    cum = (g2.sum()-g2.cumsum()+g2.iloc[0])/g2.sum()
    plt.plot(xs,cum)
    plt.xlim(0,0.5)
    plt.ylabel('Proportion SNPs ##*-geq## this freq'.replace("*-",chr(92)).replace("##",chr(36)))
    plt.xlabel('Frequency')
    plt.savefig(fname,format='png')



frm = pd.read_csv(sys.argv[1],delim_whitespace=True)
base = sys.argv[2]
pdfout = "%s.png"%base
txtout = "%s.txt"%base

getPic(frm,pdfout)
curr = getTable(frm,txtout)