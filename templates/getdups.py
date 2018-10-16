#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd
from subprocess import call

def parseArguments():
    if len(sys.argv)<=1:
        sys.argv="getdups.py $base $outfname".split()
    parser=argparse.ArgumentParser()
    parser.add_argument('base', type=str, metavar='base'),
    parser.add_argument('outfname', type=str, metavar='outfname'),
    args = parser.parse_args()
    return args

args = parseArguments()

bim_path = '%s.bim'%args.base
lmiss_path = '%s.lmiss'%args.base
plink_file = args.base

call('plink --keep-allele-order --bfile ' +plink_file+ ' --missing --out ' +plink_file, shell=True)

bim_df = pd.read_csv(bim_path,delim_whitespace=True,names=['chr','id','3','bp','allele1','allele2'])
lmiss_df = pd.read_csv(lmiss_path,delim_whitespace=True)
bim_df['F_MISS'] = lmiss_df['F_MISS']

def get_dup_snps_per_chr(CHR):
    dup_snps_per_chr = pd.DataFrame(columns=['id'])
    chr1_df = bim_df[bim_df['chr'] == CHR]
    bp_counts = chr1_df['bp'].value_counts()
    dup_bp = bp_counts[bp_counts>1].sort_index()

    for i in dup_bp.index:
        dup_df = chr1_df[chr1_df['bp'] == i]
        min_idx = dup_df['F_MISS'].idxmin()
        rm_id = dup_df.drop(labels = min_idx)
        rm_id = pd.DataFrame(rm_id['id'])
        dup_snps_per_chr = dup_snps_per_chr.append(rm_id,ignore_index=True)
    return dup_snps_per_chr

dup_snp_df = pd.DataFrame(columns=['id'])

for CHR in bim_df['chr'].value_counts().sort_index().index:
    dup_snp_df = dup_snp_df.append(get_dup_snps_per_chr(CHR),ignore_index=True)

dup_snp_df.to_csv(args.outfname,index=False,header=False)