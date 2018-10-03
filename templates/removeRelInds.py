#!/usr/bin/env python3


import pandas as pd
import sys

if len(sys.argv)<=1:
    sys.argv=["removeRelInds.py","$missing","$ibd_genome","$outfname","$pi_hat"]

outf = sys.argv[3]
TABLE = chr(9)

imissf = pd.read_csv(sys.argv[1],delim_whitespace=True,index_col=["FID","IID"])
genomef = pd.read_csv(sys.argv[2],delim_whitespace=True,usecols=["FID1","IID1","FID2","IID2","PI_HAT"])

pi_hat = float(sys.argv[4])

related_inds = genomef[genomef['PI_HAT']>pi_hat]

fail_IBD = pd.DataFrame(columns=['FID','IID'])
n = 0
for i in related_inds.index:
    if imissf.loc[related_inds.loc[i]['FID1'],related_inds.loc[i]['IID1']]['F_MISS'] >= imissf.loc[related_inds.loc[i]['FID2'],related_inds.loc[i]['IID2']]['F_MISS']:
        fail_IBD.loc[n] = [related_inds.loc[i]['FID1'],related_inds.loc[i]['IID1']]
        n += 1
    else:
        fail_IBD.loc[n] = [related_inds.loc[i]['FID2'],related_inds.loc[i]['IID2']]
        n += 1

fail_IBD.drop_duplicates().sort_values('FID').to_csv(outf,sep=TABLE,index=False)
