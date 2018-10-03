#!/usr/bin/env python3

import pandas as pd
import numpy  as np
import sys
import re

badsex = open(sys.argv[1])
fail_IBD = open(sys.argv[2])
fail_het = open(sys.argv[3])

txt = '''Remove\tDue to\tFiles
%s Individuals\tMale_XHE > %s; Female_XHE < %s\t%s
%s Individuals\tIBD > %s\t%s
%s Individuals\tHET > %s times of mean HET\t%s
'''%(len(badsex.readlines())-1,sys.argv[5],sys.argv[6],sys.argv[1],len(fail_IBD.readlines())-1,sys.argv[7],sys.argv[2],len(fail_het.readlines())-1,sys.argv[8],sys.argv[3])

def load_df(file):
    return(pd.read_csv(file,delim_whitespace=True,usecols=['FID','IID']))

pd.concat([load_df(sys.argv[1]),load_df(sys.argv[2]),load_df(sys.argv[3])],ignore_index=True).sort_values('FID').drop_duplicates().to_csv(sys.argv[4]+'.irem',sep='\t',index=False)

with open(sys.argv[4]+'-QCphase3.out','w') as fp:
    fp.write(txt)
