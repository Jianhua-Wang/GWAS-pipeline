#!/usr/bin/env python3

from __future__ import print_function


import sys
import re

f = open(sys.argv[1])
rem_name=sys.argv[2]
rem_inds = rem_snps = rem_maf = rem_hwe =  0

autosome=False
nonautosome = 0
for line in f:
    if "Options in effect" in line: 
        autosome=False
        continue
    if "--autosome" in line:
        autosome=True
        continue
    m =re.search("(\w+) out of (\w+) variants loaded from .bim file.", line)
    if m and autosome:
        nonautosome = int(m.group(2))-int(m.group(1))
        continue
    m=re.search("(\w+) (\w+) removed due to ([-A-z0-9]+ \w+)",line)
    if m:
        #print(m.group(1),m.group(2),m.group(3))
        n = int(m.group(1))
        if m.group(2)=="variants":
            if  m.group(3)=="missing genotype":
                rem_snps = rem_snps+n
            elif m.group(3)=="Hardy-Weinberg exact":
                rem_hwe = rem_hwe+n
            else:
                rem_maf = rem_maf + n
        else:
            rem_inds = rem_inds + n

txt = '''Remove\tDue to
%s Individuals(see the -nd-c.irem file)\t--mind %s
%s Markers\t--geno %s
%s Markers\t--maf %s
%s Markers\t--hwe %s
'''%(rem_inds,sys.argv[3],rem_snps,sys.argv[4],rem_maf,sys.argv[5],rem_hwe,sys.argv[6])

with open(sys.argv[7]+'-QCphase2.out','w') as fp:
    fp.write(txt)
