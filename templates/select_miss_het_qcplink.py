#!/usr/bin/env python3

# This script selects those who fail heterozygosity constraints -- previous versions also looked at missingness too
# hence the name


from pandas import read_csv
import sys

TAB = chr(9)
if len(sys.argv)<=1:
  sys.argv = ["select_miss_het.qcplink.py","$het","$outfname","$times_of_meanhet"]

hetf = sys.argv[1]

outfname = sys.argv[2]
times_of_meanhet = sys.argv[3]
times_of_meanhet = float(times_of_meanhet)
het      = read_csv(hetf,delim_whitespace=True,dtype={'IID':str})
mean_het = (het["N(NM)"]-het["O(HOM)"])/het["N(NM)"]
cut_het_high = mean_het.mean()+times_of_meanhet*mean_het.std()
cut_het_low = mean_het.mean()-times_of_meanhet*mean_het.std()
failed   = het[(mean_het<cut_het_low) | (mean_het>cut_het_high)]

failed.to_csv(outfname,columns=['FID','IID'],index=False,sep=TAB)
