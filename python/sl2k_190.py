# Weiling Li (wel4007@med.cornell.edu)
# generating feature values for each bp in 2k-TSS regions

import sys
from modules import misc_tools

import collections
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import methodtools
import more_itertools
import numpy as np
import pandas as pd
import pathlib
import scipy.spatial.distance
import scipy.signal
import time
import sys

samplename=sys.argv[1]

from scripts import cfDNA_calculate_fft_fp_any_region_herberts as funcs
samples = funcs.load_samples()
transcriptsinfo = pd.read_csv("../csv_files/finalrename.TSS.txt", sep='\t', header=None)

sl2kperbase = np.memmap(f"csv_files/results/{samplename}sl2k_window190_bytranscript.txt", mode='w+', shape=(2001,len(transcriptsinfo)), dtype='float32')

for i, (sample_name, sample) in enumerate(samples.items()):
    if sample.sample_id ==samplename:
        j=0
        for row in transcriptsinfo.iterrows():
            strand=row[1][4]
            chrom=row[1][0]
            TSS = row[1][5]
                
            #if strand=='+':
            sl_TSS=[]
            for i in range(TSS-1000,TSS+1000+1):
                sl_TSS.append(sample.get_fp_in_region(chrom=chrom,start=i-95, end=i+95))
    
            if strand=='-':
                sl_TSS=sl_TSS[::-1]
                
            #FFT_TSS[symbol]=sl_TSS
                
            sl2kperbase[:,j]= sl_TSS
    
            j=j+1

print('transcript finished\n')            
genes=[]
genenames=[]
for tran in transcriptsinfo[3]:
    if tran.split(":")[0] not in genes:
        genes.append(tran.split(":")[0])#gene_id
    if tran.split(":")[1] not in genenames:
        genenames.append(tran.split(":")[1])
# print(len(geneshg38)) #44420
# print(len(np.unique(geneshg38))) #17395

sl2ktransmat = np.memmap(f"csv_files/results/{samplename}sl2k_window190_bytranscript.txt", mode='r', shape=(2001,len(transcriptsinfo)), dtype='float32')
sl2ktransmat = pd.DataFrame(sl2ktransmat)
sl2ktransmat.columns=transcriptsinfo[3].values

sl2kgenemat=pd.DataFrame(columns=genes)

for gene in genes:
    sl2kgenemat[gene]=np.array([sl2ktransmat[col] for col in sl2ktransmat.columns if gene in col]).mean(axis=0)
    

sl2kgenemat.to_csv(f'csv_files/results/{samplename}sl2k_window190_bygene.csv', sep=',', header=True, index=False)

print('gene finished\n')