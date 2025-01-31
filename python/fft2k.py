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

fft2kperbase = np.memmap(f"csv_files/results/{samplename}fft2k_window190_bytranscript.txt", mode='w+', shape=(2001,len(transcriptsinfo)), dtype='float32')

for i, (sample_name, sample) in enumerate(samples.items()):
    if sample.sample_id ==samplename:
        j=0
        for row in transcriptsinfo.iterrows():
            strand=row[1][4]
            chrom=row[1][0]
            TSS = row[1][5]
                
            #if strand=='+':
            fft_TSS=[]
            for i in range(TSS-1000,TSS+1000+1):
                fft_TSS.append(sample.get_fft_in_region(chrom=chrom,start=i-95, end=i+95))
    
            if strand=='-':
                fft_TSS=fft_TSS[::-1]
                
            #FFT_TSS[symbol]=fft_TSS
                
            fft2kperbase[:,j]= fft_TSS
    
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

fft2ktransmat = np.memmap(f"csv_files/results/{samplename}fft2k_window190_bytranscript.txt", mode='r', shape=(2001,len(transcriptsinfo)), dtype='float32')
fft2ktransmat = pd.DataFrame(fft2ktransmat)
fft2ktransmat.columns=transcriptsinfo[3].values

fft2kgenemat=pd.DataFrame(columns=genes)

for gene in genes:
    fft2kgenemat[gene]=np.array([fft2ktransmat[col] for col in fft2ktransmat.columns if gene in col]).mean(axis=0)
    

fft2kgenemat.to_csv(f'csv_files/results/{samplename}fft2k_window190_bygene.csv', sep=',', header=True, index=False)

print('gene finished\n')