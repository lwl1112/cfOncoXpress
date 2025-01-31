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

cov2kperbase = np.memmap(f"csv_files/results/{samplename}cov2k_bytranscript.txt", mode='w+', shape=(2001,len(transcriptsinfo)), dtype='float32')

for i, (sample_name, sample) in enumerate(samples.items()):
    if sample.sample_id ==samplename:
        j=0
        for row in transcriptsinfo.iterrows():
            strand=row[1][4]
            chrom=row[1][0]
            TSS = row[1][5]
            
            if strand=='+':
                cov1k=sample.get_coverage_in_region(chrom=chrom, start=TSS-1000, end=TSS+1000+1) #-1k~1k
                # cov285=sample.get_coverage_in_region(chrom=chrom, start=TSS-285, end=TSS+285+1)
                # covndr=sample.get_coverage_in_region(chrom=chrom, start=TSS-150, end=TSS+50+1) # -150 ~ 50
            else:
                cov1k=sample.get_coverage_in_region(chrom=chrom, start=TSS-1000, end=TSS+1000+1)[::-1] #+1k~-1k
                # cov285=sample.get_coverage_in_region(chrom=chrom, start=TSS-285, end=TSS+285+1)[::-1] 
                # covndr=sample.get_coverage_in_region(chrom=chrom, start=TSS-50, end=TSS+150+1)[::-1] # +150 ~ -50

            cov2kperbase[:,j]= cov1k
            # cov570perbase[:,j]= cov285
            # covNDRperbase[:,j]= covndr
            # cov[:,j]= np.mean(cov1k)
            # gc[:,j]=sample.get_gc_content(chrom=chrom, start=TSS-1000, end=TSS+1000+1)
            j=j+1

print('transcript finished\n')            
genes=[]
# genenames=[]
for tran in transcriptsinfo[3]:
    if tran.split(":")[0] not in genes:
        genes.append(tran.split(":")[0])#gene_id
# print(len(geneshg38)) #44420
# print(len(np.unique(geneshg38))) #17395

cov2ktransmat = np.memmap(f"csv_files/results/{samplename}cov2k_bytranscript.txt", mode='r', shape=(2001,len(transcriptsinfo)), dtype='float32')
cov2ktransmat = pd.DataFrame(cov2ktransmat)
cov2ktransmat.columns=transcriptsinfo[3].values

# cov570transmat = np.memmap(f"../csv_files/finalfeatures/redo/2k/{samplename}cov570_bytranscript.txt", mode='r', shape=(571,len(transcriptsinfo)), dtype='float32')
# cov570transmat = pd.DataFrame(cov570transmat)
# cov570transmat.columns=transcriptsinfo[3].values

# covNDRtransmat = np.memmap(f"../csv_files/finalfeatures/redo/2k/{samplename}covNDR_bytranscript.txt", mode='r', shape=(201,len(transcriptsinfo)), dtype='float32')
# covNDRtransmat = pd.DataFrame(covNDRtransmat)
# covNDRtransmat.columns=transcriptsinfo[3].values

# meancovtransmat = np.memmap(f"../csv_files/finalfeatures/redo/2k/{samplename}meancov_bytranscript.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
# meancovtransmat = pd.DataFrame(meancovtransmat)
# meancovtransmat.columns=transcriptsinfo[3].values

# gctransmat = np.memmap(f"../csv_files/finalfeatures/redo/2k/{samplename}gc_bytranscript.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
# gctransmat = pd.DataFrame(gctransmat)
# gctransmat.columns=transcriptsinfo[3].values

cov2kgenemat=pd.DataFrame(columns=genes)
# cov570genemat=pd.DataFrame(columns=genes)
# covNDRgenemat=pd.DataFrame(columns=genes)
# meancovgenemat=pd.DataFrame(columns=genes)
# gcgenemat=pd.DataFrame(columns=genes)
for gene in genes:
    cov2kgenemat[gene]=np.array([cov2ktransmat[col] for col in cov2ktransmat.columns if gene in col]).mean(axis=0)
    # cov570genemat[gene]=np.array([cov570transmat[col] for col in cov570transmat.columns if gene in col]).mean(axis=0)
    # covNDRgenemat[gene]=np.array([covNDRtransmat[col] for col in covNDRtransmat.columns if gene in col]).mean(axis=0)
    # meancovgenemat[gene]=np.array([meancovtransmat[col] for col in meancovtransmat.columns if gene in col]).mean(axis=0)
    # gcgenemat[gene]=np.array([gctransmat[col] for col in gctransmat.columns if gene in col]).mean(axis=0)
    
# np.savetxt(f"../csv_files/finalfeatures/redo/cov2k/{samplename}cov2k_bygene.csv", covgenemat, delimiter=',',header= ','.join(genemat.columns.values),fmt='%s')
cov2kgenemat.to_csv(f'csv_files/results/{samplename}cov2k_bygene.csv', sep=',', header=True, index=False)
# cov570genemat.to_csv(f'../csv_files/finalfeatures/redo/2k/{samplename}cov570_bygene.csv', sep=',', header=True, index=False)
# covNDRgenemat.to_csv(f'../csv_files/finalfeatures/redo/2k/{samplename}covNDR_bygene.csv', sep=',', header=True, index=False)
# meancovgenemat.to_csv(f'../csv_files/finalfeatures/redo/2k/{samplename}meancov_bygene.csv', sep=',', header=True, index=False)
# gcgenemat.to_csv(f'../csv_files/finalfeatures/redo/2k/{samplename}gc_bygene.csv', sep=',', header=True, index=False)

print('gene finished\n')