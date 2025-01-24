# Weiling Li (wel4007@med.cornell.edu)
# generating cfDNA features except simpson entropy

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


samplename=sys.argv[1] # sample_id

from scripts import cfDNA_calculate_fft_fp_any_region_herberts as funcs
samples = funcs.load_samples()
transcriptsinfo = pd.read_csv("../../../STARref/finalrename.TSS.txt", sep='\t', header=None)

meancov1k = np.memmap(f"csv_files/{samplename}meancov1k.txt", mode='w+', shape=(1,len(transcriptsinfo)), dtype='float32')
meancov950 = np.memmap(f"csv_files/{samplename}meancov950.txt", mode='w+', shape=(1,len(transcriptsinfo)), dtype='float32')
meancovNDR= np.memmap(f"csv_files/{samplename}meancovNDR.txt", mode='w+', shape=(1,len(transcriptsinfo)), dtype='float32') 
meancov380_190= np.memmap(f"csv_files/{samplename}meancov380_190.txt", mode='w+', shape=(1,len(transcriptsinfo)), dtype='float32') 
FFT950 = np.memmap(f"csv_files/{samplename}FFT950.txt", mode='w+', shape=(1,len(transcriptsinfo)), dtype='float32')
FFTNDR = np.memmap(f"csv_files/{samplename}FFTNDR.txt", mode='w+', shape=(1,len(transcriptsinfo)), dtype='float32')
SL1k = np.memmap(f"csv_files/{samplename}SL1k.txt", mode='w+', shape=(1,len(transcriptsinfo)), dtype='float32')
SL950 = np.memmap(f"csv_files/{samplename}SL950.txt", mode='w+', shape=(1,len(transcriptsinfo)), dtype='float32')

for i, (sample_id, sample) in enumerate(samples.items()):
    if sample.sample_id ==samplename:
        j=0
        for row in transcriptsinfo.iterrows():
            strand=row[1][4]
            chrom=row[1][0]
            TSS = row[1][5]
#            print(strand+'\t'+chrom+'\t'+str(TSS))
            cov1k=[]
            cov950=[]
            covNDR=[]
            cov380_190=[]
            fft950=[]
            fftNDR=[]
            sl1k=[]
            sl950=[]
            
            if strand=='+':
                cov1k=sample.get_coverage_in_region(chrom=chrom, start=TSS-1000, end=TSS+1000+1) #-1k~1k
                cov950=sample.get_coverage_in_region(chrom=chrom, start=TSS-950, end=TSS+950) # -950~949
                covNDR=sample.get_coverage_in_region(chrom=chrom, start=TSS-150, end=TSS+50+1) # -150~50
                cov380_190=sample.get_coverage_in_region(chrom=chrom, start=TSS-380, end=TSS+190) # -380~189
                fft950=sample.get_fft_in_region(chrom=chrom,start=TSS-950, end=TSS+950) # -950~949
                fftNDR=sample.get_fft_in_region(chrom=chrom,start=TSS-380, end=TSS+190) # -380~189
                sl1k=sample.get_fp_in_region(chrom=chrom,start=TSS-1000, end=TSS+1000+1) #-1k~1k
                sl950=sample.get_fp_in_region(chrom=chrom,start=TSS-950, end=TSS+950) # -950~949
            else:
                cov1k=sample.get_coverage_in_region(chrom=chrom, start=TSS-1000, end=TSS+1000+1)#-1k~1k
                cov950=sample.get_coverage_in_region(chrom=chrom, start=TSS-949, end=TSS+950+1) #-949~950
                covNDR=sample.get_coverage_in_region(chrom=chrom, start=TSS-50, end=TSS+150+1) # -50~150
                cov380_190=sample.get_coverage_in_region(chrom=chrom, start=TSS-189, end=TSS+381) # -189~380               
                fft950=sample.get_fft_in_region(chrom=chrom,start=TSS-949, end=TSS+951)   #-949~950
                fftNDR=sample.get_fft_in_region(chrom=chrom,start=TSS-189, end=TSS+381)   #-189~380
                sl1k=sample.get_fp_in_region(chrom=chrom,start=TSS-1000, end=TSS+1000+1)   #-1k~1k
                sl950=sample.get_fp_in_region(chrom=chrom,start=TSS-949, end=TSS+951)    #-949~950
                    
            meancov1k[:,j]= np.mean(cov1k)####PM
            meancov950[:,j]= np.mean(cov950)####PM
            meancovNDR[:,j]= np.mean(covNDR)####PM
            meancov380_190[:,j]= np.mean(cov380_190)####PM
            FFT950[:,j]= fft950####PM
            FFTNDR[:,j]= fftNDR####PM
            SL1k[:,j]= sl1k####PM
            SL950[:,j]= sl950####PM
            j=j+1  