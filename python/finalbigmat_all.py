# Weiling Li (wel4007@med.cornell.edu)
# cfDNA features matrix

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

samplename=sys.argv[1] # sample_id

from scripts import cfDNA_calculate_fft_fp_any_region_herberts as funcs
samples = funcs.load_samples()
transcriptsinfo = pd.read_csv("../../../STARref/finalrename.TSS.txt", sep='\t', header=None)

##########################################
def getnewmatrix_all(bigmat):#, geneTPM):
    mynewmatrix=pd.DataFrame(columns=['gene_id','gene_name',                               'meancov1k','meancov950','meancovNDR','meancov380_190','meanFFT950','meanFFTNDR','meanSL1k','meanSL950','meansimpson1k_100_300',
'mediancov1k','mediancov950','mediancovNDR','mediancov380_190','medianFFT950','medianFFTNDR','medianSL1k','medianSL950','mediansimpson1k_100_300',
'mincov1k','mincov950','mincovNDR','mincov380_190','minFFT950','minFFTNDR','maxSL1k','maxSL950','maxsimpson1k_100_300',
                                      'TPM']) #,'Shannon','simpson']) 

#'meanshannon1k','meanshannon950','meanshannon1k_100_300','meanshannon950_100_300','meansimpson1k','meansimpson950',
#'meansimpson950_100_300','meangc1k','meannewshannon1k','meannewshannon1k_100_300','meannewsimpson1k','meannewsimpson1k_100_300',
#'medianshannon1k','medianshannon950','medianshannon1k_100_300','medianshannon950_100_300','mediansimpson1k','mediansimpson950',
#'mediansimpson950_100_300','mediangc1k','mediannewshannon1k','mediannewshannon1k_100_300','mediannewsimpson1k','mediannewsimpson1k_100_300',
#'maxshannon1k','maxshannon950','maxshannon1k_100_300','maxshannon950_100_300','maxsimpson1k','maxsimpson950',
#'maxsimpson950_100_300','maxgc1k','mingc1k','maxnewshannon1k','maxnewshannon1k_100_300','maxnewsimpson1k','maxnewsimpson1k_100_300',
    
    for gene in bigmat['gene_id'].unique():
        gene_name=bigmat[bigmat['gene_id']==gene]['gene_name'].values[0]        
        meancov1k=bigmat[bigmat['gene_id']==gene]['meancov1k'].values.mean()
        meancov950=bigmat[bigmat['gene_id']==gene]['meancov950'].values.mean()
        meancovNDR=bigmat[bigmat['gene_id']==gene]['meancovNDR'].values.mean()
        meancov380_190=bigmat[bigmat['gene_id']==gene]['meancov380_190'].values.mean()
        meanfft950=bigmat[bigmat['gene_id']==gene]['FFT950'].values.mean()
        meanfftndr=bigmat[bigmat['gene_id']==gene]['FFTNDR'].values.mean()
        meansl1k=bigmat[bigmat['gene_id']==gene]['SL1k'].values.mean()
        meansl950=bigmat[bigmat['gene_id']==gene]['SL950'].values.mean()
        # meanshannon1k=bigmat[bigmat['gene_id']==gene]['shannon1k'].values.mean()
        # meanshannon950=bigmat[bigmat['gene_id']==gene]['shannon950'].values.mean()
        # meanshannon1k_100_300=bigmat[bigmat['gene_id']==gene]['shannon1k_100_300'].values.mean()
        # meanshannon950_100_300=bigmat[bigmat['gene_id']==gene]['shannon950_100_300'].values.mean()
        # meansimpson1k=bigmat[bigmat['gene_id']==gene]['simpson1k'].values.mean()
        # meansimpson950=bigmat[bigmat['gene_id']==gene]['simpson950'].values.mean()
        meansimpson1k_100_300=bigmat[bigmat['gene_id']==gene]['simpson1k_100_300'].values.mean()
        # meansimpson950_100_300=bigmat[bigmat['gene_id']==gene]['simpson950_100_300'].values.mean()
        # meangc1k=bigmat[bigmat['gene_id']==gene]['gc1k'].values.mean()
        # meannewshannon1k=bigmat[bigmat['gene_id']==gene]['newshannon1k'].values.mean()
        # meannewshannon1k_100_300=bigmat[bigmat['gene_id']==gene]['newshannon1k_100_300'].values.mean()
        # meannewsimpson1k=bigmat[bigmat['gene_id']==gene]['newsimpson1k'].values.mean()
        # meannewsimpson1k_100_300=bigmat[bigmat['gene_id']==gene]['newsimpson1k_100_300'].values.mean()
        
        mediancov1k=np.median(bigmat[bigmat['gene_id']==gene]['meancov1k'].values)
        mediancov950=np.median(bigmat[bigmat['gene_id']==gene]['meancov950'].values)
        mediancovNDR=np.median(bigmat[bigmat['gene_id']==gene]['meancovNDR'].values)
        mediancov380_190=np.median(bigmat[bigmat['gene_id']==gene]['meancov380_190'].values)
        medianfft950=np.median(bigmat[bigmat['gene_id']==gene]['FFT950'].values)
        medianfftndr=np.median(bigmat[bigmat['gene_id']==gene]['FFTNDR'].values)
        mediansl1k=np.median(bigmat[bigmat['gene_id']==gene]['SL1k'].values)
        mediansl950=np.median(bigmat[bigmat['gene_id']==gene]['SL950'].values)
        # medianshannon1k=np.median(bigmat[bigmat['gene_id']==gene]['shannon1k'].values)
        # medianshannon950=np.median(bigmat[bigmat['gene_id']==gene]['shannon950'].values)
        # medianshannon1k_100_300=np.median(bigmat[bigmat['gene_id']==gene]['shannon1k_100_300'].values)
        # medianshannon950_100_300=np.median(bigmat[bigmat['gene_id']==gene]['shannon950_100_300'].values)
        # mediansimpson1k=np.median(bigmat[bigmat['gene_id']==gene]['simpson1k'].values)
        # mediansimpson950=np.median(bigmat[bigmat['gene_id']==gene]['simpson950'].values)
        mediansimpson1k_100_300=np.median(bigmat[bigmat['gene_id']==gene]['simpson1k_100_300'].values)
        # mediansimpson950_100_300=np.median(bigmat[bigmat['gene_id']==gene]['simpson950_100_300'].values)
        # mediangc1k=np.median(bigmat[bigmat['gene_id']==gene]['gc1k'].values)
        # mediannewshannon1k=np.median(bigmat[bigmat['gene_id']==gene]['newshannon1k'].values)
        # mediannewshannon1k_100_300=np.median(bigmat[bigmat['gene_id']==gene]['newshannon1k_100_300'].values)
        # mediannewsimpson1k=np.median(bigmat[bigmat['gene_id']==gene]['newsimpson1k'].values)
        # mediannewsimpson1k_100_300=np.median(bigmat[bigmat['gene_id']==gene]['newsimpson1k_100_300'].values)

        mincov1k=min(bigmat[bigmat['gene_id']==gene]['meancov1k'].values)
        mincov950=min(bigmat[bigmat['gene_id']==gene]['meancov950'].values)
        mincovNDR=min(bigmat[bigmat['gene_id']==gene]['meancovNDR'].values)
        mincov380_190=min(bigmat[bigmat['gene_id']==gene]['meancov380_190'].values)
        minfft950=min(bigmat[bigmat['gene_id']==gene]['FFT950'].values)
        minfftndr=min(bigmat[bigmat['gene_id']==gene]['FFTNDR'].values)
        maxsl1k=max(bigmat[bigmat['gene_id']==gene]['SL1k'].values)
        maxsl950=max(bigmat[bigmat['gene_id']==gene]['SL950'].values)
        # maxshannon1k=max(bigmat[bigmat['gene_id']==gene]['shannon1k'].values)
        # maxshannon950=max(bigmat[bigmat['gene_id']==gene]['shannon950'].values)
        # maxshannon1k_100_300=max(bigmat[bigmat['gene_id']==gene]['shannon1k_100_300'].values)
        # maxshannon950_100_300=max(bigmat[bigmat['gene_id']==gene]['shannon950_100_300'].values)
        # maxsimpson1k=max(bigmat[bigmat['gene_id']==gene]['simpson1k'].values)
        # maxsimpson950=max(bigmat[bigmat['gene_id']==gene]['simpson950'].values)
        maxsimpson1k_100_300=max(bigmat[bigmat['gene_id']==gene]['simpson1k_100_300'].values)
        # maxsimpson950_100_300=max(bigmat[bigmat['gene_id']==gene]['simpson950_100_300'].values)
        # maxgc1k=max(bigmat[bigmat['gene_id']==gene]['gc1k'].values)
        # mingc1k=min(bigmat[bigmat['gene_id']==gene]['gc1k'].values)
        # maxnewshannon1k=max(bigmat[bigmat['gene_id']==gene]['newshannon1k'].values)
        # maxnewshannon1k_100_300=max(bigmat[bigmat['gene_id']==gene]['newshannon1k_100_300'].values)
        # maxnewsimpson1k=max(bigmat[bigmat['gene_id']==gene]['newsimpson1k'].values)
        # maxnewsimpson1k_100_300=max(bigmat[bigmat['gene_id']==gene]['newsimpson1k_100_300'].values)
        
        #tpm=geneTPM[geneTPM[0]==gene][1].values[0]
       # tpm=geneTPM[geneTPM['# gene_id']==gene]['TPM'].values[0]
        newrow = {'gene_id':gene, 'gene_name': gene_name, 'meancov1k':meancov1k, 'meancov950':meancov950,'meancovNDR':meancovNDR,'meancov380_190':meancov380_190,'meanFFT950':meanfft950,'meanFFTNDR':meanfftndr,'meanSL1k':meansl1k,'meanSL950':meansl950,'meansimpson1k_100_300':meansimpson1k_100_300, 
# 'meanshannon1k':meanshannon1k,'meanshannon950':meanshannon950,'meanshannon1k_100_300':meanshannon1k_100_300,'meanshannon950_100_300':meanshannon950_100_300,'meansimpson1k':meansimpson1k,'meansimpson950':meansimpson950,'meansimpson950_100_300':meansimpson950_100_300,'meangc1k':meangc1k, 'meannewshannon1k':meannewshannon1k,'meannewshannon1k_100_300':meannewshannon1k_100_300,'meannewsimpson1k':meannewsimpson1k,'meannewsimpson1k_100_300':meannewsimpson1k_100_300,
                  'mediancov1k':mediancov1k,'mediancov950':mediancov950,'mediancovNDR':mediancovNDR,'mediancov380_190':mediancov380_190,'medianFFT950':medianfft950,'medianFFTNDR':medianfftndr,'medianSL1k':mediansl1k,'medianSL950':mediansl950,'mediansimpson1k_100_300':mediansimpson1k_100_300,
                  # 'medianshannon1k':medianshannon1k,'medianshannon950':medianshannon950,'medianshannon1k_100_300':medianshannon1k_100_300,'medianshannon950_100_300':medianshannon950_100_300,'mediansimpson1k':mediansimpson1k,'mediansimpson950':mediansimpson950,'mediansimpson950_100_300':mediansimpson950_100_300,'mediangc1k':mediangc1k,'mediannewshannon1k':mediannewshannon1k,'mediannewshannon1k_100_300':mediannewshannon1k_100_300,'mediannewsimpson1k':mediannewsimpson1k,'mediannewsimpson1k_100_300':mediannewsimpson1k_100_300,
                  'mincov1k':mincov1k,'mincov950':mincov950,'mincovNDR':mincovNDR,'mincov380_190':mincov380_190,'minFFT950':minfft950,'minFFTNDR':minfftndr,'maxSL1k':maxsl1k,'maxSL950':maxsl950,'maxsimpson1k_100_300':maxsimpson1k_100_300}#, 
                  # 'maxshannon1k':maxshannon1k,'maxshannon950':maxshannon950,'maxshannon1k_100_300':maxshannon1k_100_300,'maxshannon950_100_300':maxshannon950_100_300,'maxsimpson1k':maxsimpson1k,'maxsimpson950':maxsimpson950,'maxsimpson950_100_300':maxsimpson950_100_300,'maxgc1k':maxgc1k,'mingc1k':mingc1k, 'maxnewshannon1k':maxnewshannon1k,'maxnewshannon1k_100_300':maxnewshannon1k_100_300,'maxnewsimpson1k':maxnewsimpson1k,'maxnewsimpson1k_100_300':maxnewsimpson1k_100_300,
                  #'TPM':tpm}#,'Shannon':shannon, simpson
        new_df = pd.DataFrame([newrow])
        mynewmatrix = pd.concat([mynewmatrix, new_df], axis=0, ignore_index=True)
        
    return mynewmatrix
######################################

meancov1k = np.memmap(f"csv_files/{samplename}meancov1k.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
meancov950 = np.memmap(f"csv_files/{samplename}meancov950.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
meancovNDR= np.memmap(f"csv_files/{samplename}meancovNDR.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32') 
meancov380_190= np.memmap(f"csv_files/{samplename}meancov380_190.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32') 
FFT950 = np.memmap(f"csv_files/{samplename}FFT950.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
FFTNDR = np.memmap(f"csv_files/{samplename}FFTNDR.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
SL1k = np.memmap(f"csv_files/{samplename}SL1k.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
SL950 = np.memmap(f"csv_files/{samplename}SL950.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
# shannon1k=pd.read_csv(f"../notebooks/slurm/{samplename}finalentropy1k.fe.bed",sep='\t')['fragmentation_entropy']
# shannon950=pd.read_csv(f"../notebooks/slurm/{samplename}finalentropy950.fe.bed",sep='\t')['fragmentation_entropy']
# shannon1k_100_300=pd.read_csv(f"csv_files/entropy/{samplename}finalentropy1k.fe.bed",sep='\t')['fragmentation_entropy']
# shannon950_100_300=pd.read_csv(f"csv_files/entropy/{samplename}finalentropy950.fe.bed",sep='\t')['fragmentation_entropy']
# simpson1k=pd.read_csv(f"csv_files/entropy/finalentropy1k.simpson_entropy.{samplename}.bed",sep='\t')['simpson_entropy'] #0-500
# simpson950=pd.read_csv(f"csv_files/entropy/finalentropy950.simpson_entropy.{samplename}.bed",sep='\t')['simpson_entropy']
simpson1k_100_300=pd.read_csv(f"csv_files/finalentropy1k.simpson_entropy.{samplename}.bed",sep='\t')['simpson_entropy']
# simpson950_100_300=pd.read_csv(f"csv_files/simpson_100_300/finalentropy950.simpson_entropy.{samplename}.bed",sep='\t')['simpson_entropy']
# gc1k = np.memmap(f"csv_files/{samplename}gc1k.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
# # method2 for simpson: notebooks/finalfeatures_entropy_gc.py
# newshannon1k = np.memmap(f"csv_files/{samplename}shannon1k.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
# newsimpson1k= np.memmap(f"csv_files/{samplename}simpson1k.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32') 
# newshannon1k_100_300 = np.memmap(f"csv_files/{samplename}shannon1k_100_300.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32')
# newsimpson1k_100_300= np.memmap(f"csv_files/{samplename}simpson1k_100_300.txt", mode='r', shape=(1,len(transcriptsinfo)), dtype='float32') 
# #probability = np.memmap(f"csv_files/{samplename}probability.txt", mode='r', shape=(3,len(transcriptsinfo)), dtype='float32') 

genes=[]
genenames=[]
for tran in transcriptsinfo[3]:
    genes.append(tran.split(":")[0])#gene_id
    genenames.append(tran.split(":")[1])#gene_name

bigmat=pd.DataFrame()
bigmat['gene_id']=genes
bigmat['gene_name']=genenames
# bigmat['transcript_id']=trans
bigmat['meancov1k']=meancov1k.tolist()[0]
bigmat['meancov950']=meancov950.tolist()[0]
bigmat['meancovNDR']=meancovNDR.tolist()[0]
bigmat['meancov380_190']=meancov380_190.tolist()[0]
bigmat['FFT950']=FFT950.tolist()[0]
bigmat['FFTNDR']=FFTNDR.tolist()[0]
bigmat['SL1k']=SL1k.tolist()[0]
bigmat['SL950']=SL950.tolist()[0]
# bigmat['shannon1k']=shannon1k
# bigmat['shannon950']=shannon950
# bigmat['shannon1k_100_300']=shannon1k_100_300
# bigmat['shannon950_100_300']=shannon950_100_300
# bigmat['simpson1k']=simpson1k
# bigmat['simpson950']=simpson950
bigmat['simpson1k_100_300']=simpson1k_100_300
# bigmat['simpson950_100_300']=simpson950_100_300
# bigmat['gc1k']=gc1k.tolist()[0]
# bigmat['newshannon1k']=newshannon1k.tolist()[0]
# bigmat['newsimpson1k']=newsimpson1k.tolist()[0]
# bigmat['newshannon1k_100_300']=newshannon1k_100_300.tolist()[0]
# bigmat['newsimpson1k_100_300']=newsimpson1k_100_300.tolist()[0]
# bigmat

#geneTPM = pd.read_csv(f"../proteincodingTPM/csv/{samplename}.csv",sep=',')#rnaseqname ######################
bigmat_all = getnewmatrix_all(bigmat)#, geneTPM)
#bigmat_all
np.savetxt(f"csv_files/{samplename}bigmat.csv", bigmat_all, delimiter=',',header=','.join(bigmat_all.columns.values),fmt='%s')
# only for PM1070.1358 metastatic sites
# np.savetxt(f"csv_files/{rnaseqname}bigmat_all.csv", bigmat_all, delimiter=',',header=','.join(bigmat_all.columns.values),fmt='%s')
