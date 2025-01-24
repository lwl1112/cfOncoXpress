# Weiling Li (wel4007@med.cornell.edu)
# epicseq correlations matrix for heatmap

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

from sklearn.preprocessing import scale
from pandas import DataFrame

tpms = pd.read_csv('rna.csv',sep='\t') #slurm/herberts
#names=pd.read_csv('/athena/khuranalab/scratch/wel4007/notebooks/slurm/herberts/minimeta.txt',sep=' ',header=None)[1]

names=['DTB-205-Baseline-cfDNA','DTB-097-Progression-cfDNA','DTB-102-Progression-cfDNA','DTB-119-Progression-cfDNA','DTB-149-Baseline-cfDNA','DTB-156-Baseline-cfDNA','DTB-183-Baseline-cfDNA','DTB-210-Baseline-cfDNA','DTB-214-Baseline-cfDNA','DTB-216-Progression-cfDNA','DTB-258-Baseline-cfDNA','DTB-261-Baseline-cfDNA','DTB-266-Baseline-cfDNA']


#tmat = pd.read_csv(f'scratch/epicseq/groupby10/{names[0]}.withTPM.csv',sep='\t')
#tmat = tmat.dropna().reset_index(drop=True)
#tmat.to_csv('testprev.csv')

def getnewmatrix2_addTPM(epicseqfile, tpmfile):
    tpmmat= tpms[['IDENTIFIER',tpmfile]] #pd.read_csv(f"../csv_files/finalfeatures/{tpmfile}bigmat_all.csv")
    tpmmat.columns=['gene_name','TPM']
    #mynewmatrix=pd.DataFrame(columns=['geneName','PFE','OCF_Ratio','NDR_relative_2k','inferredGEP','prediction', 'TPM']) 
    epicseqmat=pd.read_csv(f'scratch/epicseq/groupby10/{epicseqfile}.seprarategene.csv',sep='\t')#slurm/herberts/
    mynewmatrix=pd.DataFrame(columns=['geneName','inferredGEP', 'TPM']) 
   # mynewmatrix_TPMremove0=pd.DataFrame(columns=['geneName','PFE','OCF_Ratio','NDR_relative_2k','inferredGEP','prediction', 'TPM']) 
    for index, row in epicseqmat.iterrows():
    #     if index>1:
    #         break
        genename=row['geneName']
        #pfe=row['PFE']
        #ocf=row['OCF_Ratio']
        #ndr=row['NDR_relative_2k']
        inferred=row['inferredGEP']
        #pred = row['prediction']
        
        s=tpmmat[tpmmat['gene_name']==genename]
        if s.shape[0]>0: #(not(np.isnan(inferred))) and (s.shape[0]>0):
            for row in s.iterrows():
                tpm=row[1]['TPM'] ####
                newrow = {'geneName':genename, 'inferredGEP':inferred, 'TPM':tpm}
                #'PFE':pfe, 'OCF_Ratio':ocf, 'NDR_relative_2k':ndr,, 'TPM':tpm
                new_df = pd.DataFrame([newrow])
                mynewmatrix = pd.concat([mynewmatrix, new_df], axis=0, ignore_index=True)
#                 if tpm>0:
#                     mynewmatrix_TPMremove0 = pd.concat([mynewmatrix_TPMremove0, new_df], axis=0, ignore_index=True)
        
    return mynewmatrix#, mynewmatrix_TPMremove0

#newmatPM1940_PM1640TPM, newmatPM1940_PM1640TPMremove0=getnewmatrix2_addTPM(matPM1940, PM1640TPM)
corrmat = []
for i in names: #epicseqGEP:
    corlist=[]
    for j in names:#[:-1]:#tissueTPM:
        newmat = getnewmatrix2_addTPM(i, j)
        #print(newmat)
        newmat = newmat.dropna().reset_index(drop=True)
        #print(newmat)
        #newmat.to_csv('testnow.csv')
        r, p = scipy.stats.pearsonr(np.log2(newmat['TPM'].astype(float)+1),newmat['inferredGEP'])
        corlist.append(r)
        #print(r)
    corrmat.append(corlist)


np.savetxt('../csv_files/epicseqcorrelation.final.csv',corrmat,delimiter=",")


