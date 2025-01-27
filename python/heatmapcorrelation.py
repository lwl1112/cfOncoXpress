# Weiling Li (wel4007@med.cornell.edu)
# correlations matrix for heatmap

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
tpms = pd.read_csv('../csv_files/rna.csv',sep='\t') #slurm/herberts
#names=pd.read_csv('/athena/khuranalab/scratch/wel4007/notebooks/slurm/herberts/minimeta.txt',sep=' ',header=None)[1]

names=['DTB-205-Baseline-cfDNA','DTB-097-Progression-cfDNA','DTB-102-Progression-cfDNA','DTB-119-Progression-cfDNA','DTB-149-Baseline-cfDNA','DTB-156-Baseline-cfDNA','DTB-183-Baseline-cfDNA','DTB-210-Baseline-cfDNA','DTB-214-Baseline-cfDNA','DTB-216-Progression-cfDNA','DTB-258-Baseline-cfDNA','DTB-261-Baseline-cfDNA','DTB-266-Baseline-cfDNA']

def getnewmatrix2(predfile, tpmfile):#   x:tpm, y:pred
    tpmmat= tpms[['IDENTIFIER',tpmfile]] #pd.read_csv(f"../csv_files/finalfeatures/{tpmfile}bigmat_all.csv")
    tpmmat.columns=['gene_name','TPM']
    Predmat=pd.read_csv(f'scratch/csv/prediction/{predfile}.csv')#(f"../csv_files/finalfeatures/{predfile}_prediction.csv")
    mynewmatrix=pd.DataFrame(columns=['gene_name', 'TPM','prediction']) 
#     mynewmatrix_TPMremove0=pd.DataFrame(columns=['geneName','PFE','OCF_Ratio','NDR_relative_2k','inferredGEP','prediction', 'TPM']) 
    for index, row in tpmmat.iterrows():
    #     if index>1:
    #         break
        genename=row['gene_name']
        tpm=row['TPM']

        s=Predmat[Predmat['gene_name']==genename]
        if s.shape[0]>0:
            for row2 in s.iterrows():
                pred=row2[1]['Prediction'] ####
                newrow = {'gene_name':genename, 'TPM':tpm, 'prediction':pred}
                new_df = pd.DataFrame([newrow])
                mynewmatrix = pd.concat([mynewmatrix, new_df], axis=0, ignore_index=True)
        
    return mynewmatrix# mynewmatrix_TPMremove0


corrmat = []
for i in names: #tissueTPM:
    corlist=[]
    for j in names:#[:-1]: #inferred: #our model prediction
        newmat = getnewmatrix2(i, j)
        newmat = newmat.dropna().reset_index(drop=True)
        r, p = scipy.stats.pearsonr(np.log2(newmat['TPM'].astype(float)+1),newmat['prediction'])
        corlist.append(r) # same tissue, a list of pred
    corrmat.append(corlist)
np.savetxt('../csv_files/results/correlation.final.csv',corrmat,delimiter=",")

