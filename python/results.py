# Weiling Li (wel4007@med.cornell.edu)
# correlation results for 13 samples with matched RNA-seq

import collections
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import scipy.spatial.distance
import scipy.signal
import time
import sys


def resmat(predmat,truthmat,samplename):
    mynewmatrix=pd.DataFrame(columns=['gene_id','gene_name', 'Prediction', 'TPM'])
    for i,row in predmat.iterrows():
        res=truthmat[(truthmat['IDENTIFIER']==row['gene_name']) & (truthmat['GENE_ID']==row['gene_id'])][samplename]
        if len(res)==1:
            tpm=res.values[0]
            newrow = {'gene_id':row['gene_id'], 'gene_name': row['gene_name'],'Prediction':row['Prediction'], 'TPM':tpm}
            new_df = pd.DataFrame([newrow])
            mynewmatrix = pd.concat([mynewmatrix, new_df], axis=0, ignore_index=True)        
    return mynewmatrix

mat=pd.read_csv('../csv_files/rna.csv',sep='\t')

outdir=f'/scratch/wel4007/cfDNA/herberts/csv_files/prediction' # 
for samplename in mat.columns[1:-1].values:
    print(samplename)
    pred=pd.read_csv(f'{outdir}/{samplename}.csv')[['gene_id','gene_name','Prediction']]
    newmat=resmat(pred,mat,samplename)
    newmat.to_csv(f'{outdir}/{samplename}_res.csv')
    r, p = scipy.stats.pearsonr(np.log2(newmat['TPM']+1),newmat['Prediction'])
    print(r,p)
    