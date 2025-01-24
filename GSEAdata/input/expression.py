# Weiling Li (wel4007@med.cornell.edu)
# generating an expression matrix for GSEA

import numpy as np
import pandas as pd

def combinedf(path,names):
    c1=pd.read_csv(path+names[0]+'.csv')
    c2=pd.read_csv(path+names[1]+'.csv')
    cc=pd.merge(c1,c2, on=['gene_id','gene_name']) #
    for i in range(2,len(names)):
        c=pd.read_csv(path+names[i]+'.csv')
        cc = pd.merge(cc,c, on=['gene_id','gene_name'])
    
    namelist=['gene_id','gene_name']
    for i in range(0,len(names)):
        namelist.append(names[i])
    cc.columns=namelist

    return cc

pbmctpm=pd.read_csv('../../csv_files/controls34TPM.csv')
pbmclog2tpm=pbmctpm.drop(columns=['gene_id','gene_name'])
pbmclog2tpm = np.log2(pbmclog2tpm+1)
pbmclog2tpm['gene_id']=pbmctpm['gene_id']
pbmclog2tpm['gene_name']=pbmctpm['gene_name']
#pbmclog2tpm

#herbertspath = 'scratch/wel4007/cfNDA/herberts/exp/csv_files/prediction/'
#herbertsnames=pd.read_csv('../../csv_files/herberts_samples.txt',header=None)[0].values
#herberts_df = combinedf(herbertspath, herbertsnames)
#herberts_df.to_csv('../../csv_files/herberts63_prediction.csv', index=False)

herberts_df = pd.read_csv('../../csv_files/herberts63_prediction.csv')

log2tpm_df=pd.merge(pbmclog2tpm, herberts_df, on=['gene_id','gene_name'])#mbc_rm1_df, mbc_rm2_df
log2tpm_df = log2tpm_df .drop('gene_id', axis=1)
log2tpm_df.index=log2tpm_df['gene_name']
log2tpm_df.drop('gene_name', axis=1, inplace=True)

df_min_max_scaled = log2tpm_df.copy() 
  
# apply normalization techniques 
for column in df_min_max_scaled.columns: 
    df_min_max_scaled[column] = (df_min_max_scaled[column] - df_min_max_scaled[column].min()) / (df_min_max_scaled[column].max() - df_min_max_scaled[column].min())     
  
# view normalized data 
#print(df_min_max_scaled)
df_min_max_scaled.insert(loc=0, column='DESCRIPTION', value=['NA']*log2tpm_df.shape[0])
df_min_max_scaled.insert(loc=0, column='NAME', value=df_min_max_scaled.index)
df_min_max_scaled.to_csv('controls34_herberts63.log2TPM.minmax.final.txt',index=False, sep='\t')


