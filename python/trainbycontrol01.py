# Weiling Li (wel4007@med.cornell.edu)

import collections
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import time

df = pd.read_csv('../csv_files/herberts_meta.txt', sep='\s+',header=None)

#### lwl: final !!
# log2 scale
# train all features for control pbmcs
from sklearn.svm import SVR
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import RobustScaler
control01bigmat_all=pd.read_csv("../csv_files/control01bigmat_all_rmdupgenes.csv")


# training
control01train1=pd.DataFrame()
control01train1['meancov1k']=control01bigmat_all['meancov1k']
control01train1['meancovNDR']=control01bigmat_all['meancovNDR']
control01train1['meanFFT950']=control01bigmat_all['meanFFT950']
control01train1['meanFFTNDR']=control01bigmat_all['meanFFTNDR']
control01train1['meanSL1k']=control01bigmat_all['meanSL1k']
#control01train1['meanshannon1k']=control01bigmat_all['meanshannon1k']
control01train1['meansimpson1k']= control01bigmat_all['meansimpson1k_100_300'] 
control01train1['TPM']=control01bigmat_all['TPM']
control01train1 = control01train1.dropna().reset_index(drop=True)#
X = control01train1.iloc[:,:-1]
#from sklearn.preprocessing import StandardScaler 
# zscore=StandardScaler().fit(X)  
# train_data=zscore.transform(X)  
#test_data=zscore.transform(raw_test_data)
#scaleX=MinMaxScaler().fit_transform(X)
y = np.log2(control01train1['TPM']+1)
#y = y.values.reshape(-1, 1)
#scaley=MinMaxScaler().fit_transform(y)
#regrlog2 = make_pipeline(StandardScaler(), SVR(C=1.0, epsilon=0.2)) #
#regrlog2 = make_pipeline(MinMaxScaler(), SVR(C=1.0, epsilon=0.2))
#regrlog2 = make_pipeline(RobustScaler(), SVR(C=1.0, epsilon=0.2))
regrlog2 = SVR()


regrlog2.fit(X,y)#

# save model
import pickle
#pickle.dump(regrlog2, open('../pretrained_model/svrlog2.pkl', 'wb')) 
#model = pickle.load(open('../pretrained_model/svrlog2.pkl', 'rb'))

# test
for samplename in df[1].values:
    print(samplename)
    
    bigmat_all=pd.read_csv(f"../csv_files/{samplename}bigmat.csv")#slurm/herberts/
    test1=pd.DataFrame()
    test1['gene_id']=bigmat_all['# gene_id']##########
    test1['gene_name']=bigmat_all['gene_name']
    test1['meancov1k']=bigmat_all['meancov1k']
    test1['meancovNDR']=bigmat_all['meancovNDR']
    test1['meanFFT950']=bigmat_all['meanFFT950']
    test1['meanFFTNDR']=bigmat_all['meanFFTNDR']
    test1['meanSL1k']=bigmat_all['meanSL1k']
    # test1['meanshannon1k']=bigmat_all['meanshannon1k']
    test1['meansimpson1k']= bigmat_all['meansimpson1k_100_300']  #
   # test1['TPM']=bigmat_all['TPM']
    test1 = test1.dropna().reset_index(drop=True)
    #test_data=zscore.transform(test1.iloc[:,:-1])
    #X_test_normalized = scaler.transform(test1.iloc[:,2:-1])

    predlog2 = regrlog2.predict(test1.iloc[:,2:])#-1])######(X_test_normalized)##########
    
    outputmat=pd.DataFrame()
    outputmat['gene_id']=test1['gene_id']############
    outputmat['gene_name']=test1['gene_name']
    #outputmat['TPM']=test1['TPM']
    outputmat['Prediction']=predlog2
    #r, p = scipy.stats.pearsonr(np.log2(outputmat['TPM']+1),outputmat['Prediction'])
    #print(r)

    outputmat.to_csv(f'scratch/wel4007/cfNDA/herberts/exp/csv_files/prediction/{samplename}.csv', index=False)#slurm/herberts/
    
    # #gt=np.log2(test1['TPM']+1)
    # #gt = gt.values.reshape(-1, 1)
    # #scalegt=MinMaxScaler().fit_transform(gt)
    # r, p = scipy.stats.pearsonr(np.log2(test1['TPM']+1),predlog2)
    # # r, p = scipy.stats.pearsonr(scalegt.ravel(),predlog2)
    # print(r)
    # l=np.where(test1['TPM'] == 0)[0]
    # # print(l)
    # a=np.delete(test1['TPM'].values, l)
    # b=np.delete(predlog2,l)
    # r, p = scipy.stats.pearsonr(np.log2(a+1),b)
    # print(r)