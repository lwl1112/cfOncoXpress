# Weiling Li (wel4007@med.cornell.edu)
# Marjorie Roskes (mcl4001@med.cornell.edu)
# power calculation

import numpy as np
import pandas as pd
import scipy.stats
from sympy.solvers import solve
import sympy as sp  # otherwise, confuse with sp.sqrt, and np.sqrt

import sys

fixedsize = 34

############## AR, or NEPC signatures #############, 
TCGA_ANDROGEN_RESPONSE = pd.read_csv('../csv_files/TCGA_ANDROGEN_RESPONSE',header=None)[0].values
Beltran_NEPC_up = pd.read_csv('../csv_files/Beltran_NEPC_up', header=None)[0].values

mat = pd.read_csv('../GSEAdata/input/controls34_herberts63.log2TPM.minmax.final.txt', sep='\t')
mat.drop('DESCRIPTION', axis=1, inplace=True)
mat=mat.set_index(['NAME'])
ARpathwaymat = mat.loc[TCGA_ANDROGEN_RESPONSE]
#NEpathwaymat = mat.loc[Beltran_NEPC_up]
valid_labels = [label for label in Beltran_NEPC_up if label in mat.index]
NEpathwaymat = mat.loc[valid_labels]

srr_columns = [col for col in mat.columns if col.startswith('SRR')]
ne_columns = ['GU-16-026-cfDNA', 'GU-17-295-cfDNA']
ar_columns = [col for col in mat.columns if not col.startswith('SRR')]
ar_columns.remove('GU-16-026-cfDNA')
ar_columns.remove('GU-17-295-cfDNA')

############
sAR=[] 
for samplename in  ar_columns:
    sAR.extend(ARpathwaymat[samplename].values)
mu_ar = np.mean(sAR)         ######################### lwl
sigma_ar = np.std(sAR)       ######################### lwl
#print(mu_ar, sigma_ar)

sARinPBMC =[] # tune
for samplename in  srr_columns:
    sARinPBMC.extend(ARpathwaymat[samplename].values) ######## lwl: AR pathway.
mu_pbmcar = np.mean(sARinPBMC) # tuning....      ######################### lwl
sigma_pbmcar = np.std(sARinPBMC)                 ######################### lwl
print(mu_ar/mu_pbmcar)
###############
sNE=[]
for samplename in  ne_columns:
    sNE.extend(NEpathwaymat[samplename].values)
mu_ne = np.mean(sNE)         ######################### lwl
sigma_ne = np.std(sNE)       ######################### lwl

sNEinPBMC =[] # tune
for samplename in  srr_columns:
    sNEinPBMC.extend(NEpathwaymat[samplename].values) ######## lwl: NE pathway.
mu_pbmcne = np.mean(sNEinPBMC) # tuning....      ######################### lwl
sigma_pbmcne = np.std(sNEinPBMC)                 ######################### lwl
print(mu_ne/mu_pbmcne)
###############

def power_calculate(mu_ar_proj, mu_ne_proj, sigma_ar, sigma_ne, n2target):
    cm_proj = (mu_ar_proj + mu_ne_proj) / 2
    prob_ar_ne = scipy.stats.norm.cdf(cm_proj, loc=mu_ar_proj, scale=sigma_ar) 
    prob_ne_ne = scipy.stats.norm.cdf(cm_proj, loc=mu_ne_proj, scale=sigma_ne)
    prob_ar_ar = 1-prob_ar_ne
    prob_ne_ar = 1-prob_ne_ne
    
    #p1, p2 = prob_ar_ar, prob_ar_ne # p1 = P(AR|s=AR), p2 = P(AR|s=NE)
    p1, p2  = prob_ar_ar, prob_ne_ar
    q1, q2 = 1 - p1, 1 - p2
    zalpha, zbeta = 1.96, 0.84
    # pbar = (p1 + (K * p2)) / (1 + K)
    # qbar = 1 - pbar = 1 - (p1 + (K * p2)) / (1 + K)

    # n1target = (((zalpha * np.sqrt(pbar * qbar * (1 + (1 / K)))) + (zbeta * np.sqrt((p1 * q1) + ((p2 * q2) / K)))) / (abs(p2 - p1))) ** 2
#     n2target = K * n1target = K *....: NE group size is constant
    K = sp.Symbol('K')
    sol = solve( K* ((((zalpha * sp.sqrt((p1 + (K * p2)) / (1 + K)  * (1 - (p1 + (K * p2)) / (1 + K)) * (1 + (1 / K)))) + (zbeta * sp.sqrt((p1 * q1) + ((p2 * q2) / K)))) / (abs(p2 - p1))) ** 2)
          - n2target , K)
    print(sol)
    return sol#complex(sol[0])


R = 1
n2target = fixedsize # 
#n1target, fc = [],[]

######## AR vs pbmc ##########
mu_ar_projs = (np.linspace(1, 5, 401)) * mu_pbmcar #mu_ne # tune mu_ar_proj = 1 to 5 (FC) * mu_ne 
# mu_ne_proj = mu_ne
mu_pbmc_ar_proj = mu_pbmcar #mu_ne

Ks=[]
for mu_ar_proj in mu_ar_projs: 
    #print(mu_ar_proj/mu_pbmc_ar_proj)
    K = power_calculate(mu_ar_proj, mu_pbmc_ar_proj, sigma_ar*R, sigma_pbmcar*R, n2target)    
    Ks.append(K)
    # n1target.append(n2target / K)
    # fc.append(mu_ar_proj / mu_ne_proj)

outdf=pd.DataFrame(columns=['K','ARsize','FoldChange'])
for i in range(0,len(Ks)):
    if len(Ks[i])>0: 
        newrow ={'K': Ks[i][0], 'ARsize': n2target / complex(Ks[i][0]), 'FoldChange': mu_ar_projs[i] / mu_pbmc_ar_proj}
        new_df = pd.DataFrame([newrow])
        outdf = pd.concat([outdf, new_df], axis=0, ignore_index=True)
outdf.to_csv(f'../csv_files/results/fc_1_5.AR.powercalcfinal.txt',header=True,index=False)

######## NE vs pbmc ##########
# mu_ne_projs = (np.linspace(1, 5, 401)) * mu_pbmcne #mu_ne # tune mu_ar_proj = 1 to 5 (FC) * mu_ne 
# mu_pbmc_ne_proj = mu_pbmcne #mu_ne

# Ks2=[]
# for mu_ne_proj in mu_ne_projs: 
#     K = power_calculate(mu_ne_proj, mu_pbmc_ne_proj, sigma_ne*R, sigma_pbmcne*R, n2target)    
#     Ks2.append(K)
#     # n1target.append(n2target / K)
#     # fc.append(mu_ar_proj / mu_ne_proj)

# outdf=pd.DataFrame(columns=['K','NEsize','FoldChange'])
# for i in range(0,len(Ks2)):
#     if len(Ks2[i])>0: 
#         newrow ={'K': Ks2[i][0], 'NEsize': n2target / complex(Ks2[i][0]), 'FoldChange': mu_ne_projs[i] / mu_pbmc_ne_proj}
#         new_df = pd.DataFrame([newrow])
#         outdf = pd.concat([outdf, new_df], axis=0, ignore_index=True)
# outdf.to_csv(f'ismb/fc_1_5.NE',header=True,index=False)





