from modules import bam_tools, chrom_tools, misc_tools, logger_tools

import argparse
import collections
import datetime
import methodtools
import numpy as np
import pandas as pd
import pathlib
import pysam
import re
import sys
import time

class Sample():
    def __init__(self, patient_id, sample_id):#, status, ar_status, tumor_fraction, plasma_source, total_coverage,sex):
        #         patient_id	sample_id
# AE-015	AE-015-Baseline-cfDNA
        self.patient_id = patient_id
        self.sample_id = sample_id
        # self.status = status
        # self.ar_status = ar_status
        # self.tumor_fraction = tumor_fraction
        # self.total_coverage = total_coverage
        # self.plasma_source = plasma_source
        # self.sex = sex
        self.ref_build = 'GRCh38'
        self.chrom_sizes = chrom_tools.chrom_sizes(ref_build=self.ref_build)
        self.chromosomes = [f'chr{i+1}' for i in range(22)] + ['chrX', 'chrY']
        # self.color = list(np.random.choice(range(256), size=3)/255)
        # self.histological_status_color = misc_tools.histological_status_color(status=self.status, 
        #                                                                       ar_status=self.ar_status, 
        #                                                                       plasma_source=self.plasma_source,
        #                                                                       sex=self.sex)

    @property
    def analysis_dir(self):
        #return pathlib.Path(f'scratch/mcl4001/cfDNA/analysis/Sample_{self.name}')   
        return pathlib.Path(f'scratch/mcl4001/herberts_data/analysis/{self.patient_id}/cfDNA')   
################    for entropy
    @property
    def data_dir(self):
        return pathlib.Path(f'scratch/mcl4001/herberts_data/data/{self.patient_id}/cfDNA')
#         return pathlib.Path(f'scratch/mcl4001/cfDNA/data/Sample_{self.name}')    
        

        
    @property
    def bam_per_chrom_dir(self): return self.data_dir/'bam'/'bam_per_chrom'
        
    def bam_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_id}_{chrom}.bam'
    def bai_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_id}_{chrom}.bai'
 ################   
    #scratch/mcl4001/cfDNA/analysis/Sample_*/wps/nps/fft_freq190/950
    @methodtools.lru_cache()
    def nps(self, chrom, thresh=None, rsync=False):
        path = self.analysis_dir/f'{self.sample_id}_wps'/'wps'/'nps'/f'wps_{chrom}.nps.memmap'
        if rsync is True: 
            path = job_lock.slurm_rsync_input(path)
        memmap = np.memmap(path, mode='r', dtype=np.float32)
        assert len(memmap) == self.chrom_sizes[chrom]
        if thresh is not None: memmap = memmap > thresh
        return memmap
    

###################################    
    def coverage(self, chrom, thresh=None, rsync=False): #the median normalization of wps.
        path = self.analysis_dir/'coverage'/f'{self.sample_id}_coverage_{chrom}.memmap'
        if rsync is True: 
            path = job_lock.slurm_rsync_input(path)
        memmap = np.memmap(path, mode='r', dtype=np.int32)
        assert len(memmap) == self.chrom_sizes[chrom]
        if thresh is not None: memmap = memmap > thresh
        return memmap

    
    @methodtools.lru_cache()
    def mlc(self, chrom):
        path = self.analysis_dir/'coverage'/f'mlc_{1000000}'/f'{self.sample_id}_mlc_1000000_{chrom}.memmap'
        memmap = np.memmap(path,mode='r',dtype=np.float32)
        return memmap
    
    @methodtools.lru_cache()
    def bp_short_counts(self, chrom): 
        path = self.analysis_dir/'fragment_lengths'/'delfi'/f'{self.sample_id}_short_fragment_counts_{chrom}.memmap'
        memmap = np.memmap(path, mode='r', dtype=np.int32)
        return memmap
    
    @methodtools.lru_cache()
    def bp_long_counts(self, chrom): 
        path = self.analysis_dir/'fragment_lengths'/'delfi'/f'{self.sample_id}_long_fragment_counts_{chrom}.memmap'
        memmap = np.memmap(path, mode='r', dtype=np.int32)
        return memmap
    
    @methodtools.lru_cache()
    def get_fft_in_region(self, chrom, start, end, mlc_normalize=True, region_size_normalize=True, abs_it=True):
        region_size = end - start
        assert region_size % 190 == 0
        nps = self.nps(chrom=chrom)[start:end]
        if mlc_normalize: 
            mlc = self.mlc(chrom=chrom)[start:end]
            nps = nps/ np.where(mlc!=0, mlc, 1)
        
        freqs = np.fft.fftfreq(region_size, d=1)
        posfreqidx = freqs>0
        posfreqs = freqs[posfreqidx]
        poi = 190
        foi = 1/poi
        ioi = np.argmin(abs(posfreqs-foi))
        
        fft = np.fft.fft(nps)[posfreqidx][ioi]
        
        if region_size_normalize: fft = fft * (2 / region_size)
        
        if abs_it: fft = np.abs(fft)

        return fft
    
    
    @methodtools.lru_cache()
    def get_fp_in_region(self, chrom, start, end): 
        short = self.bp_short_counts(chrom=chrom)[start:end]
        long = self.bp_long_counts(chrom=chrom)[start:end]
        
        numshort = np.sum(short)
        numlong = np.sum(long)
        
        #L/(L+s)
        #print(numshort, numlong)
        fp = (numshort / numlong) if numlong!=0 else np.nan
        
        return fp
    
    @methodtools.lru_cache()
    def get_coverage_in_region(self, chrom, start, end, mlc_normalize=True): 
        coverage = self.coverage(chrom=chrom)[start:end]
        if mlc_normalize: 
            mlc = self.mlc(chrom=chrom)[start:end]
            coverage = coverage/ np.where(mlc!=0, mlc, 1)
 
        return coverage




def load_samples():
    # df = pd.read_csv('scratch/mcl4001/herberts_data/data/metadata.txt', sep='\t', header=0, usecols=[0,1,2])
    # patient_id, sample_id, specimens = np.asarray(df['patient_id']), np.asarray(df['sample_id']), np.asarray(df['specimen'])
    samples = {}
    # for i, (patient_id, sample_id, specimen) in enumerate(zip(patient_id, sample_id, specimens)):
    #     if specimen != 'cfDNA': continue
    #     samples[sample_id] = Sample(patient_id=patient_id, sample_id=sample_id)
    df = pd.read_csv('../csv_files/herberts_meta.txt', sep='\s+',header=None)
    patient_id, sample_id = np.asarray(df[0]), np.asarray(df[1])
    for i, (patient_id, sample_id) in enumerate(zip(patient_id, sample_id)):
        samples[sample_id] = Sample(patient_id=patient_id, sample_id=sample_id)
    return samples
        

