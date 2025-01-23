import sys
from modules import bam_tools, chrom_tools, misc_tools, logger_tools

import job_lock

import argparse
import collections
import datetime
import methodtools
import numpy as np
import pandas as pd
import pathlib
import pysam
import re
import time


logger = logger_tools.logger

job_lock.JobLock.defaultcorruptfiletimeout = datetime.timedelta(hours=1)
#job_lock.JobLock.setdefaultminimumtimeforiterativelocks(datetime.timedelta(minutes=10))

class Sample():
    def __init__(self, patient_id, sample_id):#name, status, ar_status, tumor_fraction, plasma_source, total_coverage,sex):
        self.patient_id = patient_id
        self.sample_id = sample_id
        # self.name = name
        # self.status = status
        # self.ar_status = ar_status
        # self.tumor_fraction = tumor_fraction
        # self.total_coverage = total_coverage
        # self.plasma_source = plasma_source
        # self.sex = sex
        # self.ref_build = 'GRCh38'
        # self.chrom_sizes = chrom_tools.chrom_sizes(ref_build=self.ref_build)
        # self.chromosomes = [f'chr{i+1}' for i in range(22)] + ['chrX', 'chrY']
        # self.color = list(np.random.choice(range(256), size=3)/255)
        # self.histological_status_color = misc_tools.histological_status_color(status=self.status, 
        #                                                                       ar_status=self.ar_status, 
        #                                                                       plasma_source=self.plasma_source,
        #                                                                       sex=self.sex)

    @property
    def data_dir(self):
        return pathlib.Path(f'scratch/mcl4001/herberts_data/data/{self.patient_id}/cfDNA')
    
#     @methodtools.lru_cache()
#     def bam_path(self, chrom):
#         bamdir = self.data_dir/'bam'/'bam_per_chrom'
#         file_names = list(bamdir.glob(f'{self.name}_{chrom}.bam')) # lwl
# #         file_names = list(bamdir.glob(f'*.REF_{chrom.replace("chr", "*")}.bam')) #lwl
# #         file_names = [file_name for file_name in file_names if re.match(f'.*REF_[a-z]{{0,3}}{chrom.replace("chr", "")}.bam$', file_name.name) is not None] #lwl
#         if len(file_names) != 1: 
#             logger.info(f'{self.name}, {chrom}, {self.data_dir}.')
#             logger.info(f'filenames: {file_names}.')        
#         assert len(file_names) == 1
#         return file_names[0]
    
#     @methodtools.lru_cache()
#     def bam(self, chrom):
#         bam_path = self.bam_path(chrom=chrom)
# #         bai_path = bam_path.with_suffix('.bam.bai')
#         bai_path = bam_path.with_suffix('.bai') # lwl
# #         bam_path = job_lock.slurm_rsync_input(bam_path) # recover
# #         bai_path = job_lock.slurm_rsync_input(bai_path) # recover
#         return pysam.AlignmentFile(bam_path)
    @property
    def bam_per_chrom_dir(self): return self.data_dir/'bam'/'bam_per_chrom' #old bam_path
    
    @methodtools.lru_cache()
    def bam(self, chrom):
        bam_path = self.bam_per_chrom_dir/f'{self.sample_id}_{chrom}.bam'
        bai_path = self.bam_per_chrom_dir/f'{self.sample_id}_{chrom}.bai'

        bam_path = job_lock.slurm_rsync_input(bam_path)
        bai_path = job_lock.slurm_rsync_input(bai_path)
        return pysam.AlignmentFile(bam_path)

class Bins():
    def __init__(self, path):
        self.path = pathlib.Path(path)
        
    @methodtools.lru_cache()
    @property
    def bins(self): 
        df = pd.read_csv(self.path, sep='\t', header=0)
        columns = df.columns.tolist()
        if (columns[0]!='chr') or (columns[1]!='start') or (columns[2]!='end'):
            raise AssertionError(f'Bed file does not have appropriate header! First three column names should be (1) chr (2) start and (3) end.')
        return df
    
    @methodtools.lru_cache()
    @property
    def num_bins(self):
        return self.bins.shape[0]
        
    
    
# def load_samples(metadata_path): 
#     metadata = pd.read_csv(metadata_path, sep='\t', header=0)    
#     samples = {}
#     for idx, row in metadata.iterrows():
#         name, status, ar_status, tumor_fraction, total_coverage, sex, plasma_source = row['id'], row['histological_status'], row['ar_status'], row['ichorCNA_tumor_fraction'], row['total_coverage'], row['sex'], row['plasma_source']
#         #if 'control' in name: continue
#         samples[name] = Sample(name=name, status=status,ar_status=ar_status, tumor_fraction=tumor_fraction, total_coverage=total_coverage, sex=sex, plasma_source=plasma_source)    
#     return samples 

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

def get_se_in_bins(sample, bins, out_dir=None, min_fragment_length=0, max_fragment_length=500, quality_threshold=20): 
    if out_dir is None:
        out_dir = bins.path.parents[0]
    else: 
        out_dir = pathlib.Path(out_dir)
        
    out_path = out_dir/f'{bins.path.stem}.simpson_probability.{sample.sample_id}.{min_fragment_length}.{max_fragment_length}.lwl.bed'
    lock_path = out_path.with_suffix('.lock')
        
    if (out_path.exists()) and (not lock_path.exists()):
        logger.info(f'{sample.sample_id}: Already finished running successfully')
        return
    with job_lock.JobLock(lock_path, outputfiles=[out_path]) as lock: 
        if not lock: 
            logger.info(f'{sample.sample_id}: Another process is already running this task')
            return
        logger.info(f'{sample.sample_id}: Getting SE in bins {bins.path.name}')
        
        chroms, starts, ends = np.asarray(bins.bins['chr']), np.asarray(bins.bins['start']), np.asarray(bins.bins['end'])
        
        num_intervals = 4
        simpson_probability = pd.DataFrame(columns=['i1','i2','i3','i4'])
        if min_fragment_length >70:
            num_intervals = 3
            simpson_probability = pd.DataFrame(columns=['i1','i2','i3'])
        #simpson_probability = #np.zeros(bins.num_bins, dtype=np.float32)
        
        for i in range(bins.num_bins): 
            if i%1000==0: logger.info(f'{sample.sample_id}: {i}/{bins.num_bins}')

            chrom, start, end = chroms[i], starts[i], ends[i]

            regionstart, regionend = start-max_fragment_length-1, end+max_fragment_length+1
            
            bin_n_i = np.zeros(num_intervals, dtype=np.int32)

            region = f'{chrom}:{regionstart}-{regionend}'
            bam = sample.bam(chrom=chrom)
            for j, (read1, read2) in enumerate(bam_tools.read_pair_generator(bam=bam, region_string=region)):
                pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
                min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
                if (min_ind>end) or (max_ind<start): continue
                frag_len = max_ind - min_ind
                if (frag_len<min_fragment_length) or (frag_len>max_fragment_length): continue

                if (quality_threshold is not None) and ((read1.mapping_quality + read2.mapping_quality) < 2*quality_threshold): continue
                
#                 if (frag_len>=min_fragment_length) and (frag_len<70): bin_n_i[0] += 1
#                 elif (frag_len>=70) and (frag_len<170): bin_n_i[1] += 1
#                 elif (frag_len>=170) and (frag_len<210): bin_n_i[2] += 1
#                 elif (frag_len>=210) and (frag_len<max_fragment_length): bin_n_i[3] += 1
                if min_fragment_length <70:
                    if (frag_len>=min_fragment_length) and (frag_len<70): bin_n_i[0] += 1
                    elif (frag_len>=70) and (frag_len<170): bin_n_i[1] += 1
                    elif (frag_len>=170) and (frag_len<210): bin_n_i[2] += 1
                    elif (frag_len>=210) and (frag_len<max_fragment_length): bin_n_i[3] += 1
                else:
                    if (frag_len>=min_fragment_length) and (frag_len<170): bin_n_i[0] += 1
                    elif (frag_len>=170) and (frag_len<210): bin_n_i[1] += 1
                    elif (frag_len>=210) and (frag_len<max_fragment_length): bin_n_i[2] += 1

            # for alex's simpson entropy
            # total number of fragments (of relevant length) in the bin
            bin_n = np.sum(bin_n_i)
            # ratio of number of fragments of len in interval i to total number of fragments in bin
            # nan in no fragments in the bin
            p_hat_i = bin_n_i / bin_n            
            # compute the simpson entropy in the bin
            #simpson = 1/np.sum((p_hat_i**2))
            simpson_probability.loc[i] = p_hat_i
            
        
#         df = pd.DataFrame({'chr':chroms, 
#                            'start':starts, 
#                            'end':ends, 
#                            'simpson_probability':simpson_probability,
#                           })

        #with job_lock.slurm_rsync_output(out_path) as out_path:
        with out_path as out_path:
            simpson_probability.to_csv(out_path, sep='\t', header=False, index=False)
    return
    
def main():
    # init parser.
    p = argparse.ArgumentParser()
    p.add_argument("--sample-id", type=str, required=True)
    p.add_argument("--bed-path", type=str, required=True)
    p.add_argument("--out-dir", type=str, required=False, default=None)
    p.add_argument("--min-fragment-length", type=int, required=False, default=0)
    p.add_argument("--max-fragment-length", type=int, required=False, default=500)
    p.add_argument("--quality-threshold", type=int, required=False, default=20)
    

    args = p.parse_args()    
    
    samples = load_samples()#metadata_path='scratch/mcl4001/cfDNA/data/metadata.txt')
    if args.sample_id not in samples: 
        raise AssertionError(f'Invalid sample name provided! Please choose from: {list(samples.keys())}')
        
    sample = samples[args.sample_id]
    bins = Bins(path=args.bed_path)
    get_se_in_bins(sample=sample, 
                   bins=bins, 
                   out_dir=args.out_dir, 
                   min_fragment_length=args.min_fragment_length, 
                   max_fragment_length=args.max_fragment_length, 
                   quality_threshold=args.quality_threshold)
    
    
    return

if __name__=='__main__':
    main()