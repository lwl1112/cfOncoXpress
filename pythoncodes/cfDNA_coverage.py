from modules import logger_tools, bam_tools, misc_tools, chrom_tools

import job_lock

import datetime
import methodtools
import numpy as np
import os
import pandas as pd
import pathlib
import shutil
import subprocess
import tempfile
import time

logger = logger_tools.logger

job_lock.JobLock.defaultcorruptfiletimeout = datetime.timedelta(hours=1)
job_lock.JobLock.setdefaultminimumtimeforiterativelocks(datetime.timedelta(minutes=10))

class Sample():
    def __init__(self, sample_name): 
        self.sample_name = sample_name
    
    @property
    def data_dir(self): return pathlib.Path(f'scratch/mcl4001/cfDNA/data/Sample_{self.sample_name}')
    
    @property
    def analysis_dir(self): 
        directory = pathlib.Path(f'scratch/mcl4001/cfDNA/analysis/Sample_{self.sample_name}')
        directory.mkdir(parents=True, exist_ok=True)
        return directory
        
    @property
    def bam_per_chrom_dir(self): return self.data_dir/'bam'/'bam_per_chrom'
        
    def bam_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_name}_{chrom}.bam'
    def bai_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_name}_{chrom}.bai'
    
    def index_complete(self,chrom): return (self.bai_path(chrom=chrom).exists()) and (not self.bai_path(chrom=chrom).with_suffix('.index.lock').exists())
    
    @property
    def coverage_dir(self): 
        directory = self.analysis_dir/'coverage'
        directory.mkdir(parents=True, exist_ok=True)
        return directory
    
    def coverage_memmap_path(self, chrom): return self.coverage_dir/f'coverage_{chrom}.memmap'
    def coverage_lock_path(self, chrom): return self.coverage_dir/f'coverage_{chrom}.lock'
    def coverage_complete(self, chrom): return (self.coverage_memmap_path(chrom=chrom).exists()) and (not self.coverage_lock_path(chrom=chrom).exists())
    
    @methodtools.lru_cache()
    def mlc_dir(self, window=1000000): 
        directory = self.analysis_dir/'coverage'/f'mlc_{window}'
        directory.mkdir(parents=True, exist_ok=True)
        return directory
    
    def mlc_memmap_path(self, chrom, window=1000000): return self.mlc_dir(window=window)/f'mlc_{chrom}.memmap'
    def mlc_lock_path(self, chrom, window=1000000): return self.mlc_dir(window=window)/f'mlc_{chrom}.lock'
    def mlc_complete(self, chrom, window=1000000): return (self.mlc_memmap_path(chrom=chrom, window=window).exists()) and (not self.mlc_lock_path(chrom=chrom, window=window).exists())

    @property
    def depth_dir(self): 
        directory = self.analysis_dir/'coverage'/'depth'
        directory.mkdir(parents=True, exist_ok=True)
        return directory
    
    def depth_path(self, chrom): return self.depth_dir/f'total_depth_{chrom}.tsv'    
    def depth_lock_path(self, chrom): return self.depth_dir/f'total_depth_{chrom}.lock'
    def depth_complete(self, chrom): return (self.depth_path(chrom=chrom).exists()) and (not self.depth_lock_path(chrom=chrom).exists())

    @property
    def read_quality_dir(self): 
        directory = self.analysis_dir/'read_quality'
        directory.mkdir(parents=True, exist_ok=True)
        return directory
    
    def read_quality_count_path(self, chrom): return self.read_quality_dir/f'read_quality_{chrom}.tsv'
    def read_quality_count_lock_path(self, chrom): return self.read_quality_dir/f'read_quality_{chrom}.lock'
    def read_quality_count_complete(self, chrom): return (self.read_quality_count_path(chrom=chrom).exists()) and (not self.read_quality_count_lock_path(chrom=chrom).exists())
    
    @property
    def fragment_length_dir(self): 
        directory = self.analysis_dir/'fragment_lengths'
        directory.mkdir(parents=True, exist_ok=True)
        return directory
    
    def fragment_length_count_path(self, chrom): return self.fragment_length_dir/f'fragment_length_{chrom}.tsv'
    def fragment_length_count_lock_path(self, chrom): return self.fragment_length_dir/f'fragment_length_{chrom}.lock'
    def fragment_length_count_complete(self, chrom): return (self.fragment_length_count_path(chrom=chrom).exists()) and (not self.fragment_length_count_lock_path(chrom=chrom).exists())
        

        
#def load_samples(): 
#    return [Sample(sample_name='WCMFC01'), 
#           ]
#def load_samples():
#    metadata = pd.read_csv('scratch/mcl4001/cfDNA/data/metadata.txt', sep='\t', header=0)        
#    sample_names = list(metadata['id'])
#    return [Sample(sample_name=sample_name) for sample_name in sample_names if 'control' not in sample_name]

def load_samples(): 
    #return [Sample(sample_name='test')]
 
    return [Sample(sample_name='DTB-097-Progression-cfDNA')]

#def load_samples(): 
#    sample_names = ['PM1602_0.1X','PM3014_0.1X', 'PM3088_0.1X', 'PM3146_0.1X', 'PM3165_0.1X', 'PM3282_0.1X', #'PM794_0.1X','WCMFC01']
#    return [Sample(sample_name=sample_name) for sample_name in sample_names]


def coverage(sample, chrom, size): 
    if not sample.index_complete(chrom=chrom):
        logger.info(f'{sample.sample_name}: index not complete --> not ready for QC')
        return
    
    if sample.coverage_complete(chrom=chrom): 
        logger.info(f'{sample.sample_name}: coverage {chrom} already finished runninng successfully')
        return
    
    with job_lock.JobLock(sample.coverage_lock_path(chrom=chrom), outputfiles=[sample.coverage_memmap_path(chrom=chrom)]) as lock: 
        if not lock: return
        logger.info(f'{sample.sample_name}: coverage {chrom} running')
        # bam_path = job_lock.slurm_rsync_input(sample.bam_path(chrom=chrom))
        # bai_path = job_lock.slurm_rsync_input(sample.bai_path(chrom=chrom))
        bam_path = sample.bam_path(chrom=chrom)
        bai_path = sample.bai_path(chrom=chrom)
        # with job_lock.slurm_rsync_output(sample.coverage_memmap_path(chrom=chrom)) as coverage_memmap_path:
        coverage_memmap_path = sample.coverage_memmap_path(chrom=chrom)
        bam_tools.get_coverage_memmap(bam_path=bam_path, 
                                      chrom=chrom, 
                                      size=size, 
                                      out_path=coverage_memmap_path, 
                                      read_quality_threshold=20, 
                                      only_save_to_midpoint=False, 
                                      min_fragment_length=None, 
                                      max_fragment_length=None
                                     )
        return
        
    
def depth(sample, chrom, size): 
    if not sample.index_complete(chrom=chrom):
        logger.info(f'{sample.sample_name}: index not complete --> not ready for QC')
        return
    
    if sample.depth_complete(chrom=chrom): 
        logger.info(f'{sample.sample_name}: depth {chrom} already finished runninng successfully')
        return
    
    with job_lock.JobLock(sample.depth_lock_path(chrom=chrom), outputfiles=[sample.depth_path(chrom=chrom)]) as lock: 
        if not lock: return
        logger.info(f'{sample.sample_name}: depth {chrom} running')
        # bam_path = job_lock.slurm_rsync_input(sample.bam_path(chrom=chrom))
        # bai_path = job_lock.slurm_rsync_input(sample.bai_path(chrom=chrom))
        bam_path = sample.bam_path(chrom=chrom)
        bai_path = sample.bai_path(chrom=chrom)
        
        temp_unfiltered_path = pathlib.Path(os.environ['TMPDIR'])/f'{sample.sample_name}.depth_{chrom}.unfiltered.bed'
        bam_tools.get_depth(bam_path=bam_path,
                            out_path=temp_unfiltered_path, 
                            chrom=chrom, 
                            filtered=False,
                           )
        unfiltered_df = pd.read_csv(temp_unfiltered_path, sep='\t', header=0, names=['chr', 'pos', 'depth'])
        print(unfiltered_df)
        
        temp_filtered_path = pathlib.Path(os.environ['TMPDIR'])/f'{sample.sample_name}.depth_{chrom}.filtered.bed'
        bam_tools.get_depth(bam_path=bam_path,
                            out_path=temp_filtered_path, 
                            chrom=chrom, 
                            filtered=True,
                           )        
        filtered_df = pd.read_csv(temp_filtered_path, sep='\t', header=0, names=['chr', 'pos', 'depth'])
        print(filtered_df)
        
        df = pd.DataFrame({'chr': [chrom], 
                           'size': [size],
                           'unfiltered_total_depth':[np.sum(unfiltered_df['depth'])],
                           'filtered_total_depth':[np.sum(filtered_df['depth'])],
                          })
        # with job_lock.slurm_rsync_output(sample.depth_path(chrom=chrom)) as depth_path:
        depth_path = sample.depth_path(chrom=chrom)
        df.to_csv(depth_path, sep='\t', header=True, index=False)
        return

def read_quality_count(sample, chrom): 
    if not sample.index_complete(chrom=chrom):
        logger.info(f'{sample.sample_name}: index not complete --> not ready for QC')
        return
    
    if sample.read_quality_count_complete(chrom=chrom): 
        logger.info(f'{sample.sample_name}: read_quality_count {chrom} already finished runninng successfully')
        return
    
    with job_lock.JobLock(sample.read_quality_count_lock_path(chrom=chrom), outputfiles=[sample.read_quality_count_path(chrom=chrom)]) as lock: 
        if not lock: return
        logger.info(f'{sample.sample_name}: read_quality_count {chrom} running')
        # bam_path = job_lock.slurm_rsync_input(sample.bam_path(chrom=chrom))
        # bai_path = job_lock.slurm_rsync_input(sample.bai_path(chrom=chrom))
        bam_path = sample.bam_path(chrom=chrom)
        bai_path = sample.bai_path(chrom=chrom)
        
        # with job_lock.slurm_rsync_output(sample.read_quality_count_path(chrom=chrom)) as read_quality_count_path: 
        read_quality_count_path = sample.read_quality_count_path(chrom=chrom)
        df = bam_tools.get_read_mapping_quality_counts(bam_path=bam_path, 
                                                       out_path=read_quality_count_path,
                                                       region=chrom,
                                                      )
        return

def fragment_length_count(sample, chrom): 
    if not sample.index_complete(chrom=chrom):
        logger.info(f'{sample.sample_name}: index not complete --> not ready for QC')
        return
    
    if sample.fragment_length_count_complete(chrom=chrom): 
        logger.info(f'{sample.sample_name}: fragment_length_count {chrom} already finished runninng successfully')
        return
    
    with job_lock.JobLock(sample.fragment_length_count_lock_path(chrom=chrom), outputfiles=[sample.fragment_length_count_path(chrom=chrom)]) as lock: 
        if not lock: return
        logger.info(f'{sample.sample_name}: fragment_length_count {chrom} running')
        # bam_path = job_lock.slurm_rsync_input(sample.bam_path(chrom=chrom))
        # bai_path = job_lock.slurm_rsync_input(sample.bai_path(chrom=chrom))
        bam_path = sample.bam_path(chrom=chrom)
        bai_path = sample.bai_path(chrom=chrom)
        
        # with job_lock.slurm_rsync_output(sample.fragment_length_count_path(chrom=chrom)) as fragment_length_count_path: 
        fragment_length_count_path = sample.fragment_length_count_path(chrom=chrom)
        df = bam_tools.get_fragment_length_counts(bam_path=bam_path, 
                                                  out_path=fragment_length_count_path,
                                                  region=chrom,
                                                 )            
        return
    
def mlc(sample, chrom, size, window=1000000): 
    #if not sample.index_complete(chrom=chrom):
    #    logger.info(f'{sample.sample_name}: index not complete --> not ready for QC')
    #    return

    if not sample.coverage_complete(chrom=chrom): 
        logger.info(f'{sample.sample_name}: coverage {chrom} has not yet finished runninng successfully')
        return


    if sample.mlc_complete(chrom=chrom, window=window): 
        logger.info(f'{sample.sample_name}: mlc {chrom} already finished runninng successfully')
        return
    
    with job_lock.JobLock(sample.mlc_lock_path(chrom=chrom, window=window), outputfiles=[sample.mlc_memmap_path(chrom=chrom, window=window)]) as lock: 
        if not lock: return
        logger.info(f'{sample.sample_name}: mlc {window} {chrom} running')
        
        # cov_path = job_lock.slurm_rsync_input(sample.coverage_memmap_path(chrom=chrom))
        cov_path = sample.coverage_memmap_path(chrom=chrom)
        cov_memmap = np.memmap(cov_path,mode='r', shape=(size,), dtype=np.int32)
        
        # with job_lock.slurm_rsync_output(sample.mlc_memmap_path(chrom=chrom, window=window)) as mlc_memmap_path: 
        mlc_memmap_path = sample.mlc_memmap_path(chrom=chrom, window=window)
        mlc_memmap = np.memmap(mlc_memmap_path, mode='w+', shape=(size,), dtype=np.float32)
        mlc_memmap[:] = 0
        mlc_memmap[:] = misc_tools.get_local_mean(array=cov_memmap, window=window)
        
        mlc_memmap.flush()
        return




def main(): 
    samples = load_samples()
    #samples = [Sample(sra_run_id='SRR23629305', pdx_series='test', sample_name='test', sample_type='cfDNA')]
    chrom_sizes = chrom_tools.chrom_sizes(ref_build='GRCh38')
    for i, (chrom, size) in enumerate(chrom_sizes.items()):
        for j, sample in enumerate(samples): 
            #if sample.sample_name !='': continue
            #if sample.sample_name not in ['']: continue
            coverage(sample=sample, chrom=chrom, size=size)
            mlc(sample=sample, chrom=chrom, size=size, window=1000000)
            #mlc(sample=sample, chrom=chrom, size=size, window=120)
            #depth(sample=sample, chrom=chrom, size=size)
            read_quality_count(sample=sample, chrom=chrom)
            fragment_length_count(sample=sample, chrom=chrom)

        

if __name__=='__main__':
    main()