# Marjorie Roskes (mcl4001@med.cornell.edu)

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
    
    @property
    def delfi_dir(self): 
        directory = self.analysis_dir/'fragment_lengths/delfi'
        directory.mkdir(parents=True, exist_ok=True)
        return directory
    
    def short_memmap_path(self, chrom): return self.delfi_dir/f'short_fragment_counts_{chrom}.memmap'
    def short_lock_path(self, chrom): return self.delfi_dir/f'short_fragment_counts_{chrom}.lock'
    def short_complete(self, chrom): return (self.short_memmap_path(chrom=chrom).exists()) and (not self.short_lock_path(chrom=chrom).exists())
    
    def long_memmap_path(self, chrom): return self.delfi_dir/f'long_fragment_counts_{chrom}.memmap'
    def long_lock_path(self, chrom): return self.delfi_dir/f'long_fragment_counts_{chrom}.lock'
    def long_complete(self, chrom): return (self.long_memmap_path(chrom=chrom).exists()) and (not self.long_lock_path(chrom=chrom).exists())
    
    
    
def load_samples(): 
    return [Sample(sample_name='DTB-097-Progression-cfDNA'),
           ]


def short_fragment_counts(sample, chrom, size): 
    if sample.short_complete(chrom=chrom): 
        logger.info(f'{sample.sample_name}: short_fragment_counts {chrom} already finished runninng successfully')
        return
    
    with job_lock.JobLock(sample.short_lock_path(chrom=chrom), outputfiles=[sample.short_memmap_path(chrom=chrom)]) as lock: 
        if not lock: return
        logger.info(f'{sample.sample_name}: short_fragment_counts {chrom} running')
        bam_path = job_lock.slurm_rsync_input(sample.bam_path(chrom=chrom))
        bai_path = job_lock.slurm_rsync_input(sample.bai_path(chrom=chrom))
        with job_lock.slurm_rsync_output(sample.short_memmap_path(chrom=chrom)) as short_memmap_path:
            bam_tools.get_coverage_memmap(bam_path=bam_path, 
                                          chrom=chrom, 
                                          size=size, 
                                          out_path=short_memmap_path, 
                                          read_quality_threshold=20, 
                                          only_save_to_midpoint=True, 
                                          min_fragment_length=100, 
                                          max_fragment_length=150
                                         )
        return
    
def long_fragment_counts(sample, chrom, size): 
    if sample.long_complete(chrom=chrom): 
        logger.info(f'{sample.sample_name}: long_fragment_counts {chrom} already finished runninng successfully')
        return
    
    with job_lock.JobLock(sample.long_lock_path(chrom=chrom), outputfiles=[sample.long_memmap_path(chrom=chrom)]) as lock: 
        if not lock: return
        logger.info(f'{sample.sample_name}: long_fragment_counts {chrom} running')
        bam_path = job_lock.slurm_rsync_input(sample.bam_path(chrom=chrom))
        bai_path = job_lock.slurm_rsync_input(sample.bai_path(chrom=chrom))
        with job_lock.slurm_rsync_output(sample.long_memmap_path(chrom=chrom)) as long_memmap_path:
            bam_tools.get_coverage_memmap(bam_path=bam_path, 
                                          chrom=chrom, 
                                          size=size, 
                                          out_path=long_memmap_path, 
                                          read_quality_threshold=20, 
                                          only_save_to_midpoint=True, 
                                          min_fragment_length=151, 
                                          max_fragment_length=220
                                         )
        return

def main(): 
    samples = load_samples()
    chrom_sizes = chrom_tools.chrom_sizes(ref_build='GRCh38')
    for i, (chrom, size) in enumerate(chrom_sizes.items()):
        for j, sample in enumerate(samples): 
            short_fragment_counts(sample=sample, chrom=chrom, size=size)
            long_fragment_counts(sample=sample, chrom=chrom, size=size)
        

if __name__=='__main__':
    main()