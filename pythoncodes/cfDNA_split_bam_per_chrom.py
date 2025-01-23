from modules import logger_tools, bam_tools

import job_lock

import collections
import datetime
import numpy as np
import pandas as pd
import pathlib
import subprocess
import time


logger = logger_tools.logger

job_lock.JobLock.defaultcorruptfiletimeout = datetime.timedelta(hours=1)
job_lock.JobLock.setdefaultminimumtimeforiterativelocks(datetime.timedelta(minutes=10))
job_lock.JobLock.copyspeedlowerlimitbytespersecond = 1e6


class Sample():
    def __init__(self, sample_name): 
        self.sample_name = sample_name
    
    @property
    def data_dir(self): return pathlib.Path(f'scratch/mcl4001/cfDNA/data/Sample_{self.sample_name}/')
        
    @property
    def bam_dir(self): return self.data_dir/'bam'    
        
    @property
    def bam_path(self): return self.bam_dir/f'{self.sample_name}.final.bam'
    
    @property
    def bai_path(self): return self.bam_dir/f'{self.sample_name}.final.bai'
    
    @property
    def bam_per_chrom_dir(self):
        directory = self.bam_dir/f'bam_per_chrom'
        directory.mkdir(parents=True, exist_ok=True)
        return directory
        
    def split_bam_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_name}_{chrom}.bam'
    def split_lock_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_name}_{chrom}.split.lock'
    def split_complete(self, chrom): return (self.split_bam_path(chrom=chrom).exists()) and (not self.split_lock_path(chrom=chrom).exists())
    
    def index_bai_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_name}_{chrom}.bai'
    def index_lock_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_name}_{chrom}.index.lock'
    def index_complete(self, chrom): return (self.index_bai_path(chrom=chrom).exists()) and (not self.index_lock_path(chrom=chrom).exists())

'''
def run(sample_name, chrom): 
    bam_path = pathlib.Path(f'scratch/mcl4001/cfDNA/data/Sample_{sample_name}/bam/{sample_name}.final.bam')
    bai_path = bam_path.with_suffix('.bai')
    
    
    out_dir = bam_path.parent/'bam_per_chrom'
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_out_path = out_dir/f'{sample_name}_{chrom}.bam'
    bai_out_path = out_dir/f'{sample_name}_{chrom}.bai'
    
    lock_path = out_dir/f'{sample_name}_{chrom}.lock'
    
    if (bam_out_path.exists()) and (bai_out_path.exists()) and (not lock_path.exists()): 
        logger.info(f'{sample_name}: Already finished splitting {chrom}')
        return
    with job_lock.JobLock(lock_path, outputfiles=[bam_out_path, bai_out_path]) as lock: 
        if not lock: 
            logger.info(f'{sample_name}: Another process is already splitting {chrom}')
            return

        
        with job_lock.slurm_rsync_output(bam_out_path) as bam_out_path, open(bam_out_path, 'xb') as bam_file, job_lock.slurm_rsync_output(bai_out_path) as bai_out_path, open(bai_out_path, 'xb') as bai_file:            
            logger.info(f'{sample_name}: Splitting {chrom}')
            bam_path = job_lock.slurm_rsync_input(bam_path)
            bai_path = job_lock.slurm_rsync_input(bai_path)
            split_command = ['samtools', 'view', '-b', f'{str(bam_path)}', f'{chrom}']
            split_result = subprocess.run(split_command, check=True, stdout=bam_file)
            
            logger.info(f'{sample_name}: Indexing {chrom}')
            index_command = ['samtools', 'index', '-@', '1', f'{str(bam_out_path)}']
            index_result = subprocess.run(index_command, check=True, stdout=bai_file)
        return
    
'''

def split(sample, chrom, threads=16):
    if sample.split_complete(chrom=chrom): 
        logger.info(f'{sample.sample_name}: Split {chrom} already finished running successfully')
        
    with job_lock.JobLock(sample.split_lock_path(chrom=chrom), outputfiles=[sample.split_bam_path(chrom=chrom)]) as lock: 
        if not lock: return
        logger.info(f'{sample.sample_name}: Splitting {chrom}')        
        
        # bai_path = job_lock.slurm_rsync_input(sample.bai_path) # need these rsynced for the split
        # bam_path = job_lock.slurm_rsync_input(sample.bam_path)
        bai_path = sample.bai_path
        bam_path = sample.bam_path
        
        # with job_lock.slurm_rsync_output(sample.split_bam_path(chrom=chrom)) as split_bam_path:
        split_bam_path = sample.split_bam_path(chrom=chrom)
        bam_tools.split_bam_per_chrom(bam_path=bam_path, bai_path=bai_path, chrom=chrom, out_path=split_bam_path, threads=threads)
    return

def index(sample, chrom, threads=16):
    if not sample.split_complete(chrom=chrom):
        logger.info(f'{sample.sample_name}: split {chrom} has not finished running successfully')
        return
    if sample.index_complete(chrom=chrom):
        logger.info(f'{sample.sample_name}: index {chrom} already finished running successfully')
        return
    
    with job_lock.JobLock(sample.index_lock_path(chrom=chrom), outputfiles=[sample.index_bai_path(chrom=chrom)]) as lock: 
        if not lock: return
        logger.info(f'{sample.sample_name}: index {chrom} running')
        split_bam_path = sample.split_bam_path(chrom=chrom)
        # with job_lock.slurm_rsync_output(sample.index_bai_path(chrom=chrom)) as index_bai_path:         
        index_bai_path = sample.index_bai_path(chrom=chrom)
        bam_tools.index_bam(bam_path=split_bam_path, out_path=index_bai_path, threads=threads)
    return

    
def main():
    #sample_names = ['test', 'DTB-097-Progression-cfDNA']
 
    

    sample_names = ['DTB-097-Progression-cfDNA'] #sys.argv[1]
    samples = [Sample(sample_name=sample_name) for sample_name in sample_names]

    chromosomes = [f'chr{i}' for i in range(1,23,1)] + ['chrX', 'chrY']

    for j, sample in enumerate(samples): 
        #if j>0: break
        for i, chrom in enumerate(chromosomes):
            #if chrom!='chr21': continue
            split(sample=sample, chrom=chrom)
            index(sample=sample, chrom=chrom)
    return

if __name__=='__main__': 
    main()
        
        
            
    

