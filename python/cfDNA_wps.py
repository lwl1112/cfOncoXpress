# Marjorie Roskes (mcl4001@med.cornell.edu)

from modules import wps_tools

import argparse 
import numpy as np
import pandas as pd
import pathlib

# init parser.
p = argparse.ArgumentParser()
p.add_argument("--wps-k", type=int, required=True)
p.add_argument("--relevant-read-size", type=int, nargs=2, required=False, default=(0, np.infty))
p.add_argument("--job-region-window", type=int, required=False, default=10000000)
p.add_argument("--wps-tau", type=int, required=True)
p.add_argument("--wps-lambda", type=int, required=True)
p.add_argument("--wps-alpha", type=int, required=True)
p.add_argument("--relative-region-thresh", type=float, required=False, default=None)
p.add_argument("--relative-region-size", type=int, required=False, default=500)

args = p.parse_args()

def load_sample_names(metadata_path):
    metadata = pd.read_csv(metadata_path, sep='\t', header=0)        
    sample_names = []
    for idx, row in metadata.iterrows():
        name = row['id']
        if 'control' in name: continue
        sample_names.append(name)
    return sample_names


path = pathlib.Path(f'herberts_samples.txt')
sample_names = load_sample_names(metadata_path=path)
#sample_names = ['DTB-097-Progression-cfDNA'] 

for sample_name in sample_names: 
    #if sample_name != 'C000120': continue
    print(sample_name)
    sample = wps_tools.Sample(sample_name=sample_name, 
                              data_dir= f'scratch/mcl4001/cfDNA/data/Sample_{sample_name}/bam/bam_per_chrom/',
                              out_dir= f'scratch/mcl4001/cfDNA/analysis/Sample_{sample_name}/', 
                              ref_build='GRCh38',
                             )
    sample.chrom_wise_pipeline_run(wps_k=args.wps_k, 
                                   wps_tau=args.wps_tau, 
                                   wps_lambda=args.wps_lambda, 
                                   wps_alpha=args.wps_alpha, 
                                   relevant_read_size=args.relevant_read_size, 
                                   job_region_window=args.job_region_window, 
                                   relative_region_thresh=args.relative_region_thresh,
                                   relative_region_size=args.relative_region_size)
    #sample.get_wps_freq()
    #sample.assemble_wps_freq()
    #sample.get_wps_normal_freq()
    #sample.assemble_wps_normal_freq()
    #sample.get_wps_smooth_normal_freq()
    #sample.assemble_wps_smooth_freq()
    sample.assemble_nucleosome_peaks()
    sample.clean_up_regions_to_assemble()

'''
###################################### DANGER ZONE #######################################
# clean up
for sample_name in sample_names: 
    sample = wps_tools.Sample(sample_name=sample_name, 
                              data_dir= f'scratch/mcl4001/cfDNA/data/Sample_{sample_name}/bam/bam_per_chrom/',
                              out_dir= f'scratch/mcl4001/cfDNA/analysis/Sample_{sample_name}/',
                              test_mode=args.test_mode)
    sample.clean_up_regions_to_assemble(print_progress=args.print_progress)
##########################################################################################
'''