# Marjorie Roskes (mcl4001@med.cornell.edu)

import methodtools
import numpy as np
import pandas as pd
import pathlib
import time

acceptable_ref_builds = ['GRCh37', 'GRCh38', 'hg19', 'hg38', 'dummy']
# dummy chromosomes are chromosomes of all the same size (small) used for testing purposes. 

chrom_dir = pathlib.Path(f'/athena/khuranalab/scratch/mcl4001/chromosomes/')


def chromosomes(ref_build):
    if ref_build not in acceptable_ref_builds:
        raise ValueError(f'{ref_build} is not an acceptable build. Acceptable builds are: {acceptable_ref_builds}.')
    if ref_build == 'hg19': ref_build = 'GRCh37'
    if ref_build == 'hg38': ref_build = 'GRCh38'
    df = pd.read_csv(chrom_dir/f'{ref_build}.txt', sep='\t', header=0)
    return df


def chrom_sizes(ref_build):
    if ref_build not in acceptable_ref_builds:
        raise ValueError(f'{ref_build} is not an acceptable build. Acceptable builds are: {acceptable_ref_builds}.')
    df = chromosomes(ref_build=ref_build)
    chroms, sizes = np.asarray(df['Chromosome']), np.asarray(df['Total length (bp)'])
    chrom_sizes_dict = {f'chr{chroms[i]}': sizes[i] for i in range(df.shape[0])}
    return chrom_sizes_dict

