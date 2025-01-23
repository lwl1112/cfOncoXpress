from modules import misc_tools

import collections
import contextlib
import math
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib
import peakutils
import pysam
import re
import shutil
import subprocess
import tempfile
import time


def split_bam_per_chrom(bam_path, bai_path, chrom, out_path, threads=None, ):
    split_command = ['samtools', 'view', 
                     f'{bam_path}',
                     '-X', f'{bai_path}',
                     f'{chrom}',
                     '-b', # output in bam format
                     '-o', f'{out_path}', 
                    ]
    if threads is not None: split_command += ['-@', f'{threads}']
    subprocess.run(split_command, check=True)
    return

def index_bam(bam_path, out_path, threads=None): 
    index_command = ['samtools', 'index',
                     f'{bam_path}',
                     '-b', # output in bai format
                     f'{out_path}',
                    ]
    if threads is not None: index_command += ['-@', f'{threads}']
    subprocess.run(index_command, check=True)
    return
    
    
    
def get_coverage_memmap(bam_path, chrom, size, out_path, read_quality_threshold=20, only_save_to_midpoint=False, min_fragment_length=None, max_fragment_length=None):
    memmap = np.memmap(out_path, mode='w+', dtype=np.int32, shape=(size,))
    memmap[:] = 0
    
    bam = pysam.AlignmentFile(bam_path)
    for i, (read1, read2) in enumerate(read_pair_generator(bam, region_string=chrom)):
        if read1.mapping_quality<read_quality_threshold or read2.mapping_quality<read_quality_threshold: continue

        pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
        min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
        fragment_length = max_ind - min_ind
        
        if (min_fragment_length is not None) and (fragment_length<min_fragment_length): continue
        if (max_fragment_length is not None) and (fragment_length>max_fragment_length): continue
        
        if not only_save_to_midpoint: 
            memmap[min_ind:max_ind+1] += 1
        else:
            midpoint = (min_ind + max_ind) // 2
            memmap[midpoint] += 1
            
    memmap.flush()
    return memmap

def get_coverage_memmap_by_fragment_length(bam_path, chrom, size, out_path,  min_fragment_length, max_fragment_length, read_quality_threshold=20, only_save_to_midpoint=False):

    num_frag_lens = max_fragment_length - min_fragment_length
    
    memmap = np.memmap(out_path, mode='w+', dtype=np.uint8, shape=(size, num_frag_lens))
    memmap[:, :] = 0
    
    bam = pysam.AlignmentFile(bam_path)
    for i, (read1, read2) in enumerate(read_pair_generator(bam, region_string=chrom)):
        if read1.mapping_quality<read_quality_threshold or read2.mapping_quality<read_quality_threshold: continue

        pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
        min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
        fragment_length = max_ind - min_ind
        
        if (fragment_length<min_fragment_length) or (fragment_length>max_fragment_length): continue
        
        col = num_fragment_lens - (max_fragment_length - fragment_length) 
        
        if not only_save_to_midpoint: 
            memmap[min_ind:max_ind+1, col] += 1
        else:
            midpoint = (min_ind + max_ind) // 2
            memmap[midpoint, col] += 1
            
    memmap.flush()
    return memmap


    
def get_read_mapping_quality_counts(bam_path, out_path, region=None): 
    read_quality_counter = collections.Counter()
    
    bam = pysam.AlignmentFile(bam_path)    
    for i, (read1, read2) in enumerate(read_pair_generator(bam, region_string=region)):
        read_quality_counter[read1.mapping_quality] += 1
        read_quality_counter[read2.mapping_quality] += 1
    
    
    read_quality, counts = list(read_quality_counter.keys()), list(read_quality_counter.values())
    df = pd.DataFrame({'read_mapping_quality':read_quality, 'count': counts})
    df.to_csv(out_path, sep='\t', header=True, index=False)
    return df

def get_fragment_length_counts(bam_path, out_path, region=None): 
    fragment_length_counter = collections.Counter()
    
    bam = pysam.AlignmentFile(bam_path)    
    for i, (read1, read2) in enumerate(read_pair_generator(bam, region_string=region)):
        pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
        min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
        fragment_length = max_ind - min_ind        
        
        fragment_length_counter[fragment_length] += 1
    
    fragment_lengths, counts = list(fragment_length_counter.keys()), list(fragment_length_counter.values())
    df = pd.DataFrame({'fragment_length':fragment_lengths, 'count': counts})
    df.to_csv(out_path, sep='\t', header=True, index=False)
    return df

def get_depth(bam_path, out_path, chrom=None, filtered=False):
    command = ['samtools', 'depth', 
               f'{bam_path}',
               '-o', f'{out_path}', 
               '-r', f'{chrom}', 
              ]  
    if chrom is not None: command += ['-r', f'{chrom}'] # only look at the current chromosome
    if not filtered: command += ['-g', '1796',] # include UNMAP, SECONDARY, QCFAIL, and DUP
    else: command += ['-G', '3844']  # exclude UNMAP, SECONDARY, QCFAIL, DUP, and SUPPLEMENTARY
    subprocess.run(command, check=True)
    return

@contextlib.contextmanager
def split_bam(bam_path, chrom):
    bam_path = pathlib.Path(bam_path)
    split_cmd = ['samtools', 'view', f'{bam_path}', f'{chrom}', '-b']
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        bam_chrom_path = temp_dir/bam_path.with_suffix(f".{chrom}.bam").name
        with open(bam_chrom_path, 'w+') as f:  
            subprocess.check_call(split_cmd, stdout=f) # save output to temp file. 
        bai_chrom_path = bam_chrom_path.with_suffix('.bam.bai')
        index_cmd = ['samtools', 'index', '-@', '16', f'{bam_chrom_path}', f'{bai_chrom_path}']
        yield bam_chrom_path, bai_chrom_path
    
def read_pair_generator_with_qname(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = collections.defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            prev_read1, prev_read2 = read_dict.pop(qname)
            if read.is_read1:
                if prev_read1 is not None:
                    raise ValueError(f"Multiple read1s with the same qname: {qname}")
                yield qname, read, prev_read2
            else:
                if prev_read2 is not None:
                    raise ValueError(f"Multiple read2s with the same qname: {qname}")
                yield qname, prev_read1, read
                
def read_pair_generator(bam, region_string=None):
    for qname, start, end in read_pair_generator_with_qname(bam, region_string=region_string):
        yield start, end

def read_pair_qc(read1, read2, quality_thresh=10):
    quality = True if read1.mapping_quality + read1.mapping_quality > 2*quality_thresh else False
    not_translocated = True if read1.reference_id == read2.reference_id else False
    return quality, not_translocated


def randomly_split_bam_into_equal_parts(bam_path:str, out_paths:list, quality_threshold=0, random_seed=123456):
    num_splits = len(out_paths)
    bam = pysam.AlignmentFile(bam_path, 'rb')
    file_tracker = collections.Counter()
    # open output bam files for writing:
    out_files = [pysam.AlignmentFile(out_path, "wb", template=bam) for out_path in out_paths]
    for qname, read1, read2 in read_pair_generator_with_qname(bam, region_string=None):
        # remove low quality reads, if the user wants
        if (read1.mapping_quality + read1.mapping_quality) <= (2*quality_threshold): continue

        # gives a random number between [0,1)
        qname_float = misc_tools.str_to_float(s=qname)
        random_number = math.sin(123456*qname_float+random_seed) * 100000000 % 1000 / 1000
        foi = int(random_number // (1/num_splits))
        file_name = out_paths[foi].stem
        file = out_files[foi]
        file.write(read1)
        file.write(read2)
        file_tracker[file_name] += 2
        
    bam.close()
    for file in out_files: file.close()
    return file_tracker
        
        
        
    
    
    
    
'''
def get_coverage(in_bam_path:str, chrom:str="chr1", quality_thresh:int=10, max_iterations:int=None, print_progress:bool=None) -> dict:
    bam = pysam.AlignmentFile(in_bam_path)
    i = 0
    cov = collections.Counter()
    for read1, read2 in read_pair_generator(bam, region_string=chrom):
        i += 1
        if max_iterations is not None and  i>max_iterations: break
        if print_progress and i%10000==0: print(i, end = "\t")
        if not np.all(read_pair_qc(read1=read1, read2=read2, quality_thresh = quality_thresh)): continue
            
        pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
        min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
        cov.update(list(np.arange(min_ind, max_ind)))
    return cov
'''

def get_long_read_epsilon(in_bam_path:str, chrom:str="chr1", quality_thresh:int=10, max_iterations:int=None, print_progress:bool=True) -> dict:
    i = 0
    epsilons = collections.Counter()
    bam = pysam.AlignmentFile(in_bam_path, "rb")
    for read1, read2 in read_pair_generator(bam, region_string=chrom):
        i += 1
        if max_iterations is not None and  i>max_iterations: break
        if print_progress and i%10000==0: print(i, end = "\t")
        if not np.all(read_pair_qc(read1=read1, read2=read2, quality_thresh = quality_thresh)): continue
    
        pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
        min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
        insert_size = (max_ind-min_ind) + 1
        if insert_size > 167: 
            epsilon = insert_size-167
            epsilons.update([epsilon])
    return epsilons

def get_midpoint_avg_insert_size(in_bam_path:str, chrom:str="chr1", quality_thresh:int=10, max_iterations:int=None, print_progress:bool=True) -> dict:
    i = 0
    midpoints = collections.defaultdict(lambda: [])
    bam = pysam.AlignmentFile(in_bam_path, "rb")
    for read1, read2 in read_pair_generator(bam, region_string=chrom):
        i += 1
        if max_iterations is not None and  i>max_iterations: break
        if print_progress and i%10000==0: print(i, end = "\t")
        if not np.all(read_pair_qc(read1=read1, read2=read2, quality_thresh = quality_thresh)): continue
    
        pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
        min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
        midpoint, insert_size = (max_ind+min_ind) / 2, (max_ind-min_ind) + 1
        midpoints[midpoint].append(insert_size)
    midpoint_avg = collections.Counter()
    for midpoint in midpoints:
        avg_insert_size = np.mean(midpoints[midpoint])
        midpoint_avg[midpoint] = avg_insert_size
    return midpoint_avg

def plot_midpoint_insert_size_with_nuc_pos(nuc_pos:dict, midpoints:dict, chrom:str, start:int, end:int, sample_name:str=None, out_path:str=None) -> None:
    '''
    Goal: Plot the coverage with positions of nucleosomes identified by nucleoATAC labeled with verticle lines, and their area labeled with transparent rectangles.
    '''
    x, v_nuc, y_mid = [], [], []
    for i in range(start, end):
        x.append(i)
        if nuc_pos[i] != 0: v_nuc.append(i)
        y_mid.append(midpoints[i])
    plt.rcParams["figure.figsize"] = 10, 10
    plt.rcParams["font.size"] = 15
    plt.figure()
    plt.scatter(x, y_mid, color="black")
    # plt.axhline(167, color="blue")
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    for i in v_nuc:
        plt.axvline(x=i, color="red")
        rect = patches.Rectangle((i-(167/2), ymin),167,ymax-ymin,linewidth=0, facecolor='red', alpha=0.2)
        ax.add_patch(rect)
    plt.xlim(start, end)
    plt.ylim(0, ymax)
    plt.xlabel("Reference Position")
    plt.ylabel("Average Insert Size (bp)")
    plt.title(sample_name+" Average Insert Size at Read Midpoints Along "+chrom+":"+str(start)+"-"+str(end)+".", fontsize=12)
    plt.xticks(rotation=45)
    plt.show()
    if out_path is not None: plt.savefig(out_path+sample_name+"_midpoint_insert_size_with_nuc_pos_"+chrom+":"+str(start)+"-"+str(end)+".png")
    plt.close()


def get_coverage(bam_path:str, region:str="chr1", quality_thresh:int=0, max_iterations:int=None, print_progress:bool=True, sample_name:str=None) -> tuple:
    chrom = region.split(':')[0]
    if sample_name is not None: print("Getting coverage for sample {} at time: {}.".format(sample_name, time.ctime()))
    bam = pysam.AlignmentFile(bam_path, "rb")
    i = 0
    coverage = collections.Counter()
    for read1, read2 in read_pair_generator(bam, region_string=region):
        i += 1
        if max_iterations is not None and  i>max_iterations: break
        if not np.all(read_pair_qc(read1=read1, read2=read2, quality_thresh = quality_thresh)): continue
        pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
        min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
        coverage.update(list(np.arange(min_ind, max_ind)))
        if print_progress and i%100000==0: print('Found the {}th read at positions {}:{}-{} at time {}.'.format(i, chrom, min_ind, max_ind, time.ctime()))
    return coverage

def get_num_duplicate_reads(bam_path:str, region:str="chr1", print_progress:bool=True, sample_name:str=None) -> tuple:
    chrom = region.split(':')[0]
    bam = pysam.AlignmentFile(bam_path, "rb")
    i = 0
    dup_reads = collections.Counter()
    for read1, read2 in read_pair_generator(bam, region_string=region):
        i += 1
        if not np.all(read_pair_qc(read1=read1, read2=read2, quality_thresh = 0)): continue
        pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
        min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
        key = (min_ind,max_ind)
        dup_reads.update([key])
        if print_progress and i%100000==0: print('Found the {}th read at positions {}:{}-{} at time {}.'.format(i, chrom, min_ind, max_ind, time.ctime()))
    return dup_reads

def plot_insert_size(length_dict:dict, region:str, sample_name:str=None, out_path:str=None) -> None:
    plt.rcParams["figure.figsize"] = 10, 10
    plt.rcParams["font.size"] = 15
    
    items = sorted(length_dict.items())
    x = np.array([item[0] for item in items])
    y = np.array([item[1] for item in items])
    
    peaks_ind = list(peakutils.indexes(y, thres=0.1, min_dist=100))
    peaks_ind.sort(key=lambda ind: y[ind], reverse=True)
    peaks_ind = peaks_ind[0:4]

    plt.figure()
    plt.scatter(x, y, color="black")
    plt.scatter(x[peaks_ind], y[peaks_ind], color="orange")
    plt.xlabel("Insert Size (bp)")
    plt.ylabel("Count")
    plt.ylim(bottom=0)
    plt.title(sample_name + " Distribution Insert Size Along " + region + ".")
    for i in range(len(peaks_ind)): 
        text = "("+str(x[peaks_ind[i]])+", "+str(y[peaks_ind[i]])+")"
        plt.text(x=x[peaks_ind[i]], y=y[peaks_ind[i]], s=text) 
    if out_path is not None: plt.savefig(out_path + sample_name + "_insert_size_" + region + ".png")
    plt.show()
    plt.close()
    
    
def get_insert_size_and_quality(bam_path:str, region:str="chr1", quality_thresh:int=0, max_iterations:int=None, print_progress:bool=True) -> tuple:
    chrom = region.split(':')[0]
    bam = pysam.AlignmentFile(bam_path, "rb")
    i = 0
    length, quality = collections.Counter(), collections.Counter()
    for read1, read2 in read_pair_generator(bam, region_string=region):
        i += 1
        if max_iterations is not None and  i>max_iterations: break
        if not np.all(read_pair_qc(read1=read1, read2=read2, quality_thresh = quality_thresh)): continue
        pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
        min_ind, max_ind = min(pos1 + pos2), max(pos1 + pos2)
        length.update([(max_ind-min_ind)+1])
        if print_progress and i%100000==0: print('Found the {}th read at positions {}:{}-{} at time {}.'.format(i, chrom, min_ind, max_ind, time.ctime()))
        quality.update([read1.mapping_quality, read2.mapping_quality])
    return length, quality

def plot_coverage(coverage_dict:dict, region:str, sample_name:str=None, out_path:str=None) -> None:
    plt.rcParams["figure.figsize"] = 10, 10
    plt.rcParams["font.size"] = 15
    match = re.match(r".*:([0-9]+)-([0-9]+)", region)
    start_pos, end_pos = int(match.group(1)), int(match.group(2))
    x = np.arange(start_pos, end_pos)
    y = np.array([coverage_dict[i] for i in x])
    plt.figure()
    plt.scatter(x, y, color="black")
    plt.xlim(start_pos, end_pos)
    plt.ylim(bottom=0)
    plt.xlabel("Reference Position (bp)")
    plt.ylabel("Coverage")
    plt.ylim(bottom=0)
    plt.title(sample_name + " Coverage Along " + region + ".")
    if out_path is not None: plt.savefig(out_path + sample_name + "_coverage_" + region + ".png")
    plt.show()
    plt.close()

def plot_mapping_quality(quality_dict:dict, region:str, sample_name:str=None, out_path:str=None) -> None:
    plt.rcParams["figure.figsize"] = 10, 10
    plt.rcParams["font.size"] = 15
    plt.figure()
    plt.scatter(list(quality_dict.keys()), np.asarray(list(quality_dict.values()))/sum(quality_dict.values()), color="black")
    plt.xlabel("Mapping Quality")
    plt.ylabel("Relative Frequency")
    plt.ylim(bottom=0)
    plt.title(sample_name + " Distribution of Mapping Quality Along " + region + ".")
    if out_path is not None: plt.savefig(out_path + sample_name + "_mapping_quality_" + region + ".png")
    plt.show()
    plt.close()

def get_num_reads_per_chromosome(bam_path:str, print_progress:bool=True) -> dict:
    bam = pysam.AlignmentFile(bam_path)
    chromosomes = [f"chr{_}" for _ in range(1, 24)] + ["chrX", "chrY", "chrM"]
    num_reads_per_chromosome = collections.Counter()
    for i, region in enumerate(chromosomes):
        print(region)
        j = 0 
        for read1, read2 in read_pair_generator(bam, region_string=region):
            if print_progress and j%1000000==0: print(j, end="\t")
            num_reads_per_chromosome[region] += 1
            j+=1
        if print_progress: print()
    return num_reads_per_chromosome

def write_num_duplicates_to_bed(duplicates:dict, region:str, out_path:str, print_progress:bool=True) -> list:
    chrom = region.split(":")[0]
    dict_to_df = {'chromosome':[], 'start':[], 'end':[], 'number_duplicates':[]}
    for i, ((start, end), n_dups) in enumerate(duplicates.items()):
        if n_dups < 2: continue
        if print_progress and i%10000==0: print('{}% through at time {}.'.format(np.round((i/len(list(duplicates.keys())))*100, 2), time.ctime()))
        dict_to_df['chromosome'].append(chrom)
        dict_to_df['start'].append(start)
        dict_to_df['end'].append(end)
        dict_to_df['number_duplicates'].append(n_dups)
    dup_df = pd.DataFrame(dict_to_df)
    dup_df.to_csv(out_path+'duplicate_reads_{}.bed'.format(region), index=False, sep='\t')
    return dup_df
        