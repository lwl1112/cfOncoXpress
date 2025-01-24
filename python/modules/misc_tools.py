# Marjorie Roskes (mcl4001@med.cornell.edu)

from modules import chrom_tools, logger_tools

import job_lock 

import hashlib
import loess.loess_1d
import matplotlib.colors
import mygene
import numba as nb
import numpy as np
import pandas as pd
import pathlib
import pybedtools
import random
import scipy.signal
import time

logger = logger_tools.logger

@nb.njit
def ratio(h1, h2):
    return h1 / np.where(h2, h2, 1e-99)

def argsortbyratio(h1, h2):
    r = ratio(h1, h2)
    nonzero = np.nonzero(h1.astype(np.bool_) | h2.astype(np.bool_))
    r = r[nonzero]
    return tuple([np.take_along_axis(_, np.argsort(r, axis=0), axis=0) for _ in nonzero])

def roc_from_hist(h1, h2, reverse=False):
    nonzero = argsortbyratio(h1, h2)
    if reverse: nonzero = nonzero[::-1]
    h1roc, h2roc = [0], [0]
    for coord in np.transpose(nonzero):
        h1roc.append(h1roc[-1]+h1[tuple(coord)])
        h2roc.append(h2roc[-1]+h2[tuple(coord)])
    return np.array(h1roc), np.array(h2roc)

def roc(y1, y2):
    roc1, roc2 = [], []
    for value in sorted({-np.inf, np.inf} | set(y1) | set(y2)):
        roc1.append(np.count_nonzero(y1 < value) / len(y1))
        roc2.append(np.count_nonzero(y2 < value) / len(y2))
    return np.array(roc1), np.array(roc2)

def bed_to_memmap(in_bed_df, memmap_path, memmap_size, index, memmap_dtype=np.int8):
    if index!=0 and index!=1:
        raise ValueError('Index must be 0 or 1.')
    memmap = np.memmap(memmap_path, mode='w+', shape=memmap_size, dtype=memmap_dtype)
    memmap[:] = 0
    for i, (idx, row) in enumerate(in_bed_df.iterrows()):
        start, end = row['start']-index, row['end']-index 
        if start < 0: 
            raise ValueError(f'start = {start} is negative. Confirm index is really {index}.')
        if end >= len(memmap): 
            raise ValueError(f'end = {end} is greater than len(memmap)={len(memmap)}.')
        memmap[start:end+1] = 1
    memmap.flush()
    return memmap

def memmap_starts_ends(memmap):
    index = 0
    m = memmap
    if m.dtype.name == "bool" or "uint" in m.dtype.name:
        m = m.astype(np.int8)
    derivative = m.astype(np.int8)[1:] - m.astype(np.int8)[:-1]
    starts, = np.where(derivative==1)
    ends, = np.where(derivative==-1)
    starts = list(starts+1+index)
    ends = list(ends+index)
    #handle edge cases
    if m[0] == 1: starts = [0+index] + starts
    if m[-1] == 1: ends = ends + [len(m)-1+index]
    return np.asarray(starts), np.asarray(ends)

def cross_correlation(a, b):
    fourier = np.fft.fft(a), np.fft.fft(b)
    crosspower = fourier[0] * np.conj(fourier[1])
    invfourier = np.fft.ifft(crosspower)
    return np.real(invfourier)

def get_refseq(genesymbols):
    info = mygene.MyGeneInfo()
    result = info.querymany(genesymbols, scopes="symbol", fields="refseq", species="human")
    return {_["query"]: _.get("refseq", None) for _ in result}

def histological_status_color(status:str, tumor_fraction=None, ar_status=None, plasma_source=None,sex=None):
    if status in ['no_cancer', 'normal']: 
        if plasma_source == 'frozen': color = 'blue'
        elif plasma_source == 'fresh': color = 'dodgerblue' if sex == 'M' else 'lightskyblue'
    elif status in ['localized']: color = 'green'
    elif status in ['localized_advanced']: color = 'darkorange'
    elif status in ['CRPC-Adeno', 'CRPC-AR', 'Adeno', 'AR']:         
        if tumor_fraction == 0: color = 'palevioletred'
        elif (tumor_fraction is None) or (tumor_fraction > 0): color = 'magenta'
        if ar_status is not None:
            if ar_status in ['negative', 'low']: color = 'palevioletred'
            elif ar_status in ['medium', 'high']: color = 'magenta'
    elif status in ['CRPC-NEPC (Adeno+NED)', 'NEPC (Adeno+NED)']: color = 'crimson'
    elif status in ['CRPC-NEPC', 'CRPC-NE', 'NEPC', 'NE']: 
        if tumor_fraction == 0: color = 'rosybrown'
        elif (tumor_fraction is None) or (tumor_fraction > 0): color = 'darkred'
    else: assert False
    return color

def molecular_status_color(status:str):
    if status in ['CRPC-Adeno', 'CRPC-AR', 'Adeno', 'AR']: color = list(np.asarray([229, 92, 90])/255)
    elif status in ['CRPC-WNT', 'WNT']: color = list(np.asarray([104, 161, 12])/255)
    elif status in ['CRPC-SCL', 'SCL']: color = list(np.asarray([185, 96, 253])/255)
    elif status in ['CRPC-NEPC', 'CRPC-NE', 'NEPC', 'NE']: color = list(np.asarray([75, 181, 185])/255)
    else: assert False
    return color

def reshape_memmap(memmap, bin_size): 
    idxtochop = memmap.shape[0]%bin_size
    num_bins = memmap.shape[0]//bin_size
    resized_memmap = memmap[0:-idxtochop].reshape(num_bins, bin_size)
    return resized_memmap
    
def fast_get_bins(chrom_size, bin_size): 
    x = np.arange(chrom_size)
    idxtochop = x.shape[0]%bin_size
    choppedx = x[0:-idxtochop]
    num_bins = int(choppedx.shape[0]/bin_size)
    xreshaped = choppedx.reshape(num_bins,bin_size)
    starts, ends = xreshaped[:,0], xreshaped[:,-1]
    return starts, ends

def get_bin_centromere_overlap_bool(bins, buffer = 1000000, chrom=None): 
    centromere_path = pathlib.Path(f'/athena/khuranalab/scratch/mcl4001/hg38_regions_to_exclude/ucsc_centromeres.bed')
    centromeres = pd.read_csv(centromere_path, sep='\t', usecols=[1,2,3], names=['chr', 'start', 'end'], header=0)
    centromeres = centromeres.sort_values(['chr', 'start', 'end'])
    centromeres['start'] -= buffer
    centromeres['end'] += buffer
    if chrom is not None: 
        centromeres = centromeres[centromeres['chr']==chrom]
        bins = bins[bins['chr']==chrom]

    bedA, bedB = pybedtools.BedTool.from_dataframe(bins[['chr', 'start', 'end']]), pybedtools.BedTool.from_dataframe(centromeres)
    bedA, bedB = bedA.sort(), bedB.sort()
    
    res = bedA.intersect(bedB, c=True)
    resdf = res.to_dataframe()
    resdf.rename(columns={'chrom' : 'chr',
                          'name' : 'count',                      
                         }, inplace=True)
    counts = resdf['count']
    bool_counts = np.asarray(counts>0)
    return bool_counts
    
def get_bin_gap_overlap_bool(bins, chrom=None): 
    gap_path = pathlib.Path(f'/athena/khuranalab/scratch/mcl4001/hg38_regions_to_exclude/ucsc_gaps.bed')
    gaps = pd.read_csv(gap_path, sep='\t', usecols=[1,2,3], names=['chr', 'start', 'end'], header=0)
    gaps = gaps.sort_values(['chr', 'start', 'end'])
    if chrom is not None: 
        gaps = gaps[gaps['chr']==chrom]
        bins = bins[bins['chr']==chrom]    

    bedA, bedB = pybedtools.BedTool.from_dataframe(bins[['chr', 'start', 'end']]), pybedtools.BedTool.from_dataframe(gaps)
    bedA, bedB = bedA.sort(), bedB.sort()
    
    res = bedA.intersect(bedB, c=True)
    resdf = res.to_dataframe()
    resdf.rename(columns={'chrom' : 'chr',
                          'name' : 'count',                      
                         }, inplace=True)
    counts = resdf['count']
    bool_counts = np.asarray(counts>0)
    return bool_counts    

def get_bin_black_listed_region_overlap_bool(bins, chrom=None): 
    blr_path = pathlib.Path(f'/athena/khuranalab/scratch/mcl4001/hg38_regions_to_exclude/encode_black_list.bed')
    blrs = pd.read_csv(blr_path, sep='\t', names=['chr', 'start', 'end'], header=None)
    blrs = blrs.sort_values(['chr', 'start', 'end'])

    if chrom is not None: 
        blrs = blrs[blrs['chr']==chrom]
        bins = bins[bins['chr']==chrom]    
    
    bedA, bedB = pybedtools.BedTool.from_dataframe(bins[['chr', 'start', 'end']]), pybedtools.BedTool.from_dataframe(blrs)
    bedA, bedB = bedA.sort(), bedB.sort()
    
    
    res = bedA.intersect(bedB, c=True)
    resdf = res.to_dataframe()
    resdf.rename(columns={'chrom' : 'chr',
                          'name' : 'count',                      
                         }, inplace=True)
    counts = resdf['count']
    bool_counts = np.asarray(counts>0)
    return bool_counts    

def get_bin_unmappable_region_overlap_bool(bins, mappability_threshold=0.5, chrom=None):     
    map_bedgraph_path = pathlib.Path(f'/athena/khuranalab/scratch/mcl4001/hg38_regions_to_exclude/mappability_by_chrom/mappability.sorted.bedgraph')
    if chrom is not None: 
        map_bedgraph_path = pathlib.Path(f'/athena/khuranalab/scratch/mcl4001/hg38_regions_to_exclude/mappability_by_chrom/mappability.{chrom}.bedgraph')
        bins = bins[bins['chr']==chrom]
    bedgraph = pybedtools.BedTool(map_bedgraph_path)
    
    binbed = pybedtools.BedTool.from_dataframe(bins[['chr', 'start', 'end']])
    binbed = binbed.sort()
    
    res = binbed.map(bedgraph, c=4, o='mean')
    resdf = res.to_dataframe()
    resdf.rename(columns={'chrom' : 'chr',
                          'name' : 'mappability',                      
                         }, inplace=True)
    mappability = resdf['mappability'].copy()
    mappability[mappability == '.'] = 0
    mappability = np.asarray(mappability).astype(np.float32)
    bool_mappability = mappability<mappability_threshold
    return bool_mappability    


def str_to_float(s):
    return int(hashlib.sha256(s.encode()).hexdigest(), 16)

def get_gc_content(bins, rsync=False):
    '''
    Compute gc content in regions of the genome. Return an array the length of the number of bins input, each index corresponding the GC content of GRCh38 in that region.
    bins: pd df of chr, start, end
    '''
    logger.info(f'Getting GC content in {bins.shape[0]} genomic regions')
    fasta_path = '/athena/khuranalab/scratch/mcl4001/hg38_reference/hg38.fa'
    if rsync: 
        fasta_path = job_lock.slurm_rsync_input(fasta_path)
        fasta_index_path = job_lock.slurm_rsync_input('/athena/khuranalab/scratch/mcl4001/hg38_reference/hg38.fa.fai')
    
    binbed = pybedtools.BedTool().from_dataframe(bins)
    result = binbed.nucleotide_content(fi=str(fasta_path), seq=True, fullHeader=True)
    column_names = ['chr', 'start', 'end', 'AT_fraction', 'GC_fraction', 'num_A', 'num_C', 'num_G', 'num_T', 'num_N', 'num_other', 'sequence_length', 'sequence']
    df = result.to_dataframe(disable_auto_names=True, comment='#',names=column_names)
    return np.asarray(df['GC_fraction'])
    
def gc_correct_signal(bins, signal, signal_quantile_restriction=[0,1], rsync=False):
    '''
    Given a signal across genomic regions (bins), correct the signal for GC content in the region.
    bins: pd df of chr, start, end
    signal: numeric vector of same length and number of rows in bins
    signal_quantile_restriction: quantile boundaries indicating how much of the signal to include in the median calculation. default is to use the entire signal.
    
    Note: Only compatible for GRCh38 reference!!!!
    '''
    
    gc_content = get_gc_content(bins=bins, rsync=rsync)
    lower_limit, upper_limit = np.nanquantile(signal, signal_quantile_restriction)
    median = np.median(signal[(signal>=lower_limit) & (signal<=upper_limit)])
    logger.info(f'Fitting LOESS curve to {len(signal)} points')
    gc, expected_signal, weights = loess.loess_1d.loess_1d(x=gc_content, y=signal)
    logger.info(f'LOESS curve fit complete')
    corrected_signal = (signal-expected_signal) + median    
    logger.info(f'Signal corrected for GC content')
    return corrected_signal


    
def get_n_random_named_colors(n, exclude_white=True):
    options = list(matplotlib.colors.CSS4_COLORS.keys())
    if exclude_white: options = [_ for _ in options if 'white' not in _]
    return random.sample(options, n)
    
    
def get_local_mean(array, window):
    return scipy.signal.convolve(array, np.ones(window) / window, 'same')


def get_sliding_local_fft(array, window_size, period, *, axis=-1, print_progress=False):
    fft = np.full_like(array, np.nan)
    freq = np.fft.fftfreq(window_size)
    (idx,), = np.where(freq == 1/period)
    if window_size > array.shape[axis]:
        raise ValueError(f"Window size {window_size} is too big for axis {axis} of array with shape {array.shape}")

    for i in range(array.shape[axis]-window_size):
        if print_progress and ((i+1) % 100 == 0 or i+1 == array.shape[axis]-window_size):
            print(i, "/", array.shape[axis]-window_size, time.ctime())

        array_slice = [slice(None)] * array.ndim
        array_slice[axis] = slice(i, i+window_size)
        array_slice = tuple(array_slice)

        output_fft_slice = [slice(None)] * array.ndim
        output_fft_slice[axis] = i + window_size//2
        output_fft_slice = tuple(output_fft_slice)

        window_fft_slice = [slice(None)] * array.ndim
        window_fft_slice[axis] = idx
        window_fft_slice = tuple(window_fft_slice)

        window_fft = np.fft.fft(array[array_slice], axis=axis)  * (2 / window_size)
        fft[output_fft_slice] = abs(window_fft[window_fft_slice])

    return fft