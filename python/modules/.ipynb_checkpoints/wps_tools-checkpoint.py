from modules import bam_tools, bed_tools, chrom_tools, logger_tools

import job_lock

import collections
import datetime
import errno
import methodtools
import numpy as np
import pandas as pd
import pathlib
import pysam
import re
import scipy.signal
import shutil
import time 

job_lock.JobLock.defaultcorruptfiletimeout = datetime.timedelta(hours=1)
job_lock.JobLock.setdefaultminimumtimeforiterativelocks(datetime.timedelta(minutes=10))
job_lock.JobLockAndWait.defaultsilent = True



logger = logger_tools.logger

class WPSBamFile:
    def __init__(self, bam_chrom, chrom, relevant_read_size=(0, np.infty), chunk_size:int=700):
        self.bam = bam_chrom # pysam AlignmentFile
        self.chrom = chrom
        self.relevant_read_size = relevant_read_size
        self.fragment_tracker = set()
        self.fragment_tracker_start = np.infty
        self.fragment_tracker_end = -np.infty
        self.chunk_size = chunk_size
    
    def add_to_fragment_tracker(self, start, end):
        relevant_read_size = self.relevant_read_size
        fragment_tracker_start, fragment_tracker_end = self.fragment_tracker_start, self.fragment_tracker_end
        # first check if and where it needs to be updated
        # case 1: (start, end) does not intersect with curr region in fragment tracker
        if (end < fragment_tracker_start) or (start > fragment_tracker_end):
            # --> reset the fragment tracker
            '''
            WARNING! 
            This case should only happen in testing when rerunning from start of testing region. 
            In real code, all new regions should be continuous. 
            '''
            self.fragment_tracker.clear()
            self.fragment_tracker_start = np.infty
            self.fragment_tracker_end = -np.infty
            fragment_tracker_start, fragment_tracker_end = self.fragment_tracker_start, self.fragment_tracker_end        
        
        # case 2: (start, end) extends on both sides of curr region in fragment tracker
        elif (start < fragment_tracker_start) and (end > fragment_tracker_end):
            # recursively update the fragment tracker on either side
            self.add_to_fragment_tracker(start, fragment_tracker_start)
            self.add_to_fragment_tracker(fragment_tracker_end, end)
            return

        left_side_covered, right_side_covered = False, False

        if start >= fragment_tracker_start:
            left_side_covered = True

        if end <= fragment_tracker_end:
            right_side_covered = True

        # case 3: (start, end) is a subset of curr region in fragment tracker (doesn't extend on either side)
        if left_side_covered and right_side_covered:
            return

        # case 4: (start, end) only extends on RHS of curr region in fragment tracker (LHS covered)
        if left_side_covered:
            region_start = fragment_tracker_end + 1
        else:
            region_start = max(start - relevant_read_size[1], 1)

        # case 5: (start, end) only extends on LHS of curr region in fragment tracker (RHS covered)
        if right_side_covered:
            region_end = fragment_tracker_start - 1
        else:
            region_end = end + relevant_read_size[1]

        # update fragment tracker
        region = f'{self.chrom}:{region_start}-{region_end}'
        fragment_tracker = self.fragment_tracker # mutable - add to fragment tracker adds to self.fragment_tracker
        for qname, read1, read2 in bam_tools.read_pair_generator_with_qname(self.bam, region_string=region):
            #if not np.all(bam_tools.read_pair_qc(read1=read1, read2=read2, quality_thresh = 0)): continue
            pos1, pos2 = read1.get_reference_positions(), read2.get_reference_positions()
            min_ind, max_ind = min(pos1 + pos2)+1, max(pos1 + pos2)+1 # return 0-based pos1, pos2 to 1-based
            fragment_length = max_ind-min_ind
            if (fragment_length < relevant_read_size[0]) or (fragment_length > relevant_read_size[1]): 
                # fragment is too small or too big
                continue
            else: 
                fragment_tracker.add((qname, min_ind, max_ind))

        # reset the fragment tracker limits
        self.fragment_tracker_start = min(fragment_tracker_start, start)
        self.fragment_tracker_end = max(fragment_tracker_end, end)
        
    def get_wps(self, start, end):
        chunk_size = self.chunk_size
        self.add_to_fragment_tracker(start, ((end + chunk_size) // chunk_size) * chunk_size)

        relevant_read_size = self.relevant_read_size # create a link
        fragment_tracker = self.fragment_tracker # create a link so looping is faster
        
        fragments_to_remove = set()
        wps = 0
        for fragment in fragment_tracker: 
            qname, min_ind, max_ind = fragment
            if (max_ind-min_ind < relevant_read_size[0]) or (max_ind-min_ind > relevant_read_size[1]): 
                # fragment is too small or too big
                continue
            if (min_ind <= start) and (max_ind >= end): 
                # fragment spans the whole window
                wps += 1
            elif (max_ind >= start) and (min_ind <= end): 
                # fragment doesn't span the whole window but one endpoint overlaps 
                wps -= 1
            elif max_ind < start:
                # fragment is from previous window and doesn't intersect current window at all --> remove
                # fragment completely to the left of window
                fragments_to_remove.add(fragment)
            else:
                # fragment completely to right of window - still need to keep for subsequent windows
                continue
        fragment_tracker -= fragments_to_remove
        self.fragment_tracker_start = max(self.fragment_tracker_start, start)
        return wps
    
def get_positive_regions(wps_memmap, 
                         cov_memmap,
                         max_consec_neg=5, 
                         relative_region_thresh:float=None,
                         relative_region_size:int=500): 
    positive_regions = []
    pos_region_start = None
    num_consec_neg = 0
    for i in range(len(wps_memmap)):
        if pos_region_start is None: 
            if wps_memmap[i]>0: pos_region_start = i
        else: 
            if wps_memmap[i]>0: num_consec_neg = 0
            else: 
                num_consec_neg += 1
                if num_consec_neg > max_consec_neg:
                    positive_regions.append(WPS_Positive_Region(start=pos_region_start, 
                                                                end=i-num_consec_neg, 
                                                                wps_memmap=wps_memmap,
                                                                cov_memmap=cov_memmap, 
                                                                relative_region_thresh=relative_region_thresh,
                                                                relative_region_size=relative_region_size))
                    num_consec_neg = 0
                    pos_region_start = None
    if pos_region_start is not None: 
        positive_regions.append(WPS_Positive_Region(start=pos_region_start, 
                                                    end=len(wps_memmap)-1-num_consec_neg, 
                                                    wps_memmap=wps_memmap,
                                                    cov_memmap=cov_memmap,
                                                    relative_region_thresh=relative_region_thresh,
                                                    relative_region_size=relative_region_size))
    
    return positive_regions

class Sample:
    def __init__(self, sample_name, data_dir, out_dir, ref_build):
        self.name = sample_name
        self.data_dir = pathlib.Path(data_dir)
        self.out_dir = pathlib.Path(out_dir)
        self.ref_build = ref_build
        self.chrom_sizes = chrom_tools.chrom_sizes(ref_build=self.ref_build) # dict key: chrom, value: size 
        self.chroms = list(self.chrom_sizes.keys())

    @property
    def lwps_dir(self):
        path = self.out_dir/'wps'
        path.parent.mkdir(parents=True, exist_ok=True)
        return path
    
    @property
    def lwps_freq_assembly_path(self):
        path = self.lwps_dir/'raw'/f'frequencies'/f"wps_freq.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path        

    @property
    def lwps_freq_assembly_lock_path(self):
        return self.lwps_freq_assembly_path.with_suffix('.txt.lock')

    @property
    def lwps_normal_freq_assembly_path(self):
        path = self.lwps_dir/'normal'/f'frequencies'/f"wps_freq.normal.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path        

    @property
    def lwps_normal_freq_assembly_lock_path(self):
        return self.lwps_normal_freq_assembly_path.with_suffix('.txt.lock')

    @property
    def lwps_smooth_freq_assembly_path(self):
        path = self.lwps_dir/'smooth'/f'frequencies'/f"wps_freq.smooth.normal.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path        

    @property
    def lwps_smooth_freq_assembly_lock_path(self):
        return self.lwps_smooth_freq_assembly_path.with_suffix('.txt.lock')

    @property
    def lwps_nucleosome_peaks_assembly_path(self):
        path = self.lwps_dir/'nucleosome_peaks'/f"nucleosome_peaks.bed"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path        

    @property
    def  lwps_nucleosome_peaks_assembly_lock_path(self):
        return self.lwps_nucleosome_peaks_assembly_path.with_suffix('.bed.lock')

    @property
    def lwps_clean_up_regions_to_assemble_lock_path(self):
        path = self.lwps_dir/f"{self.name}_clean_up_regions_to_assemble.lock"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path        

    @methodtools.lru_cache()
    def lwps_memmap(self, chrom):
        for sample_chrom in self.sample_chroms():
            if sample_chrom.chrom != chrom: continue
            return sample_chrom.lwps_memmap_chrom

    @methodtools.lru_cache()
    def lwps_normal_memmap(self, chrom):
        for sample_chrom in self.sample_chroms():
            if sample_chrom.chrom != chrom: continue
            return sample_chrom.lwps_normal_memmap_chrom
    
    @methodtools.lru_cache()
    def lwps_smooth_normal_memmap(self, chrom):
        for sample_chrom in self.sample_chroms():
            if sample_chrom.chrom != chrom: continue
            return sample_chrom.lwps_smooth_normal_memmap_chrom
    
    @methodtools.lru_cache()
    def lwps_freq(self, chromosomes=None):
        if chromosomes is None: chromosomes = self.chroms
        freq_counter = collections.Counter()
        for sample_chrom in self.sample_chroms():
            if sample_chrom.chrom not in chromosomes: continue
            logger.info(f'{self.name}: Loading WPS frequency along {sample_chrom.chrom}')
            freq_counter = sample_chrom.lwps_freq_chrom(counter_to_update=freq_counter)
        return freq_counter

    @methodtools.lru_cache()
    def lwps_normal_freq(self, chromosomes=None):
        if chromosomes is None: chromosomes = self.chroms
        freq_counter = collections.Counter()
        for sample_chrom in self.sample_chroms():
            if sample_chrom.chrom not in chromosomes: continue
            logger.info(f'{self.name}: Loading median-normalized WPS frequency along {sample_chrom.chrom}')
            freq_counter = sample_chrom.lwps_normal_freq_chrom(counter_to_update=freq_counter)
        return freq_counter

    def lwps_smooth_normal_freq(self, chromosomes=None):
        if chromosomes is None: chromosomes = self.chroms
        freq_counter = collections.Counter()
        for sample_chrom in self.sample_chroms():
            if sample_chrom.chrom not in chromosomes: continue
            logger.info(f'{self.name}: Loading smoothed median-normalized WPS frequency along {sample_chrom.chrom}')
            freq_counter = sample_chrom.lwps_smooth_normal_freq_chrom(counter_to_update=freq_counter)
        return freq_counter

    def lwps_positive_regions(self, chromosomes=None):
        if chromosomes is None: chromosomes = self.chroms
        positive_region_df = pd.DataFrame()
        for sample_chrom in self.sample_chroms():
            if sample_chrom.chrom not in chromosomes: continue
            logger.info(f'{self.name}: Loading WPS positive regions along {sample_chrom.chrom}, current data frame shape is {positive_region_df.shape}')
            positive_region_df = sample_chrom.lwps_positive_regions_chrom(df_to_update=positive_region_df)
        return positive_region_df

    @methodtools.lru_cache()
    def nucleosome_peaks(self):
        if not self.lwps_nucleosome_peaks_assembly_path.exists(): return # output file hasn't been made yet.
        nuc_peaks = pd.read_csv(self.lwps_nucleosome_peaks_assembly_path, 
                                sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])
        return nuc_peaks

    @methodtools.lru_cache()
    def nucleosome_peaks_BedFile(self):
        if not self.lwps_nucleosome_peaks_assembly_path.exists(): return # output file hasn't been made yet.
        bed_file = bed_tools.BedFile(sample_name=self.name,
                                     path=self.lwps_nucleosome_peaks_assembly_path, 
                                     header=None, 
                                     column_names=['chrom', 'start', 'end', 'score'])
        return bed_file

    def sample_chroms(self):
        sample_chroms = []
        for chrom in self.chroms:
            sample_chroms.append(Sample_Chrom(sample_name=self.name, 
                                              data_dir=self.data_dir, 
                                              out_dir=self.out_dir, 
                                              ref_build=self.ref_build, 
                                              chrom=chrom))
        return sample_chroms
    
    def chrom_wise_pipeline_run(self, wps_k, wps_tau, wps_lambda, wps_alpha, relevant_read_size, job_region_window, relative_region_thresh:float=None, relative_region_size:int=500, chromosomes:list=None):
        for sample_chrom in self.sample_chroms():
            if (chromosomes is not None) and (sample_chrom.chrom not in chromosomes): continue
            #job_lock.slurm_clean_up_temp_dir()
            # STEP 1. run wps.
            sample_chrom.run_wps_all_regions(job_region_window=job_region_window, wps_k=wps_k,
                                             relevant_read_size=relevant_read_size)
            
            # STEP 2. assemble wps.
            wps_memmap = sample_chrom.assemble_wps_memmap_chrom(wps_k=wps_k, job_region_window=job_region_window)
            
            # Once WPS memmap has been assembled, remove regions_to_assemble
            sample_chrom.clean_up_regions_to_assemble_chrom()
            
            # STEP 3. normalize wps.
            normal_wps_memmap = sample_chrom.median_normalize_wps_chrom(wps_k=wps_k, wps_tau=wps_tau)
            
            # STEP 4. smooth wps.
            smooth_wps_memmap = sample_chrom.smooth_wps_chrom(wps_k=wps_k, wps_tau=wps_tau, wps_lambda=wps_lambda)

            # STEP 5. mean local wps.
            nps_memmap = sample_chrom.nps_wps_chrom(wps_alpha=wps_alpha)
            
            # STEP 6. call nucleosome peaks. 
            sample_chrom.call_nucleosome_peaks_chrom(relative_region_thresh=relative_region_thresh, 
                                                     relative_region_size=relative_region_size)

    def run_wps_all_chroms(self, wps_k, relevant_read_size, job_region_window):
        for sample_chrom in self.sample_chroms():
            #job_lock.slurm_clean_up_temp_dir()
            sample_chrom.run_wps_all_regions(job_region_window=job_region_window, wps_k=wps_k,
                                             relevant_read_size=relevant_read_size)
            
    def assemble_wps_memmap_all_chroms(self, wps_k, job_region_window):
        for sample_chrom in self.sample_chroms():
            wps_memmap = sample_chrom.assemble_wps_memmap_chrom(wps_k=wps_k, job_region_window=job_region_window)

    def median_normalize_wps_all_chroms(self, wps_k, wps_tau):
        for sample_chrom in self.sample_chroms():
            normal_wps_memmap = sample_chrom.median_normalize_wps_chrom(wps_k=wps_k, wps_tau=wps_tau)
    
    def smooth_wps_all_chroms(self, wps_k, wps_tau, wps_lambda):
        for sample_chrom in self.sample_chroms():
            smooth_wps_memmap = sample_chrom.smooth_wps_chrom(wps_k=wps_k, wps_tau=wps_tau, wps_lambda=wps_lambda)
            
    def call_nucleosome_peaks_all_chroms(self, relative_region_thresh:float=None, relative_region_size:int=500):
        for sample_chrom in self.sample_chroms():
            sample_chrom.call_nucleosome_peaks_chrom(relative_region_thresh=relative_region_thresh,
                                                     relative_region_size=relative_region_size)
            
    def get_wps_freq(self):
        for sample_chrom in self.sample_chroms():
            wps_freq = sample_chrom.get_wps_freq_chrom()

    def assemble_wps_freq(self):
        if self.lwps_freq_assembly_path.exists() and not self.lwps_freq_assembly_lock_path.exists(): return
        with job_lock.JobLock(self.lwps_freq_assembly_lock_path, outputfiles=[self.lwps_freq_assembly_path]) as lock:
            if not lock: return #another process is already running this job
            data = collections.Counter()
            for sample_chrom in self.sample_chroms():
                logger.info(f'{self.name}: Loading WPS frequency along {sample_chrom.chrom} for assembly')
                chrom_data = sample_chrom.lwps_freq_chrom()
                if chrom_data is not None:
                    data.update(chrom_data)
                else: return
            with open(self.lwps_freq_assembly_path, 'w') as f:
                for wps, freq in data.items():
                    f.write(f'{wps}\t{freq}\n')
        return
    
    def get_wps_normal_freq(self):
        for sample_chrom in self.sample_chroms():
            normal_wps_freq = sample_chrom.get_wps_normal_freq_chrom()

    def assemble_wps_normal_freq(self):
        if self.lwps_normal_freq_assembly_path.exists() and not self.lwps_normal_freq_assembly_lock_path.exists(): return
        with job_lock.JobLock(self.lwps_normal_freq_assembly_lock_path, outputfiles=[self.lwps_normal_freq_assembly_path]) as lock:
            if not lock: return #another process is already running this job
            data = collections.Counter()
            for sample_chrom in self.sample_chroms():
                logger.info(f'{self.name}: Loading median-normalized WPS frequency along {sample_chrom.chrom} for assembly')
                chrom_data = sample_chrom.lwps_normal_freq_chrom()
                if chrom_data is not None:
                    data.update(chrom_data)
                else: return
            with open(self.lwps_normal_freq_assembly_path, 'w') as f:
                for wps, freq in data.items():
                    f.write(f'{wps}\t{freq}\n')
        return

    def get_wps_smooth_normal_freq(self):
        for sample_chrom in self.sample_chroms():
            smooth_normal_wps_freq = sample_chrom.get_wps_smooth_normal_freq_chrom()

    def assemble_wps_smooth_freq(self):
        if self.lwps_smooth_freq_assembly_path.exists() and not self.lwps_smooth_freq_assembly_lock_path.exists(): return
        with job_lock.JobLock(self.lwps_smooth_freq_assembly_lock_path, outputfiles=[self.lwps_smooth_freq_assembly_path]) as lock:
            if not lock: return #another process is already running this job
            data = collections.Counter()
            for sample_chrom in self.sample_chroms():
                logger.info(f'{self.name}: Loading smoothed, median-normalized WPS frequency along {sample_chrom.chrom} for assembly')
                chrom_data = sample_chrom.lwps_smooth_normal_freq_chrom()
                if chrom_data is not None:
                    data.update(chrom_data)
                else: return
            with open(self.lwps_smooth_freq_assembly_path, 'w') as f:
                for wps, freq in data.items():
                    f.write(f'{wps}\t{freq}\n')
        return

    def assemble_nucleosome_peaks(self):
        if self.lwps_nucleosome_peaks_assembly_path.exists() and not self.lwps_nucleosome_peaks_assembly_lock_path.exists():
            return
        with job_lock.JobLock(self.lwps_nucleosome_peaks_assembly_lock_path, outputfiles=[self.lwps_nucleosome_peaks_assembly_path]) as lock:
            if not lock: return #another process is already running this job
            data = pd.DataFrame()
            for sample_chrom in self.sample_chroms():
                logger.info(f'{self.name}: Loading nucleosome peaks along {sample_chrom.chrom} for assembly')
                chrom_data = sample_chrom.nucleosome_peaks_chrom(df_to_update=None, rm_na=True)
                if chrom_data is not None:
                    data = data.append(chrom_data)
                else: return
            logger.info(f'{self.name}: Writing assembled nucleosome peaks')
            with open(self.lwps_nucleosome_peaks_assembly_path, 'w') as f:
                for i, row in data.iterrows():
                    f.write(f"{row['chrom']}\t{row['nucleosome_peak_start']}\t{row['nucleosome_peak_end']}\t{row['nucleosome_peak_score']}\n")
        return
    
    def clean_up_regions_to_assemble(self): 
        if (not (self.lwps_dir/'regions_to_assemble').exists() and 
            not self.lwps_clean_up_regions_to_assemble_lock_path.exists()): return # already cleaned up
        with job_lock.JobLock(self.lwps_clean_up_regions_to_assemble_lock_path) as lock:
            if not lock: return #another process is already running this job
            logger.info(f'{self.name}: Attempting removal of regions_to_assemble')
            dir_to_remove = self.lwps_dir/'regions_to_assemble'
            try:
                dir_to_remove.rmdir()
            except FileNotFoundError:
                logger.info(f'{self.name}: Successfully removed regions_to_assemble')
            except OSError as e:
                if e.errno == errno.ENOTEMPTY:
                    logger.info(f'{self.name}: Failed to remove regions_to_assemble - directory NOT EMPTY!')
                else:
                    raise
            else:
                logger.info(f'{self.name}: Successfully removed regions_to_assemble')
        return
    
class Sample_Chrom(Sample): 
    def __init__(self, sample_name, data_dir, out_dir, ref_build, chrom):
        super().__init__(sample_name=sample_name, data_dir=data_dir, out_dir=out_dir, ref_build=ref_build)
        self.chrom = chrom

    @methodtools.lru_cache()
    @property
    def chrom_size(self):
        return self.chrom_sizes[self.chrom]
    
    #  def bam_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_name}_{chrom}.bam'
    #def bai_path(self, chrom): return self.bam_per_chrom_dir/f'{self.sample_name}_{chrom}.bai'
    # C000120.final.REF_chr21.bam #          C000120.final.REF_chrUn_JTFH01001077v1_decoy.bam    
    # CTRL-201_chr17.bam
    @property
    def bam_chrom_path(self):
        #file_names = list(self.data_dir.glob(f'*.REF_{self.chrom}.bam'))
        # file_names = list(self.data_dir.glob(f'*.REF_{self.chrom.replace("chr", "*")}.bam'))
        # file_names = [file_name for file_name in file_names if re.match(f'.*REF_[a-z]{{0,3}}{self.chrom.replace("chr", "")}.bam$', file_name.name) is not None]
        file_names = list(self.data_dir.glob(f'{self.name}_{self.chrom}.bam'))

        if len(file_names) != 1: 
            logger.info(f'{self.name}, {self.chrom}, {self.data_dir}.')
            logger.info(f'filenames: {file_names}.')
        assert len(file_names) == 1
        return file_names[0]

    @property
    def bai_chrom_path(self):
        #file_names = list(self.data_dir.glob(f'*.REF_{self.chrom}.bam.bai'))
        #if len(file_names) != 1: 
        #    logger.info(f'{self.name}, {self.chrom}, {self.data_dir}.')
        #    logger.info(f'filenames: {file_names}.')
        #assert len(file_names) == 1
        #return file_names[0]
        return self.bam_chrom_path.with_suffix('.bai')#('.bam.bai')
    
    @property
    def lwps_memmap_path(self):
        path = self.lwps_dir/f'raw'/f'wps_{self.chrom}.memmap'
        path.parent.mkdir(parents=True, exist_ok=True)
        return path
    
    @property
    def lwps_memmap_lock_path(self):
        return self.lwps_memmap_path.with_suffix('.memmap.lock')

    @property
    def lwps_normal_memmap_path(self):
        path = self.lwps_dir/f'normal'/f'wps_{self.chrom}.normal.memmap'
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def lwps_normal_memmap_lock_path(self):
        return self.lwps_normal_memmap_path.with_suffix('.memmap.lock')

    @property
    def lwps_smooth_normal_memmap_path(self):
        path =  self.lwps_dir/f'smooth'/f'wps_{self.chrom}.smooth.normal.memmap'
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def lwps_smooth_normal_memmap_lock_path(self):
        return self.lwps_smooth_normal_memmap_path.with_suffix(".memmap.lock")

    @property
    def lwps_nps_memmap_path(self):
        path =  self.lwps_dir/f'nps'/f'wps_{self.chrom}.nps.memmap'
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def lwps_nps_memmap_lock_path(self):
        return self.lwps_nps_memmap_path.with_suffix(".memmap.lock")    
    
    @property
    def lwps_positive_regions_path(self):
        path = self.lwps_dir/f'positive_regions'/f'wps_positive_regions_{self.chrom}.txt'
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def lwps_nucleosome_peaks_path(self):
        path = self.lwps_dir/f'nucleosome_peaks'/f'wps_nucleosome_peaks_{self.chrom}.txt'
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def coverage_memmap_path(self):
        path = self.out_dir/f'coverage'/f'coverage_{self.chrom}.memmap'
        return path
    
    @property
    def coverage_memmap_lock_path(self):
        return self.out_dir/f'stats_{self.chrom}.txt.lock'

    @property
    def lwps_nucleosome_peaks_lock_path(self):
        return self.lwps_nucleosome_peaks_path.with_suffix('.txt.lock')
    
    @property
    def lwps_clean_up_regions_to_assemble_chrom_lock_path(self):
        return self.lwps_dir/f"clean_up_regions_to_assemble_{self.chrom}.txt.lock"

    @property
    def lwps_freq_path(self):
        path = self.lwps_dir/'raw'/f'frequencies'/f"wps_freq_{self.chrom}.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path 

    @property
    def lwps_freq_lock_path(self):
        return self.lwps_freq_path.with_suffix('.txt.lock')

    @property
    def lwps_normal_freq_path(self):
        path = self.lwps_dir/f'normal'/f'frequencies'/f"wps_freq_{self.chrom}.normal.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def lwps_normal_freq_lock_path(self):
        return self.lwps_normal_freq_path.with_suffix('.txt.lock')
    
    @property
    def lwps_smooth_normal_freq_path(self):
        path = self.lwps_dir/f'smooth'/f'frequencies'/f"wps_freq_{self.chrom}.smooth.normal.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def lwps_smooth_normal_freq_lock_path(self):
        return self.lwps_smooth_normal_freq_path.with_suffix('.txt.lock')
 
    @property
    def lwps_memmap_chrom(self):
        if not self.lwps_memmap_path.exists(): return
        return np.memmap(self.lwps_memmap_path, dtype=np.int32, mode='r')

    @property
    def lwps_normal_memmap_chrom(self):
        if not self.lwps_normal_memmap_path.exists(): return
        return np.memmap(self.lwps_normal_memmap_path, dtype=np.int32, mode='r')

    @property
    def lwps_smooth_normal_memmap_chrom(self):
        if not self.lwps_smooth_normal_memmap_path.exists(): return
        return np.memmap(self.lwps_smooth_normal_memmap_path, dtype=np.float32, mode='r')

    def lwps_freq_chrom(self, counter_to_update=None):
        if counter_to_update is None: counter_to_update = collections.Counter()
        if not self.lwps_freq_path.exists(): return counter_to_update
        with open(self.lwps_freq_path, 'r') as f:
            for line in f: 
                cov, freq = int(line.split()[0]), int(line.split()[1])
                counter_to_update[cov] += freq
        return counter_to_update

    def lwps_normal_freq_chrom(self, counter_to_update=None):
        if counter_to_update is None: counter_to_update = collections.Counter()
        if not self.lwps_normal_freq_path.exists(): return counter_to_update
        with open(self.lwps_normal_freq_path, 'r') as f:
            for line in f: 
                cov, freq = int(line.split()[0]), int(line.split()[1])
                counter_to_update[cov] += freq
        return counter_to_update

    def lwps_smooth_normal_freq_chrom(self, counter_to_update=None):
        if counter_to_update is None: counter_to_update = collections.Counter()
        if not self.lwps_smooth_normal_freq_path.exists(): return counter_to_update
        with open(self.lwps_smooth_normal_freq_path, 'r') as f:
            for line in f: 
                cov, freq = float(line.split()[0]), int(line.split()[1])
                counter_to_update[cov] += freq
        return counter_to_update

    def lwps_positive_regions_chrom(self, df_to_update=None):
        if df_to_update is None: df_to_update = pd.DataFrame()
        if not self.lwps_positive_regions_path.exists(): return df_to_update
        positive_regions_chrom_df = pd.read_csv(self.lwps_positive_regions_path, names=['start','end'], sep='\t', dtype=np.int32)
        positive_regions_chrom_df.insert(loc=0, column='chrom', value=[self.chrom]*positive_regions_chrom_df.shape[0])
        df_to_update = df_to_update.append(positive_regions_chrom_df)
        return df_to_update

    def nucleosome_peaks_chrom(self, df_to_update=None, rm_na:bool=False):
        if df_to_update is None: df_to_update = pd.DataFrame()
        #if not self.lwps_nucleosome_peaks_path.exists(): return df_to_update
        print(self.lwps_nucleosome_peaks_path)
        if not self.lwps_nucleosome_peaks_path.exists() or self.lwps_nucleosome_peaks_lock_path.exists(): return None
        nuc_peaks_chrom_df = pd.read_csv(self.lwps_nucleosome_peaks_path, 
                                         sep='\t', 
                                         header=0, 
                                         dtype={'nucleosome_peak_start':np.int32, 
                                                'nucleosome_peak_end':np.int32,
                                                'positive_region_start':np.int32, 
                                                'positive_region_end':np.int32,
                                                'nucleosome_peak_score':np.float32})
        #nuc_peaks_chrom_df.insert(loc=0, column='chrom', value=[self.chrom]*nuc_peaks_chrom_df.shape[0])
        if rm_na: nuc_peaks_chrom_df = nuc_peaks_chrom_df.dropna(axis=0)
        df_to_update = df_to_update.append(nuc_peaks_chrom_df)
        return df_to_update


    
    
    def chrom_regions(self, job_region_window):
        num_intervals = (self.chrom_size // job_region_window)
        if self.chrom_size % job_region_window != 0: num_intervals += 1
        sample_chrom_regions = []
        for i in range(num_intervals): 
            sample_chrom_region = Sample_Chrom_Region(sample_name=self.name, 
                                                      data_dir=self.data_dir, 
                                                      out_dir=self.out_dir, 
                                                      ref_build=self.ref_build,
                                                      chrom=self.chrom, 
                                                      start=i*job_region_window, 
                                                      end=(i+1)*job_region_window)
            sample_chrom_regions.append(sample_chrom_region)
        return sample_chrom_regions

    def run_wps_all_regions(self, job_region_window, wps_k, relevant_read_size):
        if self.lwps_memmap_path.exists() and not self.lwps_memmap_lock_path.exists(): return # already assembled
        for sample_chrom_region in self.chrom_regions(job_region_window=job_region_window):
            sample_chrom_region_wps = sample_chrom_region.run_wps(wps_k=wps_k, relevant_read_size=relevant_read_size)
            
    def assemble_wps_memmap_chrom(self, job_region_window, wps_k):
        if self.lwps_memmap_path.exists() and not self.lwps_memmap_lock_path.exists(): return    #the output file already exists
        for sample_chrom_region in self.chrom_regions(job_region_window=job_region_window): # all regions finished running?
            if (not sample_chrom_region.lwps_path.exists()) or (sample_chrom_region.lwps_lock_path.exists()): 
                return # not all sample_chrom_regions have finished running
        with job_lock.JobLock(self.lwps_memmap_lock_path, outputfiles=[self.lwps_memmap_path]) as lock:
            if not lock: return  #another process is already running this job
            logger.info(f'{self.name}: Assembling WPS memmap along {self.chrom}')
            with job_lock.slurm_rsync_output(self.lwps_memmap_path) as lwps_memmap_path:
                memmap = np.memmap(lwps_memmap_path, dtype=np.int32, shape=(self.chrom_size,), mode="w+")
                pad_beginning, pad_end = wps_k//2, (wps_k//2)-1
                memmap[:pad_beginning] = 0
                memmap[pad_end:] = 0
                for sample_chrom_region in self.chrom_regions(job_region_window=job_region_window):
                    start = sample_chrom_region.start + pad_beginning
                    end = np.min([sample_chrom_region.end + pad_beginning, self.chrom_size-pad_end])
                    logger.info(f'{self.name}: Assembling WPS memmap of region {sample_chrom_region.region}')
                    wps_region = np.memmap(sample_chrom_region.lwps_path, 
                                           mode='r',
                                           shape=(end-start,), 
                                           dtype=np.int32)
                    memmap[start:end] = wps_region
                memmap.flush()
                return memmap
               
    def median_normalize_wps_chrom(self, wps_k, wps_tau):
        if (not self.lwps_memmap_path.exists()) or (self.lwps_memmap_lock_path.exists()): 
            return # previous step (assembly) hasn't finished running
        if self.lwps_normal_memmap_path.exists() and not self.lwps_normal_memmap_lock_path.exists(): 
            return #the output file already exists
        with job_lock.JobLock(self.lwps_normal_memmap_lock_path, outputfiles=[self.lwps_normal_memmap_path]) as lock:
            if not lock: return #another process is already running this job
            logger.info(f'{self.name}: Median Normalizing WPS along {self.chrom}') 
            with job_lock.slurm_rsync_output(self.lwps_normal_memmap_path) as lwps_normal_memmap_path:

                wps_memmap = np.memmap(job_lock.slurm_rsync_input(self.lwps_memmap_path), dtype=np.int32, shape=(self.chrom_size,), mode="r")
                normal_memmap = np.memmap(lwps_normal_memmap_path, dtype=np.int32, shape=(self.chrom_size,), mode="w+")
                normal_memmap[:wps_tau//2] = wps_memmap[:wps_tau//2]-np.median(wps_memmap[:wps_tau])
                normal_memmap[-(wps_tau//2):] = wps_memmap[-(wps_tau//2):]-np.median(wps_memmap[-wps_tau:])
                irange = len(normal_memmap)-wps_tau +1
                interval = irange // 10
                for i in range(irange):
                    if i % interval == 0: 
                        way_through = (i/irange) * 100
                        logger.info(f'{np.round(way_through, 2)}% of the way through {self.chrom}')
                    start, end = i, i+wps_tau
                    midpoint = (start+end) // 2
                    normal_memmap[midpoint] = wps_memmap[midpoint] - np.median(wps_memmap[start:end+1])
                normal_memmap.flush()
                return normal_memmap 
            
    def smooth_wps_chrom(self, wps_k, wps_tau, wps_lambda):
        if (not self.lwps_normal_memmap_path.exists()) or (self.lwps_normal_memmap_lock_path.exists()): 
            return # previous step (normalizing) hasn't finished running
        if self.lwps_smooth_normal_memmap_path.exists() and not self.lwps_smooth_normal_memmap_lock_path.exists(): 
            return #the output file already exists
        with job_lock.JobLock(self.lwps_smooth_normal_memmap_lock_path, outputfiles=[self.lwps_smooth_normal_memmap_path]) as lock:
            if not lock: return #another process is already running this job
            logger.info(f'{self.name}: Smoothing Median-Normalized WPS along {self.chrom}')
            normal_memmap = np.memmap(job_lock.slurm_rsync_input(self.lwps_normal_memmap_path), dtype=np.int32, shape=(self.chrom_size,), mode="r")
            with job_lock.slurm_rsync_output(self.lwps_smooth_normal_memmap_path) as lwps_smooth_normal_memmap_path:
                smooth_memmap = np.memmap(lwps_smooth_normal_memmap_path, dtype=np.float32, shape=(len(normal_memmap),), mode="w+")
                smooth_memmap[:] = scipy.signal.savgol_filter(x=normal_memmap, window_length=wps_lambda, polyorder=2)
                smooth_memmap.flush()
                return smooth_memmap

    def nps_wps_chrom(self, wps_alpha):
        # mean local wps
        if (not self.lwps_smooth_normal_memmap_path.exists()) or (self.lwps_smooth_normal_memmap_lock_path.exists()): 
            return # previous step (smoothing) hasn't finished running
        if self.lwps_nps_memmap_path.exists() and not self.lwps_nps_memmap_lock_path.exists(): 
            return #the output file already exists
        with job_lock.JobLock(self.lwps_nps_memmap_lock_path, outputfiles=[self.lwps_nps_memmap_path]) as lock:
            if not lock: return #another process is already running this job
            logger.info(f'{self.name}: Generating NPS along {self.chrom}')
            smooth_memmap = np.memmap(job_lock.slurm_rsync_input(self.lwps_smooth_normal_memmap_path), 
                                      dtype=np.float32, 
                                      shape=(self.chrom_size,), 
                                      mode="r")
            with job_lock.slurm_rsync_output(self.lwps_nps_memmap_path) as lwps_nps_memmap_path:
                nps_memmap = np.memmap(lwps_nps_memmap_path, dtype=np.float32, shape=(len(smooth_memmap),), mode="w+")
                nps_memmap[:] = 0
                pad = wps_alpha//2
                nps_memmap[pad:-pad] = np.convolve(smooth_memmap, np.ones(wps_alpha) / wps_alpha, "valid")
                start_pad_value, end_pad_value = nps_memmap[pad], nps_memmap[-pad-1]
                nps_memmap[0:pad] = start_pad_value
                nps_memmap[-pad:] = end_pad_value
                nps_memmap.flush()
                return nps_memmap
            
    def call_nucleosome_peaks_chrom(self, relative_region_thresh:float=None, relative_region_size:int=500):
        if (not self.coverage_memmap_path.exists()) or (self.coverage_memmap_lock_path.exists()): return
        if (not self.lwps_smooth_normal_memmap_path.exists()) or (self.lwps_smooth_normal_memmap_lock_path.exists()): return
        if (self.lwps_positive_regions_path.exists() 
            and self.lwps_nucleosome_peaks_path.exists() 
            and not self.lwps_nucleosome_peaks_lock_path.exists()): return #the output file already exists
        with job_lock.JobLock(self.lwps_nucleosome_peaks_lock_path, outputfiles=[self.lwps_positive_regions_path, self.lwps_nucleosome_peaks_path]) as lock:
            if not lock: return #another process is already running this job
            logger.info(f'{self.name}: Calling nucleosome peaks along {self.chrom}')
            smooth_normal_memmap = np.memmap(job_lock.slurm_rsync_input(self.lwps_smooth_normal_memmap_path), 
                                             dtype=np.float32, 
                                             mode="r")
            cov_memmap = np.memmap(job_lock.slurm_rsync_input(self.coverage_memmap_path), 
                                   shape=(self.chrom_size,), 
                                   dtype=np.int32, 
                                   mode="r")
            with job_lock.slurm_rsync_output(self.lwps_positive_regions_path) as lwps_positive_regions_path,\
                                  job_lock.slurm_rsync_output(self.lwps_nucleosome_peaks_path) as lwps_nucleosome_peaks_path: 
                with open(lwps_positive_regions_path, 'w')  as positive_regions_file, \
                                                  open(lwps_nucleosome_peaks_path, 'w') as nucleosome_peaks_file:
                    nucleosome_peaks_file.write(f'chrom\t'\
                                                f'nucleosome_peak_start\t'\
                                                f'nucleosome_peak_end\t'\
                                                f'positive_region_start\t'\
                                                f'positive_region_end\t'\
                                                f'nucleosome_peak_score\n')
                    positive_regions = get_positive_regions(wps_memmap=smooth_normal_memmap, 
                                                            max_consec_neg=5, 
                                                            cov_memmap=cov_memmap,
                                                            relative_region_thresh=relative_region_thresh,
                                                            relative_region_size=relative_region_size)
                    for i, positive_region in enumerate(positive_regions):
                        positive_region.positive_region_left = positive_regions[i-1] if i>0 else None
                        positive_region.positive_region_right = positive_regions[i+1] if i<len(positive_regions)-1 else None
                        positive_regions_file.write(f'{self.chrom}\t'\
                                                    f'{positive_region.start}\t'\
                                                    f'{positive_region.end}\n')
                        nucleosome_peaks_file.write(f'{self.chrom}\t'\
                                                    f'{positive_region.nucleosome_peak[0]}\t'\
                                                    f'{positive_region.nucleosome_peak[1]}\t'\
                                                    f'{positive_region.start}\t'\
                                                    f'{positive_region.end}\t'\
                                                    f'{positive_region.nucleosome_peak_score}\n')
    def clean_up_regions_to_assemble_chrom(self): 
        if (not (self.lwps_dir/'regions_to_assemble'/self.chrom).exists() and 
            not self.lwps_clean_up_regions_to_assemble_chrom_lock_path.exists()): return # already cleaned up
        with job_lock.JobLock(self.lwps_clean_up_regions_to_assemble_chrom_lock_path) as lock:
            if not lock: return #another process is already running this job
            logger.info(f'{self.name}: Confirming raw WPS memmap file for {self.chrom} has been generated successfully')
            if not self.lwps_memmap_path.exists() or self.lwps_memmap_lock_path.exists():
                logger.info(f'{self.name}: Could not confirm successful generation of all raw WPS memmap files')
                return
            logger.info(f'{self.name}: Attempting removal of regions_to_assemble/{self.chrom}')
            shutil.rmtree(self.lwps_dir/'regions_to_assemble'/self.chrom)
            logger.info(f'{self.name}: Successfully removed regions_to_assemble/{self.chrom}')
            
        return


    def get_wps_freq_chrom(self):
        if self.lwps_memmap_chrom is None: return
        self.lwps_freq_path.parent.mkdir(parents=True, exist_ok=True)
        if self.lwps_freq_path.exists() and not self.lwps_freq_lock_path.exists(): return #the output file already exists
        with job_lock.JobLock(self.lwps_freq_lock_path, outputfiles=[self.lwps_freq_path]) as lock:
            if not lock: return #another process is already running this job
            logger.info(f'{self.name}: Getting WPS frequency along {self.chrom}')
            counter_to_update = collections.Counter()
            counter_to_update.update(self.lwps_memmap_chrom)
            logger.info(f'{self.name}: Writing WPS frequency along {self.chrom}')
            with open(self.lwps_freq_path, 'w') as f: 
                for wps, freq in counter_to_update.items():
                    f.write(f'{wps}\t{freq}\n')
            return counter_to_update

    def get_wps_normal_freq_chrom(self):
        if self.lwps_normal_memmap_chrom is None: return
        self.lwps_normal_freq_path.parent.mkdir(parents=True, exist_ok=True)
        if self.lwps_normal_freq_path.exists() and not self.lwps_normal_freq_lock_path.exists(): return #the output file already exists
        with job_lock.JobLock(self.lwps_normal_freq_lock_path, outputfiles=[self.lwps_normal_freq_path]) as lock:
            if not lock: return #another process is already running this job
            logger.info(f'{self.name}: Getting median-normalized WPS frequency along {self.chrom}')
            counter_to_update = collections.Counter()
            counter_to_update.update(self.lwps_normal_memmap_chrom)
            logger.info(f'{self.name}: Writing median-normalized WPS frequency along {self.chrom}')
            with open(self.lwps_normal_freq_path, 'w') as f: 
                for wps, freq in counter_to_update.items():
                    f.write(f'{wps}\t{freq}\n')
            return counter_to_update


    def get_wps_smooth_normal_freq_chrom(self):
        if self.lwps_smooth_normal_memmap_chrom is None: return
        self.lwps_smooth_normal_freq_path.parent.mkdir(parents=True, exist_ok=True)
        if self.lwps_smooth_normal_freq_path.exists() and not self.lwps_smooth_normal_freq_lock_path.exists(): return #the output file already exists
        with job_lock.JobLock(self.lwps_smooth_normal_freq_lock_path, outputfiles=[self.lwps_smooth_normal_freq_path]) as lock:
            if not lock: return #another process is already running this job
            logger.info(f'{self.name}: Getting smoothed median-normalized WPS frequency along {self.chrom}')
            counter_to_update = collections.Counter()
            counter_to_update.update(np.rint(self.lwps_smooth_normal_memmap_chrom).astype(np.int32))
            logger.info(f'{self.name}: Writing smoothed median-normalized WPS frequency along {self.chrom}')
            with open(self.lwps_smooth_normal_freq_path, 'w') as f: 
                for wps, freq in counter_to_update.items():
                    f.write(f'{wps}\t{freq}\n')
            return counter_to_update


class Sample_Chrom_Region(Sample_Chrom):
    def __init__(self, sample_name, data_dir, out_dir, ref_build, chrom, start, end):
        super().__init__(sample_name=sample_name, data_dir=data_dir, out_dir=out_dir, ref_build=ref_build, chrom=chrom)
        self.start = start
        self.end = end
        
    @property
    def region(self):
        return f"{self.chrom}:{self.start}-{self.end}"
        
    @property
    def lwps_path(self):
        return self.lwps_dir/f'regions_to_assemble'/(self.chrom)/(f"wps_{self.region}.memmap")

    @property
    def lwps_lock_path(self):
        return self.lwps_dir/f'regions_to_assemble'/(self.chrom)/(f"wps_{self.region}.memmap.lock")

    def run_wps(self, wps_k, relevant_read_size):
        if self.lwps_path.exists() and not self.lwps_lock_path.exists(): return    #the output file already exists
        if not self.bam_chrom_path.exists() or not self.bai_chrom_path.exists(): return
        self.lwps_path.parent.mkdir(parents=True, exist_ok=True)
        with job_lock.JobLock(self.lwps_lock_path, outputfiles=[self.lwps_path]) as lock:
            if not lock: return  #another process is already running this job
            logger.info(f'{self.name}: Getting WPS of region {self.region}')
            job_lock.slurm_rsync_input(self.bai_chrom_path)
            bam = pysam.AlignmentFile(job_lock.slurm_rsync_input(self.bam_chrom_path))
            wps_bam_file = WPSBamFile(bam_chrom=bam, 
                                      chrom=self.chrom, 
                                      relevant_read_size=relevant_read_size)
            region_start, region_end = self.start, min(self.end, self.chrom_size-wps_k+1)
            if region_start > region_end: raise ValueError("ERROR: Region end is smaller than region start!!") 
            region_size = region_end-region_start
            wps = np.memmap(self.lwps_path, shape=(region_size,), mode='w+', dtype=np.int32)
            interval = region_size//10
            for i in range(region_size):
                if i%interval==0: 
                    logger.info(f'{self.name}: {(i/region_size):.0%} of the way through getting WPS of region {self.region}')
                start, end = region_start+i+1, region_start+i+1+wps_k
                wps[i] = wps_bam_file.get_wps(start=start, end=end)
            wps.flush()
            return wps
        
      
        
class WPS_Positive_Region(Sample_Chrom):
    def __init__(self, start, end, wps_memmap, cov_memmap, relative_region_thresh:float=None, relative_region_size:int=500):
        self.start = start
        self.end = end
        self.wps_memmap = wps_memmap
        self.relative_region_thresh = relative_region_thresh # if not None, accounts for coverage in relative region when calling nuc peaks.
        self.cov_memmap = cov_memmap 
        self.relative_region_size = relative_region_size
        self.positive_region_left = None
        self.positive_region_right = None
        
    @property
    def positive_region_size(self):
        return self.end-self.start+1
    @property
    def positive_region_wps(self):
        return self.wps_memmap[self.start:self.end+1]

    @property
    def valley_left_wps(self):
        if self.positive_region_left is None: return self.wps_memmap[:self.start]
        else: return self.wps_memmap[self.positive_region_left.end:self.start]
    @property
    def valley_right_wps(self):
        if self.positive_region_right is None: return self.wps_memmap[self.end:]
        else: return self.wps_memmap[self.end:self.positive_region_right.start]

    @methodtools.lru_cache()
    @property
    def nucleosome_peak(self): 
        # return (-1, -1) if there does not exist a nucleosome peak in this positive region. 
        if (self.positive_region_size < 50) or (self.positive_region_size > 450): return (-1, -1)
        renormalized_wps = self.positive_region_wps - np.median(self.positive_region_wps)
        sub_positive_regions = get_positive_regions(wps_memmap=renormalized_wps, 
                                                    max_consec_neg=0,
                                                    cov_memmap=self.cov_memmap,
                                                    relative_region_size=self.relative_region_size,
                                                    relative_region_thresh=self.relative_region_thresh)
        contiguous_sums = []
        for sub_positive_region in sub_positive_regions:
            potential_nuc_peak_start = sub_positive_region.start + self.start
            potential_nuc_peak_end = sub_positive_region.end + self.start
            if (self.positive_region_size >= 50) and (self.positive_region_size <= 150):
                if np.max(sub_positive_region.positive_region_wps) <= 0.9: 
                    contiguous_sum = 0
                else:
                    contiguous_sum = np.sum(sub_positive_region.positive_region_wps)
            elif (self.positive_region_size > 150) and (self.positive_region_size <= 450):
                if (sub_positive_region.positive_region_size >= 50) and (sub_positive_region.positive_region_size <= 150):
                    if np.max(sub_positive_region.positive_region_wps) <= 0.9:
                        contiguous_sum = 0
                    else: 
                        contiguous_sum = np.sum(sub_positive_region.positive_region_wps) 
                else:
                    contiguous_sum = 0
            else:
                assert False # we should never be able to go here. No silent bugs!
            potential_nuc_peak_coverage = np.sum(self.cov_memmap[potential_nuc_peak_start:potential_nuc_peak_end+1])
            if potential_nuc_peak_coverage == 0: 
                contiguous_sum = 0
            elif self.relative_region_thresh is not None:
                potential_nuc_peak_midpoint = (potential_nuc_peak_start + potential_nuc_peak_end) // 2
                outer_region_start = potential_nuc_peak_midpoint - (self.relative_region_size//2)
                outer_region_end = potential_nuc_peak_midpoint + (self.relative_region_size//2)
                outer_region_coverage = np.sum(self.cov_memmap[outer_region_start:outer_region_end+1])
                relative_coverage = potential_nuc_peak_coverage / outer_region_coverage
                if relative_coverage > self.relative_region_thresh: 
                    contiguous_sum = 0
                
            # if we made it here then the relative coverage in the region is sufficient. 
            contiguous_sums.append(contiguous_sum)
        if contiguous_sums and np.any(contiguous_sums):
            nucleosome_peak_start = sub_positive_regions[np.argmax(contiguous_sums)].start + self.start
            nucleosome_peak_end = sub_positive_regions[np.argmax(contiguous_sums)].end + self.start
            if nucleosome_peak_start==nucleosome_peak_end: 
                return (-1,-1) # don't allow peaks of len==1
            else: 
                return (nucleosome_peak_start, nucleosome_peak_end)
        else: return (-1, -1)
            
    @methodtools.lru_cache()
    @property
    def nucleosome_peak_score(self): 
        # return NaN if there does not exist a nucleosome peak in this positive region. 
        if np.all(np.asarray(self.nucleosome_peak)==-1): return np.nan
        if self.positive_region_left is None and self.positive_region_right is None:
            raise ValueError('Both left and right positive region neighbors are None!')
        nucleosome_peak_start, nucleosome_peak_end = self.nucleosome_peak
        nucleosome_peak_max_wps = np.max(self.wps_memmap[nucleosome_peak_start:nucleosome_peak_end+1])
        valley_left_min = np.min(self.valley_left_wps)
        valley_right_min = np.min(self.valley_right_wps)
        return nucleosome_peak_max_wps - ((valley_left_min+valley_right_min)/2)
    
        