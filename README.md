# cfOncoXpress
<img width="1298" height="904" alt="image" src="https://github.com/user-attachments/assets/35dfe397-a334-4fdb-965b-7cb4e8085b69" />

## Requirements
python3
- job_lock
```
pip install slurm-python-utils
```

## Alignment
<a target='_blank' href='https://github.com/GavinHaLab/fastq_to_bam_paired_snakemake'>fastq_to_bam_paired_snakemake</a>
## Processing bam files to features per bp
```
cfDNA_split_bam_per_chrom.sh
cfDNA_coverage.sh
cfNDA_wps.sh
cfDNA_short_long_fragments.sh
```

## Generating cfDNA features in promoter regions
```
cfDNA_simpson_DTB-097-Progression.sh
cfDNA_features_DTB-097-Progression.sh
cfDNA_mat_DTB-097-Progression.sh
```
## Running the model: SVR regression + copy number correction
```
experiments.sh
```
## Correlation results for samples with matched RNA-seq
```
results.sh
```
## Figures
<a href='ismb2025figs.ipynb' target='_new'>ismb2025figs.ipynb</a>
```
gsea_barplot.R
```

![License](https://img.shields.io/badge/license-CC--BY--NC%204.0-lightgrey.svg)




