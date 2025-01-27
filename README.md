# cfOncoPath
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
## Running the model
```
experiments.sh
```
## Correlation results for samples with matched RNA-seq
```
results.sh
```
## Figures for ISMB 2025
```
<a href='ismb2025figs.ipynb'>ismb2025figs.ipynb</a>
gsea_barplot.R
```
