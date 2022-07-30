[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![ohmeta-badge](https://img.shields.io/badge/install%20with-ohmeta-brightgreen.svg?style=flat)](http://anaconda.org/ohmeta)
[![PyPI version](https://badge.fury.io/py/quantpi.svg)](https://badge.fury.io/py/quantpi)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/quantpi/badges/downloads.svg)](https://anaconda.org/bioconda/quantpi)

# Microbiome profiling pipeline

## Overview
<div align=center><img width="800" height="120" src="docs/workflow.svg"></div>

## Installation
```
mamba install quantpi=0.2.0
# or
pip install quantpi=0.2.0
```

## Run
```
➤ quantpi --help

    ██████  ██    ██  █████  ███    ██ ████████ ██████  ██ 
   ██    ██ ██    ██ ██   ██ ████   ██    ██    ██   ██ ██ 
   ██    ██ ██    ██ ███████ ██ ██  ██    ██    ██████  ██ 
   ██ ▄▄ ██ ██    ██ ██   ██ ██  ██ ██    ██    ██      ██ 
    ██████   ██████  ██   ██ ██   ████    ██    ██      ██ 
       ▀▀                                                  

           Omics for All, Open Source for All

A general profiling system focus on robust microbiome research

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  print software version and exit

available subcommands:
  
    init         init project
    profiling_wf
                 metagenome-profiling pipeline
    sync         quantpi sync project
```

## Workflow list
```
➤ quantpi profiling_wf --list

Running quantpi profiling_wf:
snakemake --snakefile /home/jiezhu/toolkit/quantpi/quantpi/snakefiles/profiling_wf.smk --configfile ./config.yaml --cores 240 --rerun-incomplete --keep-going --printshellcmds --re
ason --until all --list

simulate_all
prepare_short_reads
prepare_short_reads_all
prepare_long_reads_all
prepare_reads_all
raw_fastqc_all
raw_report
raw_report_merge
raw_report_all
raw_all
trimming_oas1_all
trimming_sickle_all
trimming_fastp
trimming_fastp_multiqc
trimming_fastp_all
trimming_report
trimming_report_merge
trimming_report_all
trimming_all
rmhost_soap_all
rmhost_bowtie2_index
rmhost_bowtie2
rmhost_kraken2_all
rmhost_kneaddata_all
rmhost_alignment_report
rmhost_bwa_all
rmhost_bowtie2_all
rmhost_minimap2_all
rmhost_report
rmhost_report_merge
rmhost_report_all
rmhost_all
qcreport_summary
qcreport_plot
qcreport_all
profiling_kraken2
profiling_kraken2_krona_report
profiling_kraken2_combine_kreport
profiling_kraken2_combine_kreport_mpa
profiling_kraken2_all
profiling_bracken
profiling_bracken_merge
profiling_bracken_all
profiling_kmcp_search
profiling_kmcp_search_merge
profiling_kmcp_profile
profiling_kmcp_profile_merge
profiling_kmcp_all
profiling_metaphlan2_all
profiling_metaphlan3
profiling_metaphlan3_merge
profiling_metaphlan3_all
profiling_alignment_bowtie2
profiling_alignment_bam_postprocess
profiling_genomecov_gen_bed
profiling_genomecov_gen_cov
profiling_genomecov_gen_cov_merge
profiling_genomecov_all
profiling_genome_coverm
profiling_genome_coverm_merge
profiling_genome_coverm_all
profiling_custom_bgi_soap_all
profiling_custom_bowtie2_all
profiling_custom_jgi_all
profiling_humann2_all
profiling_humann3_config
profiling_humann3
profiling_humann3_postprocess
profiling_humann3_join
profiling_humann3_split_stratified
profiling_humann3_all
profiling_all
all
```

## Workflow
### profiling_kraken2_all
<div align=center><img width="800" height="500" src="docs/workflow_kraken2.svg"></div>

### profiling_bracken_all
<div align=center><img width="800" height="500" src="docs/workflow_bracken.svg"></div>

### profiling_kmcp_all
<div align=center><img width="800" height="500" src="docs/workflow_kmcp.svg"></div>

### profiling_genomecov_all
<div align=center><img width="800" height="500" src="docs/workflow_genomecov.svg"></div>

### profiling_genome_coverm_all
<div align=center><img width="800" height="500" src="docs/workflow_coverm.svg"></div>

### profiling_metaphlan3_all
<div align=center><img width="800" height="500" src="docs/workflow_metaphlan3.svg"></div>

### profiling_humann3_all
<div align=center><img width="800" height="500" src="docs/workflow_humann3.svg"></div>