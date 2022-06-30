#!/usr/bin/env snakemake

import sys
from pprint import pprint

import pandas as pd
import quantpi

from snakemake.utils import min_version
min_version("6.0")

shell.executable("bash")

QUANTPI_DIR = quantpi.__path__[0]
WRAPPER_DIR = os.path.join(QUANTPI_DIR, "wrappers")


IS_PE = True \
    if config["params"]["reads_layout"] == "pe" \
       else False


IS_INTERLEAVED = True \
    if config["params"]["interleaved"] \
       else False


HAVE_LONG = True \
    if IS_PE and config["params"]["have_long"] \
       else False


TRIMMING_DO = True \
    if config["params"]["trimming"]["oas1"]["do"] or \
       config["params"]["trimming"]["sickle"]["do"] or \
       config["params"]["trimming"]["fastp"]["do"] \
       else False


RMHOST_DO = True \
    if config["params"]["rmhost"]["soap"]["do"] or \ 
       config["params"]["rmhost"]["bwa"]["do"] or \
       config["params"]["rmhost"]["bowtie2"]["do"] or \
       config["params"]["rmhost"]["minimap2"]["do"] or \
       config["params"]["rmhost"]["kraken2"]["do"] or \
       config["params"]["rmhost"]["kneaddata"]["do"] \
       else False


if config["params"]["simulate"]["do"]:
    SAMPLES = quantpi.parse_genomes(config["params"]["samples"],
                                   config["output"]["simulate"])
else:
    SAMPLES = quantpi.parse_samples(config["params"]["samples"],
                                   config["params"]["interleaved"],
                                   config["params"]["reads_layout"],
                                   config["params"]["begin"])

SAMPLES_ID_LIST = SAMPLES.index.unique()

READS_FORMAT = "sra" \
    if "sra" in SAMPLES.columns \
       else "fastq"


KMCP_DBS = []
KMCP_TAXIDMAP = []
if config["params"]["profiling"]["kmcp"]["do"]["bacteriome"]:
    KMCP_DBS.append("bacteriome")
    KMCP_TAXIDMAP.append(os.path.join(config["params"]["profiling"]["kmcp"]["database"]["bacteriome"], "taxid.map"))
if config["params"]["profiling"]["kmcp"]["do"]["mycobiome"]:
    KMCP_DBS.append("mycobiome")
    KMCP_TAXIDMAP.append(os.path.join(config["params"]["profiling"]["kmcp"]["database"]["mycobiome"], "taxid.map"))
if config["params"]["profiling"]["kmcp"]["do"]["virome"]:
    KMCP_DBS.append("virome")
    KMCP_TAXIDMAP.append(os.path.join(config["params"]["profiling"]["kmcp"]["database"]["virome"], "taxid.map"))


include: "../rules/simulate.smk"
include: "../rules/raw.smk"
include: "../rules/trimming.smk"
include: "../rules/rmhost.smk"
include: "../rules/qcreport.smk"

include: "../rules/profiling_dna.smk"
include: "../rules/profiling_kmer.smk"
include: "../rules/profiling_marker.smk"
include: "../rules/profiling_genomecov.smk"

include: "../rules/profiling_custom.smk"

include: "../rules/profiling_function.smk"


rule all:
    input:
        rules.simulate_all.input,
        rules.raw_all.input,
        rules.trimming_all.input,
        rules.rmhost_all.input,
        rules.qcreport_all.input,
        rules.profiling_all.input