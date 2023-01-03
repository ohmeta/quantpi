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


SAMPLES = quantpi.parse_genomes(config["params"]["samples"],
                                config["output"]["simulate"])

SAMPLES_ID_LIST = SAMPLES.index.unique()


include: "../rules/simulate.smk"


rule all:
    input:
        rules.simulate_all.input