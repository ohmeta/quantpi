#!/usr/bin/env python

from quantpi.configer import metaconfig
from quantpi.configer import parse_yaml
from quantpi.configer import update_config
from quantpi.configer import custom_help_formatter

from quantpi.sampler import parse_samples
from quantpi.sampler import get_reads

from quantpi.simulator import parse_genomes
from quantpi.simulator import get_simulate_info
from quantpi.simulator import simulate_short_reads

from quantpi.tooler import parse
from quantpi.tooler import merge

from quantpi.qcer import change
from quantpi.qcer import compute_host_rate
from quantpi.qcer import qc_bar_plot
from quantpi.qcer import parse_fastp_json

from quantpi.aligner import flagstats_summary

from quantpi.profiler import profiler_init
from quantpi.profiler import get_all_abun_df
from quantpi.profiler import get_abun_df_bgi_soap
from quantpi.profiler import get_abun_df_bowtie2
from quantpi.profiler import get_abun_df_jgi
from quantpi.profiler import get_profile
from quantpi.profiler import metaphlan_init
from quantpi.profiler import merge_metaphlan_tables

from quantpi.profiler import genomecov_parse
from quantpi.profiler import genomecov_merge

from quantpi.__about__ import __version__, __author__

name = "quantpi"