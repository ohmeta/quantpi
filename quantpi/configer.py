#!/usr/bin/env python

import argparse
import os
import sys
import shutil

from ruamel.yaml import YAML


def parse_yaml(yaml_file):
    yaml = YAML()
    with open(yaml_file, "r") as f:
        return yaml.load(f)


def update_config(yaml_file_old, yaml_file_new, yaml_content, remove=True):
    yaml = YAML()
    yaml.default_flow_style = False
    if remove:
        os.remove(yaml_file_old)
    with open(yaml_file_new, "w") as f:
        yaml.dump(yaml_content, f)


class metaconfig:
    """
    config project directory
    """

    sub_dirs = [
        "envs",
        "profiles",
        "results",
        "logs/00.simulate_short_reads",
        "logs/00.prepare_short_reads",
        "logs/00.prepare_long_reads",
        "logs/00.raw_fastqc",
        "logs/00.raw_fastqc_multiqc",
        "logs/00.raw_report",
        "logs/00.raw_report_refine",
        "logs/00.raw_report_merge",
        "logs/01.trimming_sickle",
        "logs/01.trimming_fastp",
        "logs/01.trimming_fastp_multiqc",
        "logs/01.trimming_trimmomatic",
        "logs/01.trimming_trimmomatic_multiqc",
        "logs/01.trimming_report",
        "logs/01.trimming_report_refine",
        "logs/01.trimming_report_merge",
        "logs/02.rmhost_bwa_index",
        "logs/02.rmhost_bwa",
        "logs/02.rmhost_bowtie2_index",
        "logs/02.rmhost_bowtie2",
        "logs/02.rmhost_minimap2_index",
        "logs/02.rmhost_minimap2",
        "logs/02.rmhost_kraken2",
        "logs/02.rmhost_kneaddata",
        "logs/02.rmhost_alignment_report",
        "logs/02.rmhost_report",
        "logs/02.rmhost_report_refine",
        "logs/02.rmhost_report_merge",
        "logs/04.profiling_custom_bgi_soap",
        "logs/04.profiling_custom_bgi_soap_merge",
        "logs/04.profiling_custom_bowtie2",
        "logs/04.profiling_custom_bowtie2_merge",
        "logs/04.profiling_custom_jgi",
        "logs/04.profiling_custom_jgi_merge",
        "logs/04.profiling_kraken2",
        "logs/04.profiling_bracken",
        "logs/04.profiling_bracken_merge",
        "logs/04.profiling_kmcp_search",
        "logs/04.profiling_kmcp_search_merge",
        "logs/04.profiling_kmcp_profile",
        "logs/04.profiling_kmcp_profile_merge",
        "logs/04.profiling_metaphlan2",
        "logs/04.profiling_metaphlan2_merge",
        "logs/04.profiling_metaphlan3",
        "logs/04.profiling_metaphlan3_merge",
        "logs/04.profiling_metaphlan4",
        "logs/04.profiling_metaphlan4_merge",
        "logs/04.profiling_strainphlan3_sample2markers",
        "logs/04.profiling_strainphlan3_extract_markers",
        "logs/04.profiling_strainphlan3",
        "logs/04.profiling_strainphlan4_sample2markers",
        "logs/04.profiling_strainphlan4_extract_markers",
        "logs/04.profiling_strainphlan4",
        "logs/04.profiling_humann2_config",
        "logs/04.profiling_humann2_build_chocophlan_pangenome_db",
        "logs/04.profiling_humann2",
        "logs/04.profiling_humann2_postprocess",
        "logs/04.profiling_humann2_join",
        "logs/04.profiling_humann2_split_stratified",
        "logs/04.profiling_humann3_config",
        "logs/04.profiling_humann3_build_chocophlan_pangenome_db",
        "logs/04.profiling_humann3",
        "logs/04.profiling_humann3_postprocess",
        "logs/04.profiling_humann3_join",
        "logs/04.profiling_humann3_split_stratified",
        "logs/04.profiling_alignment_bowtie2",
        "logs/04.profiling_alignment_bam_postprocess",
        "logs/04.profiling_genomecov_gen_bed",
        "logs/04.profiling_genomecov_gen_cov",
        "logs/04.profiling_genomecov_gen_cov_merge",
        "logs/04.profiling_genome_coverm",
        "logs/04.profiling_genome_coverm_merge"
    ]

    def __init__(self, work_dir):
        self.work_dir = os.path.realpath(work_dir)
        self.quantpi_dir = os.path.dirname(os.path.abspath(__file__))

        self.config_file = os.path.join(self.quantpi_dir, "config", "config.yaml")
        self.envs_dir = os.path.join(self.quantpi_dir, "envs")
        self.profiles_dir = os.path.join(self.quantpi_dir, "profiles")
        self.new_config_file = os.path.join(self.work_dir, "config.yaml")

    def __str__(self):
        message = """

     ██████  ██    ██  █████  ███    ██ ████████ ██████  ██ 
    ██    ██ ██    ██ ██   ██ ████   ██    ██    ██   ██ ██ 
    ██    ██ ██    ██ ███████ ██ ██  ██    ██    ██████  ██ 
    ██ ▄▄ ██ ██    ██ ██   ██ ██  ██ ██    ██    ██      ██ 
     ██████   ██████  ██   ██ ██   ████    ██    ██      ██ 
        ▀▀                                                  

           Omics for All, Open Source for All

A general profiling system focus on robust microbiome research.

Thanks for using quantpi. 

A metagenomics project has been created at %s


if you want to create fresh conda environments:

        quantpi --conda-create-envs-only
        quantpi --conda-create-envs-only

Looking for help:

        quantpi --help
""" % (
            self.work_dir
        )

        return message

    def create_dirs(self):
        """
        create project directory
        """
        if not os.path.exists(self.work_dir):
            os.mkdir(self.work_dir)

        for sub_dir in metaconfig.sub_dirs:
            os.makedirs(os.path.join(self.work_dir, sub_dir), exist_ok=True)

        for i in os.listdir(self.envs_dir):
            dest_file = os.path.join(self.work_dir, "envs", i)
            if os.path.exists(dest_file):
                print(f"{dest_file} exists, please remove or backup it first")
                sys.exit(-1)
            else:
                shutil.copyfile(os.path.join(self.envs_dir, i), dest_file)

        for i in os.listdir(self.profiles_dir):
            dest_dir = os.path.join(self.work_dir, "profiles", i)
            if os.path.exists(dest_dir):
                print(f"{dest_dir} exists, please remove or backup it first")
                sys.exit(-1)
            else:
                shutil.copytree(os.path.join(self.profiles_dir, i), dest_dir)

    def get_config(self):
        """
        get default configuration
        """
        return parse_yaml(self.config_file)


# https://github.com/Ecogenomics/CheckM/blob/master/checkm/customHelpFormatter.py
class custom_help_formatter(argparse.HelpFormatter):
    """Provide a customized format for help output.
    http://stackoverflow.com/questions/9642692/argparse-help-without-duplicate-allcaps
    """

    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if "%(default)" not in action.help:
            if (
                action.default != ""
                and action.default != []
                and action.default != None
                and action.default != False
            ):
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [
                        argparse.OPTIONAL, argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:
                        if "\n" in h:
                            lines = h.splitlines()
                            lines[0] += " (default: %(default)s)"
                            h = "\n".join(lines)
                        else:
                            h += " (default: %(default)s)"
            return h

    def _fill_text(self, text, width, indent):
        return "".join([indent + line for line in text.splitlines(True)])

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            (metavar,) = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(option_string)

                return "%s %s" % (", ".join(parts), args_string)

            return ", ".join(parts)

    def _get_default_metavar_for_optional(self, action):
        return action.dest.upper()

    def _get_default_metavar_for_positional(self, action):
        return action.dest