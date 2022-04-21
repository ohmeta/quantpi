#!/usr/bin/env python

import glob
import os
import sys

import pandas as pd


def parse_samples(
    samples_tsv,
    interleaved="false",
    reads_layout="pe",
    begin_point="trimming",
    check_samples=False,
):
    samples_df = pd.read_csv(
        samples_tsv, sep="\t", dtype={
            "sample_id": str,
            "fq1": str, "fq2": str})
    samples_df = samples_df.set_index(["sample_id"])

    cancel = False
    if "fq1" in samples_df.columns:
        for sample_id in samples_df.index.unique():
            sample_id = str(sample_id)
            if "." in sample_id:
                print(f"{sample_id} contain '.', please remove '.', now quiting :)")
                cancel = True

            fq1_list = samples_df.loc[[sample_id]]["fq1"].dropna().tolist()
            fq2_list = samples_df.loc[[sample_id]]["fq2"].dropna().tolist()
            for fq_file in fq1_list:
                if not fq_file.endswith(".gz"):
                    print(f"{fq_file} need gzip format")
                    cancel = True
                if check_samples:
                    if not os.path.exists(fq_file):
                        print(f"{fq_file} not exists")
                        cancel = True
                    if (reads_layout == "pe") and (not interleaved):
                        if len(fq2_list) == 0:
                            print(f"{sample_id} fq2 not exists")
                            cancel = True
    elif "sra" in samples_df.columns:
        for sample_id in samples_df.index.unique():
            sample_id = str(sample_id)
            if "." in sample_id:
                print(f"{sample_id} contain '.', please remove '.', now quiting :)")
                cancel = True

            if check_samples:
                sra_list = samples_df.loc[[sample_id]]["sra"].dropna().tolist()
                for sra_file in sra_list:
                    if not os.path.exists(sra_file):
                        print(f"{sra_file} not exists")
                        cancel = True
    else:
        print("wrong header: {header}".format(header=samples_df.columns))
        cancel = True

    if cancel:
        sys.exit(-1)
    else:
        return samples_df


def get_reads(sample_df, wildcards, col):
    return sample_df.loc[[wildcards.sample]][col].dropna().tolist()