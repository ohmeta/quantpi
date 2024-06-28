#!/usr/bin/env python

import pandas as pd
import concurrent.futures
import os
import sys
import argparse
import gzip
#import pysam
import re
import numpy as np
from natsort import index_natsorted
from pprint import pprint


def kmcp_gen_species_tax_profile(df, target, sample_name):
    df_ = df.loc[:, ["taxpath", target]]
    if not df_.empty:
        df_["taxa_name"] = df_.apply(lambda x: ";".join(x["taxpath"].split(";")[0:7]), axis=1)
        df_ = df_.loc[:, ["taxa_name", target]].groupby("taxa_name").sum(target).rename(columns={target: sample_name})
    else:
        df_ = df_.rename(columns={"taxpath": "taxa_name", target: sample_name}).set_index("taxa_name")
    return df_


def kmcp_profile_parse_species(profile_file):
    sample_name = os.path.basename(profile_file).split(".")[0]

    df = pd.read_csv(profile_file, sep="\t")

    df_percentage_s = kmcp_gen_species_tax_profile(df, "percentage", sample_name)
    df_coverage_s = kmcp_gen_species_tax_profile(df, "coverage", sample_name)
    df_reads_s = kmcp_gen_species_tax_profile(df, "reads", sample_name)
    df_ureads_s = kmcp_gen_species_tax_profile(df, "ureads", sample_name)
    df_hicureads_s = kmcp_gen_species_tax_profile(df, "hicureads", sample_name)

    return { "percentage_s": df_percentage_s,
             "coverage_s": df_coverage_s,
             "reads_s": df_reads_s,
             "ureads_s": df_ureads_s,
             "hicureads_s": df_hicureads_s
             }


def kmcp_profile_parse(profile_file):
    sample_name = os.path.basename(profile_file).split(".")[0]
    df = pd.read_csv(profile_file, sep="\t")

    df_percentage_t = df.query('rank=="strain"').loc[:, ["taxpath", "percentage"]]\
                        .set_index("taxpath").rename(columns={"percentage": sample_name})
    df_percentage_s = df.query('rank=="species"').loc[:, ["taxpath", "percentage"]]\
                        .set_index("taxpath").rename(columns={"percentage": sample_name})

    df_coverage_t = df.query('rank=="strain"').loc[:, ["taxpath", "coverage"]]\
                      .set_index("taxpath").rename(columns={"coverage": sample_name})
    df_coverage_s = df.query('rank=="species"').loc[:, ["taxpath", "coverage"]]\
                      .set_index("taxpath").rename(columns={"coverage": sample_name})

    df_reads_t = df.query('rank=="strain"').loc[:, ["taxpath", "reads"]]\
                   .set_index("taxpath").rename(columns={"reads": sample_name})
    df_reads_s = df.query('rank=="species"').loc[:, ["taxpath", "reads"]]\
                   .set_index("taxpath").rename(columns={"reads": sample_name})

    df_ureads_t = df.query('rank=="strain"').loc[:, ["taxpath", "ureads"]]\
                    .set_index("taxpath").rename(columns={"ureads": sample_name})
    df_ureads_s = df.query('rank=="species"').loc[:, ["taxpath", "ureads"]]\
                    .set_index("taxpath").rename(columns={"ureads": sample_name})

    df_hicureads_t = df.query('rank=="strain"').loc[:, ["taxpath", "hicureads"]]\
                       .set_index("taxpath").rename(columns={"hicureads": sample_name})
    df_hicureads_s = df.query('rank=="species"').loc[:, ["taxpath", "hicureads"]]\
                       .set_index("taxpath").rename(columns={"hicureads": sample_name})

    return { "percentage_t": df_percentage_t,
             "percentage_s": df_percentage_s,
             "coverage_t": df_coverage_t,
             "coverage_s": df_coverage_s,
             "reads_t": df_reads_t,
             "reads_s": df_reads_s,
             "ureads_t": df_ureads_t,
             "ureads_s": df_ureads_s,
             "hicureads_t": df_hicureads_t,
             "hicureads_s": df_hicureads_s
             }


def kmcp_profile_merge_species(profile_file_list, workers, output_dict):
    df_dict = {
        "percentage_s": [],
        "coverage_s": [],
        "reads_s": [],
        "ureads_s": [],
        "hicureads_s": []
    }

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for dfs in executor.map(kmcp_profile_parse_species, profile_file_list):
            for target in df_dict.keys():
                df_dict[target].append(dfs[target])

    for target in df_dict:
        df = pd.concat(df_dict[target], axis=1).fillna(0).reset_index()
        df = df.sort_values(by="taxa_name", key=lambda x: np.argsort(index_natsorted(df["taxa_name"])))
        df.to_csv(output_dict[target], sep="\t", index=False)


def kmcp_profile_merge(profile_file_list, workers, output_dict):
    df_dict = {
        "percentage_t": [],
        "percentage_s": [],
        "coverage_t": [],
        "coverage_s": [],
        "reads_t": [],
        "reads_s": [],
        "ureads_t": [],
        "ureads_s": [],
        "hicureads_t": [],
        "hicureads_s": []
    }

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for dfs in executor.map(kmcp_profile_parse, profile_file_list):
            for target in df_dict.keys():
                df_dict[target].append(dfs[target])

    for target in df_dict:
        df = pd.concat(df_dict[target], axis=1).fillna(0).reset_index()
        df = df.sort_values(by="taxpath", key=lambda x: np.argsort(index_natsorted(df["taxpath"])))
        df.to_csv(output_dict[target], sep="\t", index=False)


def coverm_parse(table_file):
    sample_name = os.path.basename(table_file).split(".")[0]
    df = pd.read_csv(table_file, sep="\t")
    names = df.columns.tolist()
    df_list = []
    for method_name in names[1:]:
        df_ = df.loc[:, [names[0], method_name]].rename(columns={method_name: sample_name}).set_index(names[0])
        df_list.append(df_)
    return df_list


def coverm_merge(table_files, method_list, workers, output_list):
    df_dict = {}
    for method_name in method_list:
        df_dict[method_name] = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df_list in executor.map(coverm_parse, table_files):
            for i in range(len(df_list)):
                df_dict[method_list[i]].append(df_list[i])

    for i in range(len(method_list)):
        df_out = pd.concat(df_dict[method_list[i]], axis=1).reset_index()
        df_out = df_out.sort_values(by="Genome",
                                    key=lambda x: np.argsort(index_natsorted(df_out["Genome"])))
        df_out.to_csv(output_list[i], sep="\t", index=False)


def genomecov_parse(tsv_file):
    sample_name = os.path.basename(tsv_file).split(".")[0]
    df = pd.read_csv(tsv_file, sep="\t")
    df = df.rename(columns={
        "cov_mean_sample_0": f"cov_mean_{sample_name}",
        "percentage_covered_sample_0": f"percentage_covered_{sample_name}"})
    df = df.set_index(["contig", "length", "GC"])
    #print(df.head())

    return df.loc[:, [f"cov_mean_{sample_name}"]], df.loc[:, [f"percentage_covered_{sample_name}"]]


def genomecov_merge(tsv_files, workers, **kwargs):
    df_list_cov = []
    df_list_per = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df_cov, df_per in executor.map(genomecov_parse, tsv_files):
            df_list_cov.append(df_cov)
            df_list_per.append(df_per)

    cov_df = pd.concat(df_list_cov, axis=1).reset_index()
    per_df = pd.concat(df_list_per, axis=1).reset_index()

    if "output_cov" in kwargs:
        cov_df.to_csv(kwargs["output_cov"], sep="\t", index=False)
    if "output_per" in kwargs:
        per_df.to_csv(kwargs["output_per"], sep="\t", index=False)


def metaphlan_init(version):
    global METAPHLAN_VERSION
    METAPHLAN_VERSION = version


def read_metaphlan_table(table):
    sample_id = os.path.basename(table).split(".")[0]

    if METAPHLAN_VERSION == 2:
        dict_ = {"clade_name": [], sample_id: []}
        with open(table, "r") as ih:
            for line in ih:
                if not line.startswith("#"):
                    clade_name, abun = line.strip().split("\t")[0:2]
                    dict_["clade_name"].append(clade_name)
                    dict_[sample_id].append(abun)
        df = pd.DataFrame(dict_).set_index("clade_name")
        return df

    elif (METAPHLAN_VERSION == 3) or (METAPHLAN_VERSION == 4):
        dict_ = {"clade_name": [], "clade_taxid": [], sample_id: []}
        with open(table, "r") as ih:
            for line in ih:
                if not line.startswith("#"):
                    clade_name, clade_taxid, abun = line.strip().split("\t")[
                        0:3]
                    dict_["clade_name"].append(clade_name)
                    dict_["clade_taxid"].append(clade_taxid)
                    dict_[sample_id].append(abun)
        df = pd.DataFrame(dict_).set_index(["clade_name", "clade_taxid"])
        return df
    else:
        return None


def merge_metaphlan_tables(table_files, workers, **kwargs):
    abun_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for abun_df in executor.map(read_metaphlan_table, table_files):
            if abun_df is not None:
                abun_list.append(abun_df)

    abun_df_ = pd.concat(abun_list, axis=1).fillna(0)\
                 .reset_index().set_index("clade_name")
    if "output" in kwargs:
        abun_df_.reset_index().to_csv(kwargs["output"], sep="\t", index=False)

    df_list = []
    for i in ["t", "s", "g", "f", "o", "c", "p", "k"]:
        profile_df = pd.DataFrame()
        if METAPHLAN_VERSION == 2:
            regex_pattern = rf"\|{i}__[^\|]*$"
            if i == "k":
                regex_pattern = rf"k__[^\|]*$"
            profile_df = abun_df_.filter(regex=regex_pattern, axis=0).reset_index() 

        elif METAPHLAN_VERSION == 3:
            regex_pattern = rf"UNKNOWN|\|{i}__[^\|]*$"
            if i == "k":
                regex_pattern = rf"UNKNOWN|k__[^\|]*$"
            profile_df = abun_df_.filter(regex=regex_pattern, axis=0).reset_index()

        elif METAPHLAN_VERSION == 4:
            regex_pattern = rf"UNCLASSIFIED|\|{i}__[^\|]*$"
            if i == "k":
                regex_pattern = rf"UNCLASSIFIED|k__[^\|]*$"
            profile_df = abun_df_.filter(regex=regex_pattern, axis=0).reset_index()

        df_list.append(profile_df)
    return [abun_df_.reset_index()] + df_list


def profiler_init(index_metadata):
    global INDEX_METADATA__
    INDEX_METADATA__ = pd.read_csv(index_metadata, sep="\t")


def set_lineages_to(row, key, level):
    LINEAGES = [
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain",
    ]
    LINEAGES_DICT = {
        "strain": LINEAGES,
        "species": LINEAGES[0:7],
        "genus": LINEAGES[0:6],
        "family": LINEAGES[0:5],
        "order": LINEAGES[0:4],
        "class": LINEAGES[0:3],
        "phylum": LINEAGES[0:2],
        "superkingdom": LINEAGES[0:1],
    }

    LEVEL_DICT = {
        "strain": "t",
        "species": "s",
        "genus": "g",
        "family": "f",
        "order": "o",
        "class": "c",
        "phylum": "p",
        "superkingdom": "k",
    }

    lineages_dict = {
        "k": "k__unclassified_" + row["mgs_id"],
        "p": "p__unclassified_" + row["mgs_id"],
        "c": "c__unclassified_" + row["mgs_id"],
        "o": "o__unclassified_" + row["mgs_id"],
        "f": "f__unclassified_" + row["mgs_id"],
        "g": "g__unclassified_" + row["mgs_id"],
        "s": "s__unclassified_" + row["mgs_id"],
        "t": "t__unclassified_" + row["mgs_id"],
    }
    for line in row[key].split(";"):
        lev, tax = line.split("__")
        if lev == "d":
            lev = "k"
        if tax != "":
            lineages_dict[lev] = lev + "__" + tax
            if lev == "s" or lev == "t":
                lineages_dict[lev] = lineages_dict[lev] + "_" + row["mgs_id"]

    lineages = []
    for i in LINEAGES_DICT[level]:
        lineages.append(lineages_dict[LEVEL_DICT[i]])

    return "|".join(lineages)


def get_mgs_id(row):
    return "_".join(row["ID"].split("_")[0:-1])


def get_abun_df_bgi_soap(soap_file):
    sample_id = os.path.basename(soap_file).split(".")[0]

    reads_count_dict = {}
    with gzip.open(soap_file, 'rt') as h:
        for line in h:
            ref_name = line.split("\t")[7]
            if ref_name in reads_count_dict:
                reads_count_dict[ref_name] += 1
            else:
                reads_count_dict[ref_name] = 1
    reads_count_df = pd.DataFrame(list(reads_count_dict.items()), columns=[
                                  "reference_name", "reads_count"])
    abun = reads_count_df.merge(INDEX_METADATA__)

    '''
    abun["count_by_len"] = abun["reads_count"] / abun["reference_length"]
    abun["count_by_len_rate"] = abun["count_by_len"] / \
        sum(abun["count_by_len"])

    # geneset method 1
    abun_df = abun.groupby("lineages_full")\
                  .agg({"count_by_len_rate": "sum"})\
                  .rename(columns={"count_by_len_rate": sample_id})
    # geneset method 2
    count_by_len_df = abun.groupby("lineages_full").agg(
        {"count_by_len": "sum"})
    count_by_len_df["count_by_len_rate"] = count_by_len_df["count_by_len"] / \
        sum(count_by_len_df["count_by_len"])
    abun_df = count_by_len_df.loc[:, ["count_by_len_rate"]]\
                             .rename(columns={"count_by_len_rate": sample_id})
    '''

    # count method
    abun_count = abun.groupby("lineages_full").agg(
        {"reads_count": "sum"})
    abun_count["count_rate"] = abun_count["reads_count"] / \
        sum(abun_count["reads_count"])

    abun_df = abun_count.loc[:, ["count_rate"]].rename(
        columns={"count_rate": sample_id})
    count_df = abun_count.loc[:, ["reads_count"]].rename(
        columns={"reads_count": sample_id})
    return count_df, abun_df


#def get_abun_df_bowtie2(bam):
#    sample_id = os.path.basename(bam).split(".")[0]
#
#    reads_count_dict = {}
#    sam = pysam.AlignmentFile(bam, "rb")
#    for record in sam:
#        if not record.is_unmapped:
#            tmp_nm = record.get_tag("NM:i")
#            tmp_len = sum([int(i) for i in re.findall(
#                r"(\d+)(?:M|I|D)", record.cigarstring)])
#
#            if ((1 - tmp_nm / tmp_len)) >= 0.95:
#                if record.reference_name in reads_count_dict:
#                    reads_count_dict[record.reference_name] += 1
#                else:
#                    reads_count_dict[record.reference_name] = 1
#
#    reads_count_df = pd.DataFrame(list(reads_count_dict.items()), columns=[
#                                  "reference_name", "reads_count"])
#    abun = reads_count_df.merge(INDEX_METADATA__)
#
#    abun_count = abun.groupby("lineages_full").agg(
#        {"reads_count": "sum"})
#    abun_count["count_rate"] = abun_count["reads_count"] / \
#        sum(abun_count["reads_count"])
#
#    abun_df = abun_count.loc[:, ["count_rate"]].rename(
#        columns={"count_rate": sample_id})
#    count_df = abun_count.loc[:, ["reads_count"]].rename(
#        columns={"reads_count": sample_id})
#    return count_df, abun_df


def get_abun_df_hsx(abun_file):
    sample_id = os.path.basename(abun_file).split(".")[0]

    try:
        if os.path.exists(abun_file):
            abun = pd.read_csv(abun_file, sep="\t")
        else:
            print("%s is not exists" % abun_file)
            return None, None
    except pd.errors.EmptyDataError:
        print("%s is empty" % abun_file)
        return None, None

    abun["mgs_id"] = abun.apply(get_mgs_id, axis=1)

    count_df = (
        abun.loc[:, ["mgs_id", "reads_pairs"]]
        .groupby("mgs_id")
        .agg({"reads_pairs": "sum"})
        .rename(columns={"reads_pairs": sample_id})
    )
    abun_df = (
        abun.loc[:, ["mgs_id", "gene_abundance"]]
        .groupby("mgs_id")
        .agg({"gene_abundance": "sum"})
        .rename(columns={"gene_abundance": sample_id})
    )
    return count_df, abun_df


def get_abun_df_jgi(depth_file):
    sample_id = os.path.basename(depth_file).split(".")[0]

    try:
        if os.path.exists(depth_file):
            depth = pd.read_csv(depth_file, sep="\t")
        else:
            print("%s is not exists" % depth_file)
            return None, None
    except pd.errors.EmptyDataError:
        print("%s is empty" % depth_file)
        return None, None

    depth = (
        depth.rename(columns={"contigName": "contig_name"})
        .merge(INDEX_METADATA__)
        .groupby("mgs_id")
        .agg({"totalAvgDepth": "mean"})
    )
    depth[sample_id] = depth["totalAvgDepth"] / sum(depth["totalAvgDepth"])
    depth_df = depth.loc[:, ["totalAvgDepth"]].rename(
        columns={"totalAvgDepth": sample_id}
    )
    abun_df = depth.loc[:, [sample_id]]
    return depth_df, abun_df


def get_all_abun_df(abun_files, workers, method):
    count_list = []
    abun_list = []
    if method == "jgi":
        func = get_abun_df_jgi
    elif method == "hsx":
        func = get_abun_df_hsx
    elif method == "bgi_soap":
        func = get_abun_df_bgi_soap
    elif method == "bowtie2":
        func = get_abun_df_bowtie2
    else:
        print("unspoort method %s" % method)

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for count_df, abun_df in executor.map(func, abun_files):
            if (count_df is not None) and (abun_df is not None):
                count_list.append(count_df)
                abun_list.append(abun_df)

    count_df_ = pd.concat(count_list, axis=1).reset_index()
    abun_df_ = pd.concat(abun_list, axis=1).reset_index()

    return count_df_, abun_df_


def get_profile(abun_tax_df, samples_list, key, profile_tsv):
    # level_ = "lineages_" + level + "_new"
    # abun_tax_df[level_] = abun_tax_df.apply(lambda x: set_lineages_to(x, level), axis=1)
    level_ = key
    _profile = (
        abun_tax_df.loc[:, [level_] + samples_list]
        .melt(id_vars=[level_])
        .groupby(["variable", level_])
        .agg({"value": "sum"})
        .reset_index()
        .pivot(index=level_, columns="variable", values="value")
    )
    profile_ = _profile.reset_index().loc[:, [level_] + samples_list]
    profile_.to_csv(profile_tsv, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser("metagenomics species abundance profiler")
    parser.add_argument("-l", "--abundance_list",
                        type=str, help="abundance list")
    parser.add_argument(
        "--method", default="hsx", choices=["hsx", "jgi"], help="compute method"
    )
    parser.add_argument(
        "--database", default=None, help="contig and genome relationships"
    )
    parser.add_argument(
        "--taxonomy", default=None, help="genome database taxonomy information"
    )
    parser.add_argument("--threads", default=8, type=int, help="threads")
    parser.add_argument("--out_prefix", default="./",
                        type=str, help="output prefix")

    args = parser.parse_args()

    abun_files = pd.read_csv(args.abundance_list, names=[
                             "path"]).loc[:, "path"].values

    if args.method == "jgi" and args.database is None:
        print("pleas supply database when parse jgi depth file")
        sys.exit(1)

    if args.method == "hsx":
        count_df, abun_df = get_all_abun_df(abun_files, args.threads, "hsx")
    elif args.method == "jgi":
        profiler_init(args.database)
        count_df, abun_df = get_all_abun_df(abun_files, args.threads, "jgi")
    else:
        print("unsupport method: %s" % args.method)

    outdir = os.path.dirname(args.out_prefix)
    if (outdir == "") or (outdir == "."):
        outdir = "."
    else:
        os.makedirs(outdir, exist_ok=True)
    outprefix = os.path.basename(args.out_prefix)
    count_profile = os.path.join(outdir, outprefix + "_count_profile.tsv")
    abun_profile = os.path.join(outdir, outprefix + "_abundance_profile.tsv")
    count_df.to_csv(count_profile, sep="\t", index=False)
    abun_df.to_csv(abun_profile, sep="\t", index=False)

    if args.taxonomy is not None:
        taxonomy_df = pd.read_csv(args.taxonomy, sep="\t")
        abun_tax_df = abun_df.merge(taxonomy_df)
        samples_list = sorted(abun_df.columns[1:].to_list())

        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_superkingdom_new",
            os.path.join(outdir, outprefix +
                         "_abundance_profile_superkingdom.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_phylum_new",
            os.path.join(outdir, outprefix + "_abundance_profile_phylum.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_order_new",
            os.path.join(outdir, outprefix + "_abundance_profile_order.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_class_new",
            os.path.join(outdir, outprefix + "_abundance_profile_class.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_family_new",
            os.path.join(outdir, outprefix + "_abundance_profile_family.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_genus_new",
            os.path.join(outdir, outprefix + "_abundance_profile_genus.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_species_new",
            os.path.join(outdir, outprefix + "_abundance_profile_species.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_strain_new",
            os.path.join(outdir, outprefix + "_abundance_profile_strain.tsv"),
        )


if __name__ == "__main__":
    main()
