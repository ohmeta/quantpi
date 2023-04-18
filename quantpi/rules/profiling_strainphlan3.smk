## reference
## https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-3


if config["params"]["profiling"]["strainphlan"]["do_v3"]:
    rule profiling_strainphlan3_sample2markers:
        input:
            sam = os.path.join(config["output"]["profiling"],
                               "profile/metaphlan3/{sample}/{sample}.sam.bz2"),
            aln = os.path.join(config["output"]["profiling"],
                               "profile/metaphlan3/{sample}/{sample}.bowtie2.bz2")
        output:
            os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/consensus_markers/{sample}.pkl")
        log:
            os.path.join(
                config["output"]["profiling"],
                "logs/strainphlan3_sample2markers/{sample}.strainphlan3_sample2markers.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"],
                "benchmark/strainphlan3_sample2markers/{sample}.strainphlan3_sample2markers.benchmark.txt")
        params:
            outdir =  os.path.join(config["output"]["profiling"], "profile/strainphlan3/consensus_markers")
        conda:
            config["envs"]["biobakery3"]
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            sample2markers.py \
            --input {input.sam} \
            --output_dir {params.outdir}/ \
            --nprocs {threads} \
            > {log} 2>&1
            '''


    if config["params"]["profiling"]["strainphlan"]["reference_genome"]["use"]:

        STRAINPHLAN_CLADES_V3 = \
            pd.read_csv(config["params"]["profiling"]["strainphlan"]["clades_tsv_v3"], sep="\t")\
            .set_index("clade")
        STRAINPHLAN_CLADES_LIST_V3 = STRAINPHLAN_CLADES_V3.index.unique()


        rule profiling_strainphlan3_prepare_reference_genome:
            input:
                reference_genome = lambda wildcards: STRAINPHLAN_CLADES_V3.loc[wildcards.clade, "fna_path"]
            output:
                reference_genome = os.path.join(
                    config["output"]["profiling"],
                    "databases/strainphlan3/reference_genomes/{clade}.fna")
            run:
                if input.reference_genome.endswith(".gz"):
                    shell(f'''pigz -fkdc {input.reference_genome} > {output.reference_genome}''')
                else:
                    fna_path = os.path.realpath(input.reference_genome)
                    fna_path_ = os.path.realpath(output.reference_genome)
                    shellf('''ln -s {fna_path} {fna_path_}''')

    else:
        checkpoint profiling_strainphlan3_print_clades:
            input:
                database_pkl = expand(os.path.join(
                    config["params"]["profiling"]["metaphlan"]["bowtie2db"], "{index}.pkl"),
                    index = config["params"]["profiling"]["metaphlan"]["index_v3"]),
                consensus_markers = expand(os.path.join(
                    config["output"]["profiling"],
                    "profile/strainphlan3/consensus_markers/{sample}.pkl"),
                    sample=SAMPLES_ID_LIST)
            output:
                markers_txt = os.path.join(config["output"]["profiling"], "databases/strainphlan3/clade_markers.txt"),
                markers_dir = os.path.join(config["output"]["profiling"], "databases/strainphlan3/clade_markers")
            params:
                marker_in_n_samples = config["params"]["profiling"]["strainphlan"]["marker_in_n_samples"]
            conda:
                config["envs"]["biobakery3"]
            shell:
                '''
                strainphlan \
                --database {input.database_pkl} \
                --samples {input.consensus_markers} \
                --marker_in_n_samples {params.marker_in_n_samples} \
                -o /dev/null \
                --print_clades_only \
                >{output.markers_txt} 2>{log}

                mkdir -p {output.markers_dir}

                cat {output.markers_txt} | \
                grep "s__" | \
                awk -F'[\t:]' '{{print $5}}' | \
                xargs -I XXX mkdir -p {output.markers_dir}/XXX
                '''


        def aggregate_profiling_strainphlan3_extract_markers(wildcards):
            checkpoint_output = checkpoints.profiling_strainphlan3_print_clades.get(**wildcards).output.markers_dir

            return expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/clade_markers/{clade}/done"),
                clade=list(set([i.split("/")[0] \
                    for i in glob_wildcards(os.path.join(
                        checkpoint_output,
                        "{{clade}}",
                        f"{{clade}}.fna")).clade])))


    rule profiling_strainphlan3_extract_markers:
        input:
            database_pkl = expand(os.path.join(
                config["params"]["profiling"]["metaphlan"]["bowtie2db"], "{index}.pkl"),
                index = config["params"]["profiling"]["metaphlan"]["index_v3"])
        output:
            clade_marker = os.path.join(
                config["output"]["profiling"],
                "databases/strainphlan3/clade_markers/{clade}/{clade}.fna")
        log:
            os.path.join(
                config["output"]["profiling"], 
                "logs/strainphlan3_extract_markers/{clade}.strainphlan3_extract_markers.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"],
                "benchmark/strainphlan3_extract_markers/{clade}.strainphlan3_extract_markers.benchmark.txt")
        conda:
            config["envs"]["biobakery3"]
        params:
            clade = "{clade}",
            outdir = os.path.join(config["output"]["profiling"], "databases/strainphlan3/clade_markers/{clade}")
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            mkdir -p {params.outdir}

            extract_markers.py \
            --database {input.database_pkl} \
            --clade {params.clade} \
            --output_dir {params.outdir}/ \
            > {log} 2>&1
            '''


    rule profiling_strainphlan3:
        input:
            database_pkl = expand(os.path.join(
                config["params"]["profiling"]["metaphlan"]["bowtie2db"], "{index}.pkl"),
                index = config["params"]["profiling"]["metaphlan"]["index_v3"]),
            clade_marker = os.path.join(
                config["output"]["profiling"],
                "databases/strainphlan3/clade_markers/{clade}/{clade}.fna"),
            consensus_markers = expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/consensus_markers/{sample}.pkl"),
                sample=SAMPLES_ID_LIST)
        output:
            done = os.path.join(config["output"]["profiling"], "profile/strainphlan3/clade_markers/{clade}/done")
        log:
            os.path.join(
                config["output"]["profiling"],
                "logs/strainphlan3/{clade}.strainphlan3.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"],
                "benchmark/strainphlan3/{clade}.strainphlan3.benchmark.txt")
        conda:
            config["envs"]["biobakery3"]
        params:
            clade = "{clade}",
            outdir = os.path.join(config["output"]["profiling"], "profile/strainphlan3/clade_markers/{clade}"),
            trim_sequences = config["params"]["profiling"]["strainphlan"]["trim_sequences"],
            marker_in_n_samples = config["params"]["profiling"]["strainphlan"]["marker_in_n_samples"],
            sample_with_n_markers = config["params"]["profiling"]["strainphlan"]["sample_with_n_markers"],
            secondary_sample_with_n_markers = config["params"]["profiling"]["strainphlan"]["secondary_sample_with_n_markers"],
            #sample_with_n_markers_after_filt = config["params"]["profiling"]["strainphlan"]["sample_with_n_markers_after_filt"],
            phylophlan_mode = config["params"]["profiling"]["strainphlan"]["phylophlan_mode"],
            reference_opts = "--references %s" % \
            os.path.join(config["output"]["profiling"], "databases/strainphlan3/reference_genomes/{clade}.fna") \
            if config["params"]["profiling"]["strainphlan"]["reference_genome"]["use"] \
            else "",
            opts = config["params"]["profiling"]["strainphlan"]["external_opts_v3"]
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            rm -rf {params.outdir}
            mkdir -p {params.outdir}

            strainphlan \
            --database {input.database_pkl} \
            --samples {input.consensus_markers} \
            --clade_markers {input.clade_marker} \
            {params.reference_opts} \
            --output_dir {params.outdir}/ \
            --nprocs {threads} \
            --clade {params.clade} \
            --trim_sequences {params.trim_sequences} \
            --marker_in_n_samples {params.marker_in_n_samples} \
            --sample_with_n_markers {params.sample_with_n_markers} \
            --secondary_sample_with_n_markers {params.secondary_sample_with_n_markers} \
            --phylophlan_mode {params.phylophlan_mode} \
            --mutation_rates \
            {params.opts} \
            >{log} 2>&1

            touch {output.done}
            '''


    if config["params"]["profiling"]["strainphlan"]["reference_genome"]["use"]:
        rule profiling_strainphlan3_all:
            input:
                expand(os.path.join(
                    config["output"]["profiling"],
                    "profile/strainphlan3/clade_markers/{clade}/done"),
                    clade=STRAINPHLAN_CLADES_LIST_V3)

    else:
        rule profiling_strainphlan3_all:
            input:
                aggregate_profiling_strainphlan3_extract_markers


else:
    rule profiling_strainphlan3_all:
        input: