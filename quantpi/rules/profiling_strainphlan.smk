rule profiling_strainphlan3_sample2markers:
    input:
        sam = os.path.join(config["output"]["profiling"], "profile/metaphlan3/{sample}/{sample}.sam.bz2"),
        aln = os.path.join(config["output"]["profiling"], "profile/metaphlan3/{sample}/{sample}.bowtie2.bz2")
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
            reference_genome = os.path.join(config["output"]["profiling"], "databases/strainphlan3/reference_genomes/{clade}.fna")
        threads:
            1
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
                config["params"]["profiling"]["metaphlan"]["bowtie2db_v3"], "{index}.pkl"),
                index = config["params"]["profiling"]["metaphlan"]["index_prefix_v3"]),
            consensus_markers = expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/consensus_markers/{sample}.pkl"),
                sample=SAMPLES_ID_LIST)
        output:
            markers_txt = os.path.join(config["output"]["profiling"], "databases/strainphlan3/clade_markers.txt"),
            markers_txt_gt = expand(
                os.path.join(config["output"]["profiling"], "databases/strainphlan3/clade_markers.txt.gt{num}"),
                num=config["params"]["profiling"]["strainphlan"]["clade_detected_in_min_samples_num"]),
            markers_dir = directory(os.path.join(config["output"]["profiling"], "databases/strainphlan3/clade_markers"))
        log:
            os.path.join(config["output"]["profiling"], "logs/strainphlan3_print_clades/strainphlan3_print_clades.log")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/strainphlan3_print_clades.txt")
        params:
            marker_in_n_samples = config["params"]["profiling"]["strainphlan"]["marker_in_n_samples"],
            clade_detected_in_min_samples_num = config["params"]["profiling"]["strainphlan"]["clade_detected_in_min_samples_num"]
        conda:
            config["envs"]["biobakery3"]
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            strainphlan \
            --nprocs {threads} \
            --database {input.database_pkl} \
            --samples {input.consensus_markers} \
            --marker_in_n_samples {params.marker_in_n_samples} \
            -o /dev/null \
            --print_clades_only \
            >{output.markers_txt} 2>{log}

            cat {output.markers_txt} | \
            grep "s__" | \
            awk -F'[\t]' '{{print $2}}' | \
            awk -F'[: ]' '{{print $1  "\t" $4}}' \
            awk -v samnum={params.clade_detected_in_min_samples_num} '$2>=samnum{{print $1}}' \
            > {output.markers_txt_gt}

            mkdir -p {output.markers_dir}

            cat {output.markers_txt_gt} | \
            xargs -I XXX mkdir -p {output.markers_dir}/XXX
            '''


    def aggregate_profiling_strainphlan3_extract_markers(wildcards):
        checkpoint_output = checkpoints.profiling_strainphlan3_print_clades.get(**wildcards).output.markers_dir

        clades = []
        for i in glob_wildcards(os.path.join(checkpoint_output, "{clade}")).clade:
            if "s__" in i:
                clades.append(i.split("/")[0])
        clades = list(set(clades))

        from pprint import pprint
        pprint(clades)

        clades_done = expand(os.path.join(
            config["output"]["profiling"],
            "profile/strainphlan3/clade_markers/{clade}/done"),
            clade = clades)

        return clades_done


rule profiling_strainphlan3_extract_markers:
    input:
        database_pkl = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v3"], "{index}.pkl"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v3"])
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
        1
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
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v3"], "{index}.pkl"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v3"]),
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
        marker_in_n_samples = config["params"]["profiling"]["strainphlan"]["marker_in_n_samples"],
        sample_with_n_markers = config["params"]["profiling"]["strainphlan"]["sample_with_n_markers"],
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
        --marker_in_n_samples {params.marker_in_n_samples} \
        --sample_with_n_markers {params.sample_with_n_markers} \
        --phylophlan_mode {params.phylophlan_mode} \
        --mutation_rates \
        {params.opts} \
        >{log} 2>&1

        touch {output.done}
        '''


if config["params"]["profiling"]["strainphlan"]["do_v3"]:
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



rule profiling_strainphlan40_sample2markers:
    input:
        database_pkl = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v40"], "{index}.pkl"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v40"]),
        sam = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan40/{sample}/{sample}.sam.bz2"),
        aln = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan40/{sample}/{sample}.bowtie2.bz2")
    output:
        os.path.join(
            config["output"]["profiling"],
            "profile/strainphlan40/consensus_markers/{sample}.pkl")
    log:
        os.path.join(
            config["output"]["profiling"],
            "logs/strainphlan40_sample2markers/{sample}.strainphlan40_sample2markers.log")
    benchmark:
        os.path.join(
            config["output"]["profiling"],
            "benchmark/strainphlan40_sample2markers/{sample}.strainphlan40_sample2markers.benchmark.txt")
    params:
        outdir =  os.path.join(config["output"]["profiling"], "profile/strainphlan40/consensus_markers")
    conda:
        config["envs"]["biobakery40"]
    priority:
        20
    threads:
        config["params"]["profiling"]["threads"]
    shell:
        '''
        sample2markers.py \
        --database {input.database_pkl} \
        --input {input.sam} \
        --output_dir {params.outdir}/ \
        --nprocs {threads} \
        > {log} 2>&1
        '''


if config["params"]["profiling"]["strainphlan"]["reference_genome"]["use"]:

    STRAINPHLAN_CLADES_V40 = \
        pd.read_csv(config["params"]["profiling"]["strainphlan"]["clades_tsv_v40"], sep="\t")\
        .set_index("clade")
    STRAINPHLAN_CLADES_LIST_V40 = STRAINPHLAN_CLADES_V40.index.unique()


    rule profiling_strainphlan40_prepare_reference_genome:
        input:
            reference_genome = lambda wildcards: STRAINPHLAN_CLADES_V40.loc[wildcards.clade, "fna_path"]
        output:
            reference_genome = os.path.join(
                config["output"]["profiling"],
                "databases/strainphlan40/reference_genomes/{clade}.fna")
        threads:
            1
        run:
            if input.reference_genome.endswith(".gz"):
                shell(f'''pigz -fkdc {input.reference_genome} > {output.reference_genome}''')
            else:
                fna_path = os.path.realpath(input.reference_genome)
                fna_path_ = os.path.realpath(output.reference_genome)
                shellf('''ln -s {fna_path} {fna_path_}''')

else:
    checkpoint profiling_strainphlan40_print_clades:
        input:
            database_pkl = expand(os.path.join(
                config["params"]["profiling"]["metaphlan"]["bowtie2db_v40"], "{index}.pkl"),
                index = config["params"]["profiling"]["metaphlan"]["index_prefix_v40"]),
            consensus_markers = expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan40/consensus_markers/{sample}.pkl"),
                sample=SAMPLES_ID_LIST)
        output:
            markers_txt = os.path.join(config["output"]["profiling"], "databases/strainphlan40/clade_markers.txt"),
            markers_txt_gt = expand(
                os.path.join(config["output"]["profiling"], "databases/strainphlan40/clade_markers.txt.gt{num}"),
                num=config["params"]["profiling"]["strainphlan"]["clade_detected_in_min_samples_num"]),
            markers_dir = directory(os.path.join(config["output"]["profiling"], "databases/strainphlan40/clade_markers"))
        log:
            os.path.join(config["output"]["profiling"], "logs/strainphlan40_print_clades/strainphlan40_print_clades.log")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/strainphlan40_print_clades.txt")
        params:
            marker_in_n_samples = config["params"]["profiling"]["strainphlan"]["marker_in_n_samples"],
            clade_detected_in_min_samples_num = config["params"]["profiling"]["strainphlan"]["clade_detected_in_min_samples_num"]
        conda:
            config["envs"]["biobakery40"]
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            OUTDIR=$(dirname {output.markers_txt})
            rm -rf $OUTDIR
            mkdir -p $OUTDIR

            strainphlan \
            --nprocs {threads} \
            --database {input.database_pkl} \
            --samples {input.consensus_markers} \
            --marker_in_n_samples {params.marker_in_n_samples} \
            --output_dir $OUTDIR \
            --print_clades_only \
            >{output.markers_txt} 2>{log}

            cat {output.markers_txt} | \
            grep "t__" | \
            awk -F'[\t]' '{{print $2}}' | \
            awk -F'[: ]' '{{print $1  "\t" $4}}' \
            awk -v samnum={params.clade_detected_in_min_samples_num} '$2>=samnum{{print $1}}' \
            > {output.markers_txt_gt}

            mkdir -p {output.markers_dir}

            cat {output.markers_txt_gt} | \
            xargs -I XXX mkdir -p {output.markers_dir}/XXX
            '''


    def aggregate_profiling_strainphlan40_extract_markers(wildcards):
        checkpoint_output = checkpoints.profiling_strainphlan40_print_clades.get(**wildcards).output.markers_dir

        clades = []
        for i in glob_wildcards(os.path.join(checkpoint_output, "{clade}")).clade:
            if "t__SGB" in i:
                clades.append(i.split("/")[0])
        clades = list(set(clades))

        from pprint import pprint
        pprint(clades)

        clades_done = expand(os.path.join(
            config["output"]["profiling"],
            "profile/strainphlan40/clade_markers/{clade}/done"),
            clade = clades)

        return clades_done


rule profiling_strainphlan40_extract_markers:
    input:
        database_pkl = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v40"], "{index}.pkl"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v40"])
    output:
        clade_marker = os.path.join(
            config["output"]["profiling"],
            "databases/strainphlan40/clade_markers/{clade}/{clade}.fna")
    log:
        os.path.join(
            config["output"]["profiling"],
            "logs/strainphlan40_extract_markers/{clade}.strainphlan40_extract_markers.log")
    benchmark:
        os.path.join(
            config["output"]["profiling"],
            "benchmark/strainphlan40_extract_markers/{clade}.strainphlan40_extract_markers.benchmark.txt")
    conda:
        config["envs"]["biobakery40"]
    params:
        clade = "{clade}",
        outdir = os.path.join(config["output"]["profiling"], "databases/strainphlan40/clade_markers/{clade}")
    priority:
        20
    threads:
        1
    shell:
        '''
        mkdir -p {params.outdir}

        extract_markers.py \
        --database {input.database_pkl} \
        --clade {params.clade} \
        --output_dir {params.outdir}/ \
        > {log} 2>&1
        '''


rule profiling_strainphlan40:
    input:
        database_pkl = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v40"], "{index}.pkl"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v40"]),
        clade_marker = os.path.join(
            config["output"]["profiling"],
            "databases/strainphlan40/clade_markers/{clade}/{clade}.fna"),
        consensus_markers = expand(os.path.join(
            config["output"]["profiling"],
            "profile/strainphlan40/consensus_markers/{sample}.pkl"),
            sample=SAMPLES_ID_LIST)
    output:
        done = os.path.join(config["output"]["profiling"], "profile/strainphlan40/clade_markers/{clade}/done")
    log:
        os.path.join(
            config["output"]["profiling"],
            "logs/strainphlan40/{clade}.strainphlan40.log")
    benchmark:
        os.path.join(
            config["output"]["profiling"],
            "benchmark/strainphlan40/{clade}.strainphlan40.benchmark.txt")
    conda:
        config["envs"]["biobakery40"]
    params:
        clade = "{clade}",
        outdir = os.path.join(config["output"]["profiling"], "profile/strainphlan40/clade_markers/{clade}"),
        marker_in_n_samples = config["params"]["profiling"]["strainphlan"]["marker_in_n_samples"],
        sample_with_n_markers = config["params"]["profiling"]["strainphlan"]["sample_with_n_markers"],
        breadth_thres = config["params"]["profiling"]["strainphlan"]["breadth_thres"],
        phylophlan_mode = config["params"]["profiling"]["strainphlan"]["phylophlan_mode"],
        reference_opts = "--references %s" % \
        os.path.join(config["output"]["profiling"], "databases/strainphlan40/reference_genomes/{clade}.fna") \
        if config["params"]["profiling"]["strainphlan"]["reference_genome"]["use"] \
        else "",
        opts = config["params"]["profiling"]["strainphlan"]["external_opts_v40"]
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
        --marker_in_n_samples {params.marker_in_n_samples} \
        --sample_with_n_markers {params.sample_with_n_markers} \
        --breadth_thres {params.breadth_thres} \
        --phylophlan_mode {params.phylophlan_mode} \
        --mutation_rates \
        {params.opts} \
        >{log} 2>&1

        touch {output.done}
        '''


if config["params"]["profiling"]["strainphlan"]["do_v40"]:
    if config["params"]["profiling"]["strainphlan"]["reference_genome"]["use"]:
        rule profiling_strainphlan40_all:
            input:
                expand(os.path.join(
                    config["output"]["profiling"],
                    "profile/strainphlan40/clade_markers/{clade}/done"),
                    clade=STRAINPHLAN_CLADES_LIST_V40)
    else:
        rule profiling_strainphlan40_all:
            input:
                aggregate_profiling_strainphlan40_extract_markers
else:
    rule profiling_strainphlan40_all:
        input:



rule profiling_strainphlan41_sample2markers:
    input:
        database_pkl = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v41"], "{index}.pkl"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v41"]),
        sam = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan41/{sample}/{sample}.sam.bz2"),
        aln = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan41/{sample}/{sample}.bowtie2.bz2")
    output:
        os.path.join(
            config["output"]["profiling"],
            "profile/strainphlan41/consensus_markers/{sample}.pkl")
    log:
        os.path.join(
            config["output"]["profiling"],
            "logs/strainphlan41_sample2markers/{sample}.strainphlan41_sample2markers.log")
    benchmark:
        os.path.join(
            config["output"]["profiling"],
            "benchmark/strainphlan41_sample2markers/{sample}.strainphlan41_sample2markers.benchmark.txt")
    params:
        outdir =  os.path.join(config["output"]["profiling"], "profile/strainphlan41/consensus_markers")
    conda:
        config["envs"]["biobakery41"]
    priority:
        20
    threads:
        config["params"]["profiling"]["threads"]
    shell:
        '''
        sample2markers.py \
        --database {input.database_pkl} \
        --input {input.sam} \
        --output_dir {params.outdir}/ \
        --nprocs {threads} \
        > {log} 2>&1
        '''


if config["params"]["profiling"]["strainphlan"]["reference_genome"]["use"]:

    STRAINPHLAN_CLADES_V41 = \
        pd.read_csv(config["params"]["profiling"]["strainphlan"]["clades_tsv_v41"], sep="\t")\
        .set_index("clade")
    STRAINPHLAN_CLADES_LIST_V41 = STRAINPHLAN_CLADES_V41.index.unique()


    rule profiling_strainphlan41_prepare_reference_genome:
        input:
            reference_genome = lambda wildcards: STRAINPHLAN_CLADES_V41.loc[wildcards.clade, "fna_path"]
        output:
            reference_genome = os.path.join(
                config["output"]["profiling"],
                "databases/strainphlan41/reference_genomes/{clade}.fna")
        threads:
            1
        run:
            if input.reference_genome.endswith(".gz"):
                shell(f'''pigz -fkdc {input.reference_genome} > {output.reference_genome}''')
            else:
                fna_path = os.path.realpath(input.reference_genome)
                fna_path_ = os.path.realpath(output.reference_genome)
                shellf('''ln -s {fna_path} {fna_path_}''')

else:
    checkpoint profiling_strainphlan41_print_clades:
        input:
            database_pkl = expand(os.path.join(
                config["params"]["profiling"]["metaphlan"]["bowtie2db_v41"], "{index}.pkl"),
                index = config["params"]["profiling"]["metaphlan"]["index_prefix_v41"]),
            consensus_markers = expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan41/consensus_markers/{sample}.pkl"),
                sample=SAMPLES_ID_LIST)
        output:
            markers_txt = os.path.join(config["output"]["profiling"], "databases/strainphlan41/clade_markers.txt"),
            markers_txt_gt = expand(
                os.path.join(config["output"]["profiling"], "databases/strainphlan41/clade_markers.txt.gt{num}"),
                num=config["params"]["profiling"]["strainphlan"]["clade_detected_in_min_samples_num"]),
            markers_dir = directory(os.path.join(config["output"]["profiling"], "databases/strainphlan41/clade_markers"))
        log:
            os.path.join(config["output"]["profiling"], "logs/strainphlan41_print_clades/strainphlan41_print_clades.log")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/strainphlan41_print_clades.txt")
        params:
            marker_in_n_samples = config["params"]["profiling"]["strainphlan"]["marker_in_n_samples"],
            clade_detected_in_min_samples_num = config["params"]["profiling"]["strainphlan"]["clade_detected_in_min_samples_num"]
        conda:
            config["envs"]["biobakery41"]
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            OUTDIR=$(dirname {output.markers_txt})
            rm -rf $OUTDIR
            mkdir -p $OUTDIR

            strainphlan \
            --nprocs {threads} \
            --database {input.database_pkl} \
            --samples {input.consensus_markers} \
            --marker_in_n_samples {params.marker_in_n_samples} \
            --output_dir $OUTDIR \
            --print_clades_only \
            >{output.markers_txt} 2>{log}

            cat {output.markers_txt} | \
            grep "t__" | \
            awk -F'[\t]' '{{print $2}}' | \
            awk -F'[: ]' '{{print $1  "\t" $4}}' \
            awk -v samnum={params.clade_detected_in_min_samples_num} '$2>=samnum{{print $1}}' \
            > {output.markers_txt_gt}

            mkdir -p {output.markers_dir}

            cat {output.markers_txt_gt} | \
            xargs -I XXX mkdir -p {output.markers_dir}/XXX
            '''


    def aggregate_profiling_strainphlan41_extract_markers(wildcards):
        checkpoint_output = checkpoints.profiling_strainphlan41_print_clades.get(**wildcards).output.markers_dir

        clades = []
        for i in glob_wildcards(os.path.join(checkpoint_output, "{clade}")).clade:
            if "t__SGB" in i:
                clades.append(i.split("/")[0])
        clades = list(set(clades))

        from pprint import pprint
        pprint(clades)

        clades_done = expand(os.path.join(
            config["output"]["profiling"],
            "profile/strainphlan41/clade_markers/{clade}/done"),
            clade = clades)

        return clades_done


rule profiling_strainphlan41_extract_markers:
    input:
        database_pkl = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v41"], "{index}.pkl"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v41"])
    output:
        clade_marker = os.path.join(
            config["output"]["profiling"],
            "databases/strainphlan41/clade_markers/{clade}/{clade}.fna")
    log:
        os.path.join(
            config["output"]["profiling"],
            "logs/strainphlan41_extract_markers/{clade}.strainphlan41_extract_markers.log")
    benchmark:
        os.path.join(
            config["output"]["profiling"],
            "benchmark/strainphlan41_extract_markers/{clade}.strainphlan41_extract_markers.benchmark.txt")
    conda:
        config["envs"]["biobakery41"]
    params:
        clade = "{clade}",
        outdir = os.path.join(config["output"]["profiling"], "databases/strainphlan41/clade_markers/{clade}")
    priority:
        20
    threads:
        1
    shell:
        '''
        mkdir -p {params.outdir}

        extract_markers.py \
        --database {input.database_pkl} \
        --clade {params.clade} \
        --output_dir {params.outdir}/ \
        > {log} 2>&1
        '''


rule profiling_strainphlan41:
    input:
        database_pkl = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v41"], "{index}.pkl"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v41"]),
        clade_marker = os.path.join(
            config["output"]["profiling"],
            "databases/strainphlan41/clade_markers/{clade}/{clade}.fna"),
        consensus_markers = expand(os.path.join(
            config["output"]["profiling"],
            "profile/strainphlan41/consensus_markers/{sample}.pkl"),
            sample=SAMPLES_ID_LIST)
    output:
        done = os.path.join(config["output"]["profiling"], "profile/strainphlan41/clade_markers/{clade}/done")
    log:
        os.path.join(
            config["output"]["profiling"],
            "logs/strainphlan41/{clade}.strainphlan41.log")
    benchmark:
        os.path.join(
            config["output"]["profiling"],
            "benchmark/strainphlan41/{clade}.strainphlan41.benchmark.txt")
    conda:
        config["envs"]["biobakery41"]
    params:
        clade = "{clade}",
        outdir = os.path.join(config["output"]["profiling"], "profile/strainphlan41/clade_markers/{clade}"),
        marker_in_n_samples = config["params"]["profiling"]["strainphlan"]["marker_in_n_samples"],
        sample_with_n_markers = config["params"]["profiling"]["strainphlan"]["sample_with_n_markers"],
        breadth_thres = config["params"]["profiling"]["strainphlan"]["breadth_thres"],
        phylophlan_mode = config["params"]["profiling"]["strainphlan"]["phylophlan_mode"],
        reference_opts = "--references %s" % \
        os.path.join(config["output"]["profiling"], "databases/strainphlan41/reference_genomes/{clade}.fna") \
        if config["params"]["profiling"]["strainphlan"]["reference_genome"]["use"] \
        else "",
        opts = config["params"]["profiling"]["strainphlan"]["external_opts_v41"]
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
        --marker_in_n_samples {params.marker_in_n_samples} \
        --sample_with_n_markers {params.sample_with_n_markers} \
        --breadth_thres {params.breadth_thres} \
        --phylophlan_mode {params.phylophlan_mode} \
        --mutation_rates \
        {params.opts} \
        >{log} 2>&1

        touch {output.done}
        '''


if config["params"]["profiling"]["strainphlan"]["do_v41"]:
    if config["params"]["profiling"]["strainphlan"]["reference_genome"]["use"]:
        rule profiling_strainphlan41_all:
            input:
                expand(os.path.join(
                    config["output"]["profiling"],
                    "profile/strainphlan41/clade_markers/{clade}/done"),
                    clade=STRAINPHLAN_CLADES_LIST_V41)
    else:
        rule profiling_strainphlan41_all:
            input:
                aggregate_profiling_strainphlan41_extract_markers
else:
    rule profiling_strainphlan41_all:
        input:
