rule profiling_metaphlan2:
    input:
        reads = profiling_input_with_short_reads,
        index = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v2"], "{index}.{suffix}"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v2"],
            suffix = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "pkl"])
    output:
        profile = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"),
        samout = os.path.join(config["output"]["profiling"], "profile/metaphlan2/{sample}/{sample}.sam.bz2"),
        mapout = os.path.join(config["output"]["profiling"], "profile/metaphlan2/{sample}/{sample}.bowtie2.bz2")
    log:
        os.path.join(config["output"]["profiling"], "logs/metaphlan2/{sample}.metaphlan2.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/metaphlan2/{sample}.metaphlan2.benchmark.txt")
    params:
        sample_id = "{sample}",
        bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db_v2"],
        index = config["params"]["profiling"]["metaphlan"]["index_prefix_v2"],
        bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
        outdir = os.path.join(config["output"]["profiling"], "profile/metaphlan2/{sample}"),
        read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
        min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
        stat = config["params"]["profiling"]["metaphlan"]["stat"],
        taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
        analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
        external_opts = config["params"]["profiling"]["metaphlan"]["external_opts_v2"]
    priority:
        20
    threads:
        config["params"]["profiling"]["threads"]
    resources:
        mem_mb=config["params"]["profiling"]["metaphlan"]["mem_mb"]
    conda:
        config["envs"]["biobakery2"]
    shell:
        '''
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        reads=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.reads})

        metaphlan2.py \
        $reads \
        --nproc {threads} \
        --input_type fastq \
        --bowtie2db {params.bowtie2db} \
        --index {params.index} \
        --bt2_ps {params.bowtie2_presets} \
        --read_min_len {params.read_min_len} \
        --min_cu_len {params.min_cu_len} \
        --stat {params.stat} \
        --tax_lev {params.taxonomic_level} \
        -t {params.analysis_type} \
        --sample_id {params.sample_id} \
        --sample_id_key {params.sample_id} \
        {params.external_opts} \
        --bowtie2out {output.mapout} \
        --samout {output.samout} \
        --output_file {output.profile} \
        >{log} 2>&1
        '''


rule profiling_metaphlan2_merge:
    input:
        expand(os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"),
            sample=SAMPLES_ID_LIST)
    output:
        expand(os.path.join(
            config["output"]["profiling"],
            "report/metaphlan2/metaphlan2.merged.abundance.profile.{level}.tsv"),
            level=["all", "strain", "species", "genus", "family",
                    "order", "class", "phylum", "superkingdom"])
    threads:
        config["params"]["profiling"]["threads"]
    resources:
        mem_mb=config["params"]["profiling"]["metaphlan"]["mem_mb"]
    priority:
        20
    run:
        quantpi.metaphlan_init(2)
        df_list = quantpi.merge_metaphlan_tables(input, threads)
        for i in range(0, len(df_list)):
            df_list[i].to_csv(output[i], sep='\t', index=False)


if config["params"]["profiling"]["metaphlan"]["do_v2"]:
    rule profiling_metaphlan2_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "report/metaphlan2/metaphlan2.merged.abundance.profile.{level}.tsv"),
                level=["all", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"])
else:
    rule profiling_metaphlan2_all:
        input:



rule profiling_metaphlan3:
    input:
        reads = profiling_input_with_short_reads,
        index = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v3"],
            "{index}.{suffix}"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v3"],
            suffix = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "pkl"])
    output:
        profile = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv"),
        samout = os.path.join(config["output"]["profiling"], "profile/metaphlan3/{sample}/{sample}.sam.bz2"),
        mapout = os.path.join(config["output"]["profiling"], "profile/metaphlan3/{sample}/{sample}.bowtie2.bz2")
    log:
        os.path.join(config["output"]["profiling"], "logs/metaphlan3/{sample}.metaphlan3.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/metaphlan3/{sample}.metaphlan3.benchmark.txt")
    params:
        sample_id = "{sample}",
        bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db_v3"],
        index = config["params"]["profiling"]["metaphlan"]["index_prefix_v3"],
        bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
        outdir = os.path.join(config["output"]["profiling"], "profile/metaphlan3/{sample}"),
        read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
        min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
        stat = config["params"]["profiling"]["metaphlan"]["stat"],
        taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
        analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
        external_opts = config["params"]["profiling"]["metaphlan"]["external_opts_v3"]
    priority:
        20
    threads:
        config["params"]["profiling"]["threads"]
    resources:
        mem_mb=config["params"]["profiling"]["metaphlan"]["mem_mb"]
    conda:
        config["envs"]["biobakery3"]
    shell:
        '''
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        reads=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.reads})

        metaphlan \
        $reads \
        --nproc {threads} \
        --input_type fastq \
        --bowtie2db {params.bowtie2db} \
        --index {params.index} \
        --bt2_ps {params.bowtie2_presets} \
        --read_min_len {params.read_min_len} \
        --min_cu_len {params.min_cu_len} \
        --stat {params.stat} \
        --tax_lev {params.taxonomic_level} \
        -t {params.analysis_type} \
        --sample_id {params.sample_id} \
        --sample_id_key {params.sample_id} \
        {params.external_opts} \
        --bowtie2out {output.mapout} \
        --samout {output.samout} \
        --output_file {output.profile} \
        >{log} 2>&1
        '''


rule profiling_metaphlan3_merge:
    input:
        abuns = expand(os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv"),
            sample=SAMPLES_ID_LIST)
    output:
        profiles = expand(os.path.join(
            config["output"]["profiling"],
            "report/metaphlan3/metaphlan3.merged.abundance.profile.{level}.tsv"),
            level=["all", "strain", "species", "genus", "family", "order", "class", "phylum", "superkingdom"])
    threads:
        config["params"]["profiling"]["threads"]
    resources:
        mem_mb=config["params"]["profiling"]["metaphlan"]["mem_mb"]
    priority:
        20
    run:
        quantpi.metaphlan_init(3)
        profile_list = quantpi.merge_metaphlan_tables(input.abuns, threads)
        for i in range(0, len(profile_list)):
            profile_list[i].to_csv(output.profiles[i], sep='\t', index=False)


if config["params"]["profiling"]["metaphlan"]["do_v3"]:
    rule profiling_metaphlan3_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "report/metaphlan3/metaphlan3.merged.abundance.profile.{level}.tsv"),
                level=["all", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"])
else:
    rule profiling_metaphlan3_all:
        input:



rule profiling_metaphlan40:
    input:
        reads = profiling_input_with_short_reads,
        index = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v40"],
            "{index}.{suffix}"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v40"],
            suffix = ["1.bt2l", "2.bt2l", "3.bt2l", "4.bt2l", "rev.1.bt2l", "rev.2.bt2l", "pkl"])
    output:
        profile = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan40/{sample}/{sample}.metaphlan40.abundance.profile.tsv"),
        samout = os.path.join(config["output"]["profiling"], "profile/metaphlan40/{sample}/{sample}.sam.bz2"),
        mapout = os.path.join(config["output"]["profiling"], "profile/metaphlan40/{sample}/{sample}.bowtie2.bz2")
    log:
        os.path.join(config["output"]["profiling"], "logs/metaphlan40/{sample}.metaphlan40.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/metaphlan40/{sample}.metaphlan40.benchmark.txt")
    params:
        sample_id = "{sample}",
        bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db_v40"],
        index = config["params"]["profiling"]["metaphlan"]["index_prefix_v40"],
        bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
        outdir = os.path.join(config["output"]["profiling"], "profile/metaphlan40/{sample}"),
        read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
        min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
        stat = config["params"]["profiling"]["metaphlan"]["stat"],
        taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
        analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
        external_opts = config["params"]["profiling"]["metaphlan"]["external_opts_v40"]
    priority:
        20
    threads:
        config["params"]["profiling"]["threads"]
    resources:
        mem_mb=config["params"]["profiling"]["metaphlan"]["mem_mb"]
    conda:
        config["envs"]["biobakery40"]
    shell:
        '''
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        reads=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.reads})

        metaphlan \
        $reads \
        --nproc {threads} \
        --input_type fastq \
        --bowtie2db {params.bowtie2db} \
        --index {params.index} \
        --bt2_ps {params.bowtie2_presets} \
        --read_min_len {params.read_min_len} \
        --min_cu_len {params.min_cu_len} \
        --stat {params.stat} \
        --tax_lev {params.taxonomic_level} \
        -t {params.analysis_type} \
        --sample_id {params.sample_id} \
        --sample_id_key {params.sample_id} \
        {params.external_opts} \
        --bowtie2out {output.mapout} \
        --samout {output.samout} \
        --output_file {output.profile} \
        --offline \
        >{log} 2>&1
        '''


rule profiling_metaphlan40_merge:
    input:
        abuns = expand(os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan40/{sample}/{sample}.metaphlan40.abundance.profile.tsv"),
            sample=SAMPLES_ID_LIST)
    output:
        profiles = expand(os.path.join(
            config["output"]["profiling"],
            "report/metaphlan40/metaphlan40.merged.abundance.profile.{level}.tsv"),
            level=["all", "strain", "species", "genus", "family", "order", "class", "phylum", "superkingdom"])
    threads:
        config["params"]["profiling"]["threads"]
    resources:
        mem_mb=config["params"]["profiling"]["metaphlan"]["mem_mb"]
    priority:
        20
    run:
        quantpi.metaphlan_init(4)
        profile_list = quantpi.merge_metaphlan_tables(input.abuns, threads)
        for i in range(0, len(profile_list)):
            profile_list[i].to_csv(output.profiles[i], sep='\t', index=False)


if config["params"]["profiling"]["metaphlan"]["do_v40"]:
    rule profiling_metaphlan40_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "report/metaphlan40/metaphlan40.merged.abundance.profile.{level}.tsv"),
                level=["all", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"])
else:
    rule profiling_metaphlan40_all:
        input:



rule profiling_metaphlan41:
    input:
        reads = profiling_input_with_short_reads,
        index = expand(os.path.join(
            config["params"]["profiling"]["metaphlan"]["bowtie2db_v41"], "{index}.{suffix}"),
            index = config["params"]["profiling"]["metaphlan"]["index_prefix_v41"],
            suffix = ["1.bt2l", "2.bt2l", "3.bt2l", "4.bt2l", "rev.1.bt2l", "rev.2.bt2l", "pkl"])
    output:
        profile = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan41/{sample}/{sample}.metaphlan41.abundance.profile.tsv"),
        vsc = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan41/{sample}/{sample}.metaphlan41.vsc.tsv"),
        samout = os.path.join(config["output"]["profiling"], "profile/metaphlan41/{sample}/{sample}.sam.bz2"),
        mapout = os.path.join(config["output"]["profiling"], "profile/metaphlan41/{sample}/{sample}.bowtie2.bz2")
    log:
        os.path.join(config["output"]["profiling"], "logs/metaphlan41/{sample}.metaphlan41.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/metaphlan41/{sample}.metaphlan41.benchmark.txt")
    params:
        sample_id = "{sample}",
        bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db_v41"],
        index = config["params"]["profiling"]["metaphlan"]["index_prefix_v41"],
        bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
        outdir = os.path.join(config["output"]["profiling"], "profile/metaphlan41/{sample}"),
        read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
        min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
        stat = config["params"]["profiling"]["metaphlan"]["stat"],
        taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
        analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
        external_opts = config["params"]["profiling"]["metaphlan"]["external_opts_v41"]
    priority:
        20
    threads:
        config["params"]["profiling"]["threads"]
    resources:
        mem_mb=config["params"]["profiling"]["metaphlan"]["mem_mb"]
    conda:
        config["envs"]["biobakery41"]
    shell:
        '''
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        reads=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.reads})

        metaphlan \
        $reads \
        --nproc {threads} \
        --input_type fastq \
        --bowtie2db {params.bowtie2db} \
        --index {params.index} \
        --bt2_ps {params.bowtie2_presets} \
        --read_min_len {params.read_min_len} \
        --min_cu_len {params.min_cu_len} \
        --stat {params.stat} \
        --tax_lev {params.taxonomic_level} \
        -t {params.analysis_type} \
        --sample_id {params.sample_id} \
        --sample_id_key {params.sample_id} \
        {params.external_opts} \
        --profile_vsc \
        --vsc_out {output.vsc} \
        --bowtie2out {output.mapout} \
        --samout {output.samout} \
        --output_file {output.profile} \
        --offline \
        >{log} 2>&1
        '''


rule profiling_metaphlan41_merge:
    input:
        abuns = expand(os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan41/{sample}/{sample}.metaphlan41.abundance.profile.tsv"),
            sample=SAMPLES_ID_LIST)
    output:
        profiles = expand(os.path.join(
            config["output"]["profiling"],
            "report/metaphlan41/metaphlan41.merged.abundance.profile.{level}.tsv"),
            level=["all", "strain", "species", "genus", "family", "order", "class", "phylum", "superkingdom"])
    threads:
        config["params"]["profiling"]["threads"]
    resources:
        mem_mb=config["params"]["profiling"]["metaphlan"]["mem_mb"]
    priority:
        20
    run:
        quantpi.metaphlan_init(4)
        profile_list = quantpi.merge_metaphlan_tables(input.abuns, threads)
        for i in range(0, len(profile_list)):
            profile_list[i].to_csv(output.profiles[i], sep='\t', index=False)


if config["params"]["profiling"]["metaphlan"]["do_v41"]:
    rule profiling_metaphlan41_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "report/metaphlan41/metaphlan41.merged.abundance.profile.{level}.tsv"),
                level=["all", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"])
else:
    rule profiling_metaphlan41_all:
        input:
