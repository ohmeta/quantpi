if config["params"]["profiling"]["metaphlan"]["do_v2"]:
    rule profiling_metaphlan2:
        input:
            reads = profiling_input_with_short_reads,
            index = expand(
                os.path.join(
                    config["params"]["profiling"]["metaphlan"]["bowtie2db"],
                    "{index}.{suffix}"),
                index = config["params"]["profiling"]["metaphlan"]["index_v2"],
                suffix = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "pkl"])
        output:
            profile = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"),
            samout = temp(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.sam.bz2")) \
                if config["params"]["profiling"]["metaphlan"]["no_sam"] \
                else os.path.join(config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.sam.bz2"),
            mapout = temp(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.bowtie2.bz2")) \
                if config["params"]["profiling"]["metaphlan"]["no_map"] \
                else os.path.join(config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.bowtie2.bz2")
        conda:
            config["envs"]["biobakery2"]
        log:
            os.path.join(config["output"]["profiling"], "logs/metaphlan2/{sample}.metaphlan2.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/metaphlan2/{sample}.metaphlan2.benchmark.txt")
        params:
            bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            index = config["params"]["profiling"]["metaphlan"]["index_v2"],
            bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
            sample_id = "{sample}",
            outdir = os.path.join(config["output"]["profiling"], "profile/metaphlan2/{sample}"),
            read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
            min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
            stat_q = config["params"]["profiling"]["metaphlan"]["stat_q_v2"],
            stat = config["params"]["profiling"]["metaphlan"]["stat"],
            taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
            analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
            avoid_disqm = "--avoid_disqm" \
                if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] \
                else "",
            ignore_markers = "--ignore_markers %s" % config["params"]["profiling"]["metaphlan"]["ignore_markers"] \
                if os.path.exists(config["params"]["profiling"]["metaphlan"]["ignore_markers"]) \
                else "",
            ignore_viruses = "--ignore_viruses" \
                if not config["params"]["profiling"]["metaphlan"]["add_viruses"] \
                else "",
            ignore_eukaryotes = "--ignore_eukaryotes" \
                if config["params"]["profiling"]["metaphlan"]["ignore_eukaryotes"] \
                else "",
            ignore_bacteria = "--ignore_bacteria" \
                if config["params"]["profiling"]["metaphlan"]["ignore_bacteria"] \
                else "",
            ignore_archaea = "--ignore_archaea" \
                if config["params"]["profiling"]["metaphlan"]["ignore_archaea"] \
                else "",
            biom_out = "--biom %s" % os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.biom") \
                if config["params"]["profiling"]["metaphlan"]["biom"] \
                else "",
            external_opts = config["params"]["profiling"]["metaphlan"]["external_opts_v2"]
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
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
            --stat_q {params.stat_q} \
            --stat {params.stat} \
            --tax_lev {params.taxonomic_level} \
            -t {params.analysis_type} \
            {params.avoid_disqm} \
            {params.ignore_markers} \
            {params.ignore_viruses} \
            {params.ignore_eukaryotes} \
            {params.ignore_bacteria} \
            {params.ignore_archaea} \
            {params.biom_out} \
            {params.external_opts} \
            --sample_id {params.sample_id} \
            --sample_id_key {params.sample_id} \
            --bowtie2out {output.mapout} \
            --samout {output.samout} \
            --output_file {output} \
            2> {log}
            '''


    rule profiling_metaphlan2_merge:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"),
                   sample=SAMPLES_ID_LIST)
        output:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/metaphlan2/metaphlan2.merged.abundance.profile.{level}.tsv"),
                level=["all", "strain", "species", "genus", "family",
                       "order", "class", "phylum", "superkingdom"])
        threads:
            config["params"]["profiling"]["threads"]
        priority:
            20
        run:
           quantpi.metaphlan_init(2)
           df_list = quantpi.merge_metaphlan_tables(input, threads)
           for i in range(0, len(df_list)):
               df_list[i].to_csv(output[i], sep='\t', index=False)


    rule profiling_metaphlan2_all:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/metaphlan2/metaphlan2.merged.abundance.profile.{level}.tsv"),
                level=["all", "superkingdom", "phylum", "class",
                       "order", "family", "genus", "species", "strain"])

else:
    rule profiling_metaphlan2_all:
        input:
