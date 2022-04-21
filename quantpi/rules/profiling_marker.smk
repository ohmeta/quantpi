if config["params"]["profiling"]["metaphlan"]["do_v2"]:
    rule profiling_metaphlan2:
        input:
            reads = profiling_input_with_short_reads
        output:
           profile = protected(os.path.join(
               config["output"]["profiling"],
               "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"))
        conda:
            config["envs"]["metaphlan2"]
        log:
            os.path.join(config["output"]["profiling"], "logs/metaphlan2/{sample}.metaphlan2.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/metaphlan2/{sample}.metaphlan2.benchmark.txt")
        params:
            sample_id = "{sample}",
            wrapper_dir = WRAPPER_DIR,
            input_str = lambda wildcards: ",".join(profiling_input_with_short_reads(wildcards)),
            read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
            bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            index = config["params"]["profiling"]["metaphlan"]["index_v2"],
            bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
            min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
            taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
            avoid_disqm = "--avoid_disqm" \
                if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] \
                   else "",
            stat_q = config["params"]["profiling"]["metaphlan"]["stat_q_v2"],
            stat = config["params"]["profiling"]["metaphlan"]["stat"],
            analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
            bowtie2out = os.path.join(
               config["output"]["profiling"],
               "profile/metaphlan2/{sample}/{sample}.metaphlan2.bowtie2.bz2")
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            python {params.wrapper_dir}/metaphlan2_wrapper.py \
            --input {input.reads} \
            --analysis_type {params.analysis_type} \
            --input_type multifastq \
            --bowtie2db {params.bowtie2db} \
            --index {params.index} \
            --bt2_ps {params.bowtie2_presets} \
            --bowtie2out {params.bowtie2out} \
            --tax_lev {params.taxonomic_level} \
            --min_cu_len {params.min_cu_len} \
            --stat_q {params.stat_q} \
            --output_file {output.profile} \
            --sample_id {params.sample_id} \
            --nproc {threads} \
            --read_min_len {params.read_min_len} \
            >{log} 2>&1
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
                       "order", "family", "genus", "species", "strain"]),

            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_metaphlan2_all:
        input:


if config["params"]["profiling"]["metaphlan"]["do_v3"]:
    rule profiling_metaphlan3:
        input:
            reads = profiling_input_with_short_reads
            #index_mp3 = expand(
            #    os.path.join(
            #        config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            #        "{mpa_name}.{suffix}"),
            #    mpa_name = config["params"]["profiling"]["metaphlan"]["index_v3"],
            #    suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "pkl"])
        output:
            protected(os.path.join(
                config["output"]["profiling"],
                    "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv"))
        log:
            os.path.join(config["output"]["profiling"], "logs/metaphlan3/{sample}.metaphlan3.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/metaphlan3/{sample}.metaphlan3.benchmark.txt")
        conda:
            config["envs"]["metaphlan3"]
        params:
            sample_id = "{sample}",
            read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
            bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            index = config["params"]["profiling"]["metaphlan"]["index_v3"],
            bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
            min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
            taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
            ignore_markers = "--ignore_markers %s" % config["params"]["profiling"]["metaphlan"]["ignore_markers"] \
                if os.path.exists(config["params"]["profiling"]["metaphlan"]["ignore_markers"]) \
                else "",
            avoid_disqm = "--avoid_disqm" \
                if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] \
                else "",
            stat_q = config["params"]["profiling"]["metaphlan"]["stat_q_v3"],
            stat = config["params"]["profiling"]["metaphlan"]["stat"],
            analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
            unknown_estimation = "--unknown_estimation" \
                if config["params"]["profiling"]["metaphlan"]["unknown_estimation"] \
                else "",
            add_viruses = "--add_viruses" \
                if config["params"]["profiling"]["metaphlan"]["add_viruses"] \
                else "",
            map_out = "--no_map" if config["params"]["profiling"]["metaphlan"]["no_map"] \
                else "--bowtie2out %s" % os.path.join(
                    config["output"]["profiling"],
                    "profile/metaphlan3/{sample}/{sample}.metaphlan3.bowtie2.bz2"),
            biom_out = "--biom %s" % os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.biom") \
                if config["params"]["profiling"]["metaphlan"]["biom"] \
                else "",
            legacy_output = "--legacy-output" \
                if config["params"]["profiling"]["metaphlan"]["legacy_output"] \
                else "",
            cami_format_output = "--CAMI_format_output" \
                if config["params"]["profiling"]["metaphlan"]["cami_format_output"] \
                else ""
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            reads=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.reads})

            metaphlan \
            $reads \
            --input_type fastq \
            --read_min_len {params.read_min_len} \
            --nproc {threads} \
            {params.map_out} \
            --bowtie2db {params.bowtie2db} \
            --index {params.index} \
            --bt2_ps {params.bowtie2_presets} \
            --min_cu_len {params.min_cu_len} \
            --tax_lev {params.taxonomic_level} \
            {params.ignore_markers} \
            {params.avoid_disqm} \
            --stat_q {params.stat_q} \
            --stat {params.stat} \
            -t {params.analysis_type} \
            {params.unknown_estimation} \
            {params.add_viruses} \
            --output_file {output} \
            --sample_id {params.sample_id} \
            {params.legacy_output} \
            {params.cami_format_output} \
            {params.biom_out} \
            2> {log}
            '''


    rule profiling_metaphlan3_merge:
        input:
           abuns = expand(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv"),
                          sample=SAMPLES_ID_LIST)
        output:
            profiles = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/metaphlan3/metaphlan3.merged.abundance.profile.{level}.tsv"),
                level=["all", "strain", "species", "genus", "family",
                       "order", "class", "phylum", "superkingdom"])
        threads:
            config["params"]["profiling"]["threads"]
        priority:
            20
        run:
            quantpi.metaphlan_init(3)
            profile_list = quantpi.merge_metaphlan_tables(input.abuns, threads)
            for i in range(0, len(profile_list)):
                profile_list[i].to_csv(output.profiles[i], sep='\t', index=False)


    rule profiling_metaphlan3_all:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/metaphlan3/metaphlan3.merged.abundance.profile.{level}.tsv"),
                level=["all", "superkingdom", "phylum", "class",
                       "order", "family", "genus", "species", "strain"]),
                
            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_metaphlan3_all:
        input: