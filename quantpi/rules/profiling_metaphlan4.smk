if config["params"]["profiling"]["metaphlan"]["do_v4"]:
    rule profiling_metaphlan4:
        input:
            reads = profiling_input_with_short_reads,
            index = expand(
                os.path.join(
                    config["params"]["profiling"]["metaphlan"]["bowtie2db"],
                    "{index}.{suffix}"),
                index = config["params"]["profiling"]["metaphlan"]["index_v4"],
                suffix = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "pkl"])
        output:
            profile = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan4/{sample}/{sample}.metaphlan4.abundance.profile.tsv"),
            samout = temp(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan4/{sample}/{sample}.sam.bz2")) \
                if config["params"]["profiling"]["metaphlan"]["no_sam"] \
                else os.path.join(config["output"]["profiling"],
                "profile/metaphlan4/{sample}/{sample}.sam.bz2"),
            mapout = temp(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan4/{sample}/{sample}.bowtie2.bz2")) \
                if config["params"]["profiling"]["metaphlan"]["no_map"] \
                else os.path.join(config["output"]["profiling"],
                "profile/metaphlan4/{sample}/{sample}.bowtie2.bz2")
        log:
            os.path.join(config["output"]["profiling"], "logs/metaphlan4/{sample}.metaphlan4.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/metaphlan4/{sample}.metaphlan4.benchmark.txt")
        conda:
            config["envs"]["metaphlan4"]
        params:
            index = config["params"]["profiling"]["metaphlan"]["index_v4"],
            bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
            sample_id = "{sample}",
            outdir = os.path.join(config["output"]["profiling"], "profile/metaphlan4/{sample}"),
            read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
            min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
            stat_q = config["params"]["profiling"]["metaphlan"]["stat_q_v4"],
            stat = config["params"]["profiling"]["metaphlan"]["stat"],
            taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
            analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
            avoid_disqm = "--avoid_disqm" if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] else "",
            ignore_markers = "--ignore_markers %s" % config["params"]["profiling"]["metaphlan"]["ignore_markers"] \
                if os.path.exists(config["params"]["profiling"]["metaphlan"]["ignore_markers"]) \
                else "",
            unknown_estimation = "--unclassified_estimation" \
                if config["params"]["profiling"]["metaphlan"]["unknown_estimation"] \
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
            ignore_ksgbs = "--ignore_ksgbs" \
                if config["params"]["profiling"]["metaphlan"]["ignore_ksgbs"] \
                else "",
            ignore_usgbs = "--ignore_usgbs" \
                if config["params"]["profiling"]["metaphlan"]["ignore_usgbs"] \
                else "",
            biom_out = "--biom %s" % os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan4/{sample}/{sample}.metaphlan4.abundance.profile.biom") \
                if config["params"]["profiling"]["metaphlan"]["biom"] \
                else "",
            legacy_output = "--legacy-output" \
                if config["params"]["profiling"]["metaphlan"]["legacy_output"] \
                else "",
            cami_format_output = "--CAMI_format_output" \
                if config["params"]["profiling"]["metaphlan"]["cami_format_output"] \
                else "",
            external_opts = config["params"]["profiling"]["metaphlan"]["external_opts_v4"]
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            mkdir -p {params.outdir}

            reads=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.reads})

            metaphlan \
            $reads \
            --force \
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
            {params.unknown_estimation} \
            {params.ignore_eukaryotes} \
            {params.ignore_bacteria} \
            {params.ignore_archaea} \
            {params.ignore_ksgbs} \
            {params.ignore_usgbs} \
            {params.biom_out} \
            {params.legacy_output} \
            {params.cami_format_output} \
            {params.external_opts} \
            --sample_id_key {params.sample_id} \
            --sample_id {params.sample_id} \
            --bowtie2out {output.mapout} \
            --samout {output.samout} \
            --output_file {output} \
            2> {log}
            '''


    rule profiling_metaphlan4_merge:
        input:
           abuns = expand(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan4/{sample}/{sample}.metaphlan4.abundance.profile.tsv"),
                          sample=SAMPLES_ID_LIST)
        output:
            profiles = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/metaphlan4/metaphlan4.merged.abundance.profile.{level}.tsv"),
                level=["all", "strain", "species", "genus", "family",
                       "order", "class", "phylum", "superkingdom"])
        threads:
            config["params"]["profiling"]["threads"]
        priority:
            20
        run:
            quantpi.metaphlan_init(4)
            profile_list = quantpi.merge_metaphlan_tables(input.abuns, threads)
            for i in range(0, len(profile_list)):
                profile_list[i].to_csv(output.profiles[i], sep='\t', index=False)


    rule profiling_metaphlan4_all:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/metaphlan4/metaphlan4.merged.abundance.profile.{level}.tsv"),
                level=["all", "superkingdom", "phylum", "class",
                       "order", "family", "genus", "species", "strain"]),
                
            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_metaphlan4_all:
        input:
