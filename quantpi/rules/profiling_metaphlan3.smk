
if config["params"]["profiling"]["metaphlan"]["do_v3"]:
    rule profiling_metaphlan3:
        input:
            reads = profiling_input_with_short_reads,
            index = expand(
                os.path.join(
                    config["params"]["profiling"]["metaphlan"]["bowtie2db"],
                    "{index}.{suffix}"),
                index = config["params"]["profiling"]["metaphlan"]["index_v3"],
                suffix = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "pkl"])
        output:
            profile = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv"),
            samout = temp(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.sam")) \
                if config["params"]["profiling"]["metaphlan"]["no_sam"] \
                else os.path.join(config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.sam"),
            mapout = temp(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.bowtie2.bz2")) \
                if config["params"]["profiling"]["metaphlan"]["no_map"] \
                else os.path.join(config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.bowtie2.bz2")
        log:
            os.path.join(config["output"]["profiling"], "logs/metaphlan3/{sample}.metaphlan3.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/metaphlan3/{sample}.metaphlan3.benchmark.txt")
        conda:
            config["envs"]["metaphlan3"]
        params:
            bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            index = config["params"]["profiling"]["metaphlan"]["index_v3"],
            bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
            sample_id = "{sample}",
            outdir = os.path.join(config["output"]["profiling"], "profile/metaphlan3/{sample}"),
            read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
            min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
            stat_q = config["params"]["profiling"]["metaphlan"]["stat_q_v3"],
            stat = config["params"]["profiling"]["metaphlan"]["stat"],
            taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
            analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
            avoid_disqm = "--avoid_disqm" \
                if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] \
                else "",
            ignore_markers = "--ignore_markers %s" % config["params"]["profiling"]["metaphlan"]["ignore_markers"] \
                if os.path.exists(config["params"]["profiling"]["metaphlan"]["ignore_markers"]) \
                else "",
            unknown_estimation = "--unknown_estimation" \
                if config["params"]["profiling"]["metaphlan"]["unknown_estimation"] \
                else "",
            add_viruses = "--add_viruses" \
                if config["params"]["profiling"]["metaphlan"]["add_viruses"] \
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
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.biom") \
                if config["params"]["profiling"]["metaphlan"]["biom"] \
                else "",
            legacy_output = "--legacy-output" \
                if config["params"]["profiling"]["metaphlan"]["legacy_output"] \
                else "",
            cami_format_output = "--CAMI_format_output" \
                if config["params"]["profiling"]["metaphlan"]["cami_format_output"] \
                else "",
            external_opts = config["params"]["profiling"]["metaphlan"]["external_opts_v3"]
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
            {params.add_viruses} \
            {params.ignore_eukaryotes} \
            {params.ignore_bacteria} \
            {params.ignore_archaea} \
            {params.biom_out} \
            {params.legacy_output} \
            {params.cami_format_output} \
            {params.external_opts} \
            --sample_id {params.sample_id} \
            --sample_id_key {params.sample_id} \
            --bowtie2out {output.mapout} \
            --samout {output.samout} \
            --output_file {output.profile} \
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
