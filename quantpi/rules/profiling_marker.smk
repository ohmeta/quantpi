if config["params"]["profiling"]["metaphlan"]["do_v2"]:
    rule profiling_metaphlan2:
        input:
            reads = profiling_input_with_short_reads
        output:
           profile = os.path.join(
               config["output"]["profiling"],
               "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv")
        conda:
            config["envs"]["metaphlan2"]
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
            map_out = "--no_map" if config["params"]["profiling"]["metaphlan"]["no_map"] \
                else "--bowtie2out %s" % os.path.join(
                    config["output"]["profiling"],
                    "profile/metaphlan2/{sample}/{sample}.metaphlan2.bowtie2.bz2"),
            biom_out = "--biom %s" % os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.biom") \
                if config["params"]["profiling"]["metaphlan"]["biom"] \
                else ""
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
            {params.map_out} \
            {params.biom_out} \
            --sample_id {params.sample_id} \
            --sample_id_key {params.sample_id} \
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
            os.path.join(
                config["output"]["profiling"],
                    "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv")
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
            {params.map_out} \
            {params.biom_out} \
            {params.legacy_output} \
            {params.cami_format_output} \
            --sample_id {params.sample_id} \
            --sample_id_key {params.sample_id} \
            --output_file {output} \
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


if config["params"]["profiling"]["metaphlan"]["do_v4"]:
    rule profiling_metaphlan4:
        input:
            reads = profiling_input_with_short_reads
            #index_mp3 = expand(
            #    os.path.join(
            #        config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            #        "{mpa_name}.{suffix}"),
            #    mpa_name = config["params"]["profiling"]["metaphlan"]["index_v4"],
            #    suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "pkl"])
        output:
            os.path.join(
                config["output"]["profiling"],
                    "profile/metaphlan4/{sample}/{sample}.metaphlan4.abundance.profile.tsv")
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
            map_out = "--no_map" if config["params"]["profiling"]["metaphlan"]["no_map"] \
                else "--bowtie2out %s" % os.path.join(
                    config["output"]["profiling"],
                    "profile/metaphlan4/{sample}/{sample}.metaphlan4.bowtie2.bz2"),
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
                else ""
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
            {params.map_out} \
            {params.biom_out} \
            {params.legacy_output} \
            {params.cami_format_output} \
            --sample_id_key {params.sample_id} \
            --sample_id {params.sample_id} \
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