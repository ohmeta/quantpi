if config["params"]["profiling"]["bgi_soap"]["do"]:
    rule profiling_custom_bgi_soap:
        input:
            reads = profiling_input_with_short_reads,
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["profiling"]["bgi_soap"]["index_prefix"],
                           suffix=["amb", "ann", "bwt", "fmv", "hot", "lkt", "pac",
                                   "rev.bwt", "rev.fmv", "rev.lkt", "rev.pac", "sa", ".sai"])
        output:
            soap = os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap/{sample}/{sample}.bgi_soap.soap.gz")
        log:
            os.path.join(config["output"]["profiling"], "logs/bgi_soap/{sample}.bgi_soap.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/bgi_soap/{sample}.profiling.bgi_soap.benchmark.txt")
        params:
            index = config["params"]["profiling"]["bgi_soap"]["index_prefix"],
            minimal_insert_size = config["params"]["profiling"]["bgi_soap"]["minimal_insert_size"],
            maximal_insert_size = config["params"]["profiling"]["bgi_soap"]["maximal_insert_size"],
            report_repeat_hits = config["params"]["profiling"]["bgi_soap"]["report_repeat_hits"],
            match_model = config["params"]["profiling"]["bgi_soap"]["match_model"],
            align_seed = config["params"]["profiling"]["bgi_soap"]["align_seed"],
            max_mismatch_num = config["params"]["profiling"]["bgi_soap"]["max_mismatch_num"],
            identity = config["params"]["profiling"]["bgi_soap"]["identity"],
            soap = os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap/{sample}/{sample}.bgi_soap.soap")
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        run:
            if IS_PE:
                shell(
                    '''
                    soap2.22 -a {input.reads[0]} -b {input.reads[1]} -D {params.index} \
                    -m {params.minimal_insert_size} \
                    -x {params.maximal_insert_size} \
                    -r {params.report_repeat_hits} \
                    -l {params.align_seed} \
                    -M {params.match_model} \
                    -v {params.max_mismatch_num} \
                    -c {params.identity} \
                    -S -p {threads} \
                    -o {params.soap}.pe \
                    -2 {params.soap}.se \
                    2> {log}
                    ''')
                shell(
                    '''
                    cat {params.soap}.pe {params.soap}.se > {params.soap}
                    rm -rf {params.soap}.pe {params.soap}.se
                    ''')
            else:
                shell(
                    '''
                    soap2.22 -a {input.reads[0]} -D {params.index} \
                    -m {params.minimal_insert_size} \
                    -x {params.maximal_insert_size} \
                    -r {params.report_repeat_hits} \
                    -l {params.align_seed} \
                    -M {params.match_model} \
                    -v {params.max_mismatch_num} \
                    -c {params.identity} \
                    -S -p {threads} \
                    -o {params.soap} \
                    2> {log}
                    ''')
            shell('''pigz {params.soap}''')


    rule profiling_custom_bgi_soap_merge:
        input:
            soap_list = expand(os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap/{sample}/{sample}.bgi_soap.soap.gz"),
                   sample=SAMPLES_ID_LIST),
            taxonomy = config["params"]["profiling"]["bgi_soap"]["taxonomy"]
        output:
            abun_profile = os.path.join(
                config["output"]["profiling"],
                "report/bgi_soap/bgi_soap.merged.abundance.profile.tsv"),
            count_profile = os.path.join(
                config["output"]["profiling"],
                "report/bgi_soap/bgi_soap.merged.count.profile.tsv")
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/bgi_soap/profiling_bgi_soap_merge.log")
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        run:
            quantpi.profiler_init(input.taxonomy)

            count_df, abun_df = quantpi.get_all_abun_df(input.soap_list, threads, "bgi_soap")
            count_df.to_csv(output.count_profile, sep='\t', index=False)
            abun_df.to_csv(output.abun_profile, sep='\t', index=False)


    rule profiling_custom_bgi_soap_all:
        input:
            os.path.join(
                config["output"]["profiling"],
                "report/bgi_soap/bgi_soap.merged.abundance.profile.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "report/bgi_soap/bgi_soap.merged.count.profile.tsv"),

            rules.qcreport_all.input

else:
    rule profiling_custom_bgi_soap_all:
        input:


if config["params"]["profiling"]["bowtie2"]["do"]:
    rule profiling_custom_bowtie2:
        input:
            reads = profiling_input_with_short_reads,
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["profiling"]["bowtie2"]["index_prefix"],
                           suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]),
            taxonomy = config["params"]["profiling"]["bowtie2"]["taxonomy"]
        output:
            flagstat = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.flagstat"),
            bam = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.sorted.bam"),
            bai = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.sorted.bam.bai"),
            count_profile = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.count.tsv"),
            abun_profile = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.abun.tsv")
        log:
            os.path.join(config["output"]["profiling"], "logs/bowtie2/{sample}.bowtie2.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/bowtie2/{sample}.profiling.bowtie2.benchmark.txt")
        params:
            index_prefix = config["params"]["profiling"]["bowtie2"]["index_prefix"]
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        run:
            shell(
                '''
                bowtie2 \
                -x {params.index_prefix} \
                --end-to-end \
                --very-sensitive \
                --phred33 \
                --threads {threads} \
                --seed 0 \
                --time \
                -k 2 \
                --no-unal \
                --no-discordant \
                %s \
                2> {log} | \
                tee >(samtools flagstat \
                      -@{threads} - \
                      > {output.flagstat}) | \
                samtools sort \
                -m 3G \
                -@{threads} \
                -T {output.bam} \
                -O BAM -o {output.bam} -
                ''' % \
                "-1 {input.reads[0]} -2 {input.reads[1]}" if IS_PE \
                else "-U {input.reads[0]}")

            shell(
                '''
                samtools index -@{threads} {output.bam} {output.bai} 2>> {log}
                ''')

            quantpi.profiler_init(input.taxonomy)
            count_df, abun_df = quantpi.get_abun_df_bowtie2(output.bam)
            count_df.reset_index().to_csv(output.count_profile, sep='\t', index=False)
            abun_df.reset_index().to_csv(output.abun_profile, sep='\t', index=False)


    rule profiling_custom_bowtie2_merge:
        input:
            count_tsv_list = expand(os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.count.tsv"),
                                    sample=SAMPLES_ID_LIST),
            abun_tsv_list = expand(os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.abun.tsv"),
                                   sample=SAMPLES_ID_LIST)
        output:
            abun_profile = os.path.join(
                config["output"]["profiling"],
                "report/bowtie2/bowtie2.merged.abundance.profile.tsv"),
            count_profile = os.path.join(
                config["output"]["profiling"],
                "report/bowtie2/bowtie2.merged.count.profile.tsv")
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/bowtie2/profiling_bowtie2_merge.log")
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        run:
            def parse(tsv_file):
                df = pd.read_csv(tsv_file, sep='\t').set_index("lineages_full")
                return df

            count_list = [parse(i) for i in input.count_tsv_list]
            abun_list = [parse(i) for i in input.abun_tsv_list]

            count_df = pd.concat(count_list, axis=1).reset_index()
            abun_df = pd.concat(abun_list, axis=1).reset_index()

            count_df.to_csv(output.count_profile, sep='\t', index=False)
            abun_df.to_csv(output.abun_profile, sep='\t', index=False)


    rule profiling_custom_bowtie2_all:
        input:
            os.path.join(
                config["output"]["profiling"],
                "report/bowtie2/bowtie2.merged.abundance.profile.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "report/bowtie2/bowtie2.merged.count.profile.tsv"),

            rules.qcreport_all.input

else:
    rule profiling_custom_bowtie2_all:
        input:


if config["params"]["profiling"]["jgi"]["do"]:
    rule profiling_custom_jgi:
        input:
            reads = profiling_input_with_short_reads,
            index_database = expand(
                "{prefix}.{suffix}",
                prefix=config["params"]["profiling"]["jgi"]["index_prefix"],
                suffix=["1.bt2l", "2.bt2l", "3.bt2l", "4.bt2l",
                        "rev.1.bt2l", "rev.2.bt2l"])
        output:
            os.path.join(
                config["output"]["profiling"],
                "profile/jgi/{sample}/{sample}.jgi.coverage.gz")
        conda:
            config["envs"]["jgi"]
        log:
            os.path.join(config["output"]["profiling"], "logs/jgi/{sample}.jgi.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/jgi/{sample}.jgi.benchmark.txt")
        threads:
            config["params"]["profiling"]["threads"]
        params:
            index_prefix = config["params"]["profiling"]["jgi"]["index_prefix"],
            memory_limit = config["params"]["profiling"]["jgi"]["memory_limit"],
            fragment = config["params"]["profiling"]["jgi"]["fragment"],
            tmp_prefix = os.path.join(
                config["output"]["profiling"],
                "profile/jgi/{sample}/{sample}.temp"),
            coverage = os.path.join(
                config["output"]["profiling"],
                "profile/jgi/{sample}/{sample}.jgi.coverage")
        priority:
            20
        shell:
            '''
            reads=""
            if [ IS_PE ]
            then
                reads="-1 {input.reads[0]} -2 {input.reads[1]}"
            else
                reads="-U {input.reads[0]"
            fi


            bowtie2 \
            -x {params.index_prefix} \
            $reads \
            --end-to-end \
            --very-sensitive \
            --phred33 \
            --threads {threads} \
            --seed 0 \
            --time \
            -k 2 \
            --no-unal \
            --no-discordant \
            -X {params.fragment} \
            2> {log} | \
            samtools sort \
            -@{threads} \
            -m {params.memory_limit} \
            -T {params.tmp_prefix} -O BAM - | \
            jgi_summarize_bam_contig_depths \
            --outputDepth {params.coverage} -

            gzip {params.coverage}
            rm -rf {params.temp_prefix}*.bam
            echo "Profiling done" >> {log}
            '''


    rule profiling_custom_jgi_merge:
        input:
            coverage = expand(os.path.join(
                config["output"]["profiling"],
                "profile/jgi/{sample}/{sample}.jgi.coverage.gz"),
                              sample=SAMPLES_ID_LIST),
            taxonomy = config["params"]["profiling"]["jgi"]["taxonomy"],
            index_metadata = config["params"]["profiling"]["jgi"]["metadata"]
        output:
            abundance_profile = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.abundance.profile.tsv"),
            coverage_profile = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.coverage.profile.tsv"),
            abundance_profile_k = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.abundance.profile.superkingdom.tsv"),
            abundance_profile_p = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.abundance.profile.phylum.tsv"),
            abundance_profile_o = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.abundance.profile.order.tsv"),
            abundance_profile_c = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.abundance.profile.class.tsv"),
            abundance_profile_f = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.abundance.profile.family.tsv"),
            abundance_profile_g = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.abundance.profile.genus.tsv"),
            abundance_profile_s = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.abundance.profile.species.tsv"),
            abundance_profile_t = os.path.join(
                config["output"]["profiling"],
                "report/jgi/jgi.merged.abundance.profile.strain.tsv")
        threads:
            config["params"]["profiling"]["threads"]
        priority:
            20
        run:
            taxonomy_df = pd.read_csv(input.taxonomy, sep='\t')

            quantpi.profiler_init(input.index_metadata)

            coverage_df, abundance_df = quantpi.get_all_abun_df(
                input.coverage, threads, "jgi")

            samples_list = sorted(abundance_df.columns[1:].to_list())

            coverage_df.to_csv(output.coverage_profile, sep='\t', index=False)
            abundance_df.to_csv(output.abundance_profile, sep='\t', index=False)

            abundance_taxonomy_df = abundance_df.merge(taxonomy_df)

            quantpi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_superkingdom_new", output.abundance_profile_k)

            quantpi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_phylum_new", output.abundance_profile_p)

            quantpi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_order_new", output.abundance_profile_o)

            quantpi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_class_new", output.abundance_profile_c)

            quantpi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_family_new", output.abundance_profile_f)

            quantpi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_genus_new", output.abundance_profile_g)

            quantpi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_species_new", output.abundance_profile_s)

            quantpi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_strain_new", output.abundance_profile_t)


    rule profiling_custom_jgi_all:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "report/jgi/jgi.merged.{target}.profile.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/jgi/jgi.merged.abundance.profile.{level}.tsv")],
                   target=["abundance", "coverage"],
                   level=["superkingdom",
                          "phylum",
                          "order",
                          "class",
                          "family",
                          "genus",
                          "species",
                          "strain"]),

            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_custom_jgi_all:
        input: