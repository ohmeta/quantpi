KMCP_DB_NUMBER = len(KMCP_DBS)

if KMCP_DB_NUMBER > 0:
    if IS_PE:
        rule profiling_kmcp_search:
            input:
                reads = profiling_input_with_short_reads,
                db_dir = lambda wildcards: config["params"]["profiling"]["kmcp"]["database"][wildcards.kmcp_db]
            output:
                os.path.join(
                    config["output"]["profiling"], "search/kmcp/{sample}/{sample}.kmcp_search@{kmcp_db}.tsv.gz")
            conda:
                config["envs"]["kmcp"]
            log:
                os.path.join(
                    config["output"]["profiling"],
                    "logs/kmcp/search/{sample}.kmcp_search@{kmcp_db}.log")
            benchmark:
                os.path.join(
                    config["output"]["profiling"],
                    "benchmark/kmcp/search/{sample}.kmcp_search@{kmcp_db}.benchmark.txt")
            params:
                reads_mode = config["params"]["profiling"]["kmcp"]["search"]["reads_mode"],
                min_query_len = config["params"]["profiling"]["kmcp"]["search"]["min_query_len"],
                min_query_cov = config["params"]["profiling"]["kmcp"]["search"]["min_query_cov"],
                external_opts = config["params"]["profiling"]["kmcp"]["search"]["external_opts"]
            priority:
                20
            threads:
                config["params"]["profiling"]["kmcp"]["search"]["threads"]
            resources:
                mem_mb=config["params"]["profiling"]["kmcp"]["mem_mb"]
            shell:
                '''
                if [ "{params.reads_mode}" == "Paired-end" ]
                then
                    kmcp search \
                    --threads {threads} \
                    --load-whole-db \
                    --db-dir {input.db_dir} \
                    --min-query-len {params.min_query_len} \
                    --min-query-cov {params.min_query_cov} \
                    {params.external_opts} \
                    -1 {input.reads[0]} \
                    -2 {input.reads[1]} \
                    --out-file {output} \
                    --log {log}
                else
                    kmcp search \
                    --threads {threads} \
                    --load-whole-db \
                    --db-dir {input.db_dir} \
                    --min-query-len {params.min_query_len} \
                    --min-query-cov {params.min_query_cov} \
                    {input.reads} \
                    --out-file {output} \
                    --log {log}
                fi
                '''
    else:
        rule profiling_kmcp_search:
            input:
                reads = profiling_input_with_short_reads,
                db_dir = lambda wildcards: config["params"]["profiling"]["kmcp"]["database"][wildcards.kmcp_db]
            output:
                temp(os.path.join(config["output"]["profiling"], "search/kmcp/{sample}/{sample}.kmcp_search@{kmcp_db}.tsv.gz"))
            conda:
                config["envs"]["kmcp"]
            log:
                os.path.join(config["output"]["profiling"], "logs/kmcp/search/{sample}.kmcp_search@{kmcp_db}.log")
            benchmark:
                os.path.join(config["output"]["profiling"], "benchmark/kmcp/search/{sample}.kmcp_search@{kmcp_db}.benchmark.txt")
            params:
                reads_mode = config["params"]["profiling"]["kmcp"]["search"]["reads_mode"],
                min_query_len = config["params"]["profiling"]["kmcp"]["search"]["min_query_len"],
                min_query_cov = config["params"]["profiling"]["kmcp"]["search"]["min_query_cov"],
                external_opts = config["params"]["profiling"]["kmcp"]["search"]["external_opts"]
            priority:
                20
            threads:
                config["params"]["profiling"]["kmcp"]["search"]["threads"]
            resources:
                mem_mb=config["params"]["profiling"]["kmcp"]["mem_mb"]
            shell:
                '''
                kmcp search \
                --threads {threads} \
                --load-whole-db \
                --db-dir {input.db_dir} \
                --min-query-len {params.min_query_len} \
                --min-query-cov {params.min_query_cov} \
                {input.reads} \
                --out-file {output} \
                --log {log}
                '''


    rule profiling_kmcp_search_merge:
        input:
            expand(os.path.join(config["output"]["profiling"],
                "search/kmcp/{{sample}}/{{sample}}.kmcp_search@{kmcp_db}.tsv.gz"),
                kmcp_db=KMCP_DBS)
        output:
            os.path.join(config["output"]["profiling"], "search/kmcp/{sample}/{sample}.kmcp_search@all.tsv.gz")
        conda:
            config["envs"]["kmcp"]
        log:
            os.path.join(config["output"]["profiling"],
                "logs/kmcp/search_merge/{sample}.kmcp_search_merge.log")
        priority:
            21
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/kmcp/search_merge/{sample}.kmcp_search_merge.benchmark.txt")
        params:
            kmcp_db_number = KMCP_DB_NUMBER,
            outdir = os.path.join(config["output"]["profiling"], "search/kmcp/{sample}")
        threads:
            config["params"]["profiling"]["kmcp"]["profile"]["threads"]
        resources:
            mem_mb=config["params"]["profiling"]["kmcp"]["mem_mb"]
        shell:
            '''
            rm -rf {output}

            if [ {params.kmcp_db_number} == 1 ]
            then
                pushd {params.outdir} && \
                mv $(basename {input[0]}) $(basename {output}) && \
                popd
            else
                kmcp merge \
                --threads {threads} \
                --out-file {output} \
                --log {log} \
                {input}
            fi
            '''


    KMCP_PROFILING_MODE = {
        0: "pathogen_detection",
        1: "higher_recall",
        2: "high_recall",
        3: "default",
        4: "high_precision",
        5: "higher_precision"
    }


    KMCP_PROFILING_MODE_DO = {}
    for profile_mode in config["params"]["profiling"]["kmcp"]["profile"]["mode"]:
        if profile_mode in KMCP_PROFILING_MODE:
            KMCP_PROFILING_MODE_DO[KMCP_PROFILING_MODE[profile_mode]] = profile_mode
        else:
            print(f'''profile mode {profile_mode} is not supported by KMCP''')
            sys.exit(1)


    rule profiling_kmcp_profile:
        input:
            search = os.path.join(config["output"]["profiling"],
                "search/kmcp/{sample}/{sample}.kmcp_search@all.tsv.gz"),
            taxidmap= KMCP_TAXIDMAP,
            taxdump = config["params"]["profiling"]["kmcp"]["database"]["taxdump"]
        output:
            default_profile = os.path.join(config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.default_format.{profiling_mode}.profile"),
            metaphlan_profile = os.path.join(config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.metaphlan_format.{profiling_mode}.profile"),
            cami_profile = os.path.join(config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.CAMI_format.{profiling_mode}.profile"),
            binning_result = os.path.join(config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.{profiling_mode}.binning.gz")
        priority:
            22
        conda:
            config["envs"]["kmcp"]
        log:
            os.path.join(config["output"]["profiling"],
                "logs/kmcp/profile/{sample}.kmcp_profile_{profiling_mode}.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                "benchmark/kmcp/profile/{sample}.kmcp_profile_{profiling_mode}.benchmark.txt")
        params:
            sample_id = "{sample}",
            level = config["params"]["profiling"]["kmcp"]["profile"]["level"],
            binning_result = os.path.join(config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.{profiling_mode}"),
            profiling_mode = lambda wildcards: KMCP_PROFILING_MODE_DO[wildcards.profiling_mode],
            metaphlan_report_version = config["params"]["profiling"]["kmcp"]["profile"]["metaphlan_report_version"],
            disable_two_stage_taxonomy_assignment = "--no-amb-corr" \
                if config["params"]["profiling"]["kmcp"]["profile"]["disable_two_stage_taxonomy_assignment"] \
                else "",
            external_opts = config["params"]["profiling"]["kmcp"]["profile"]["external_opts"]
        threads:
            config["params"]["profiling"]["kmcp"]["profile"]["threads"]
        resources:
            mem_mb=config["params"]["profiling"]["kmcp"]["mem_mb"]
        shell:
            '''
            taxidmap=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.taxidmap})

            kmcp profile \
            --threads {threads} \
            --taxid-map $taxidmap \
            --taxdump {input.taxdump} \
            --level {params.level} \
            --mode {params.profiling_mode} \
            --out-prefix {output.default_profile} \
            --metaphlan-report {output.metaphlan_profile} \
            --metaphlan-report-version {params.metaphlan_report_version} \
            --cami-report {output.cami_profile} \
            --sample-id {params.sample_id} \
            --binning-result {params.binning_result} \
            --log {log} \
            {params.disable_two_stage_taxonomy_assignment} \
            {params.external_opts} \
            {input.search}
            '''


    rule profiling_kmcp_profile_merge_kmcp:
        input:
            profile = expand(os.path.join(
                config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.default_format.{{profiling_mode}}.profile"),
                sample=SAMPLES_ID_LIST)
        output:
            percentage = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/kmcp_format/{{profiling_mode}}/merged.percentage.profile.{level}.tsv"),
                level=["species"]),
            coverage = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/kmcp_format/{{profiling_mode}}/merged.coverage.profile.{level}.tsv"),
                level=["species"]),
            reads = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/kmcp_format/{{profiling_mode}}/merged.reads.profile.{level}.tsv"),
                level=["species"]),
            ureads = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/kmcp_format/{{profiling_mode}}/merged.ureads.profile.{level}.tsv"),
                level=["species"]),
            hicureads = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/kmcp_format/{{profiling_mode}}/merged.hicureads.profile.{level}.tsv"),
                level=["species"])
        log:
            os.path.join(config["output"]["profiling"],
                "logs/kmcp/profile_merge/kmcp_profile_merge_kmcp_{profiling_mode}.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                "benchmark/kmcp/profile_merge/kmcp_profile_merge_kmcp_{profiling_mode}.benchmark.txt")
        threads:
            config["params"]["profiling"]["threads"]
        resources:
            mem_mb=config["params"]["profiling"]["kmcp"]["mem_mb"]
        priority:
            24
        run:
            #output_dict = {
            #    "percentage_t": output.percentage[0],
            #    "percentage_s": output.percentage[1],
            #    "coverage_t": output.coverage[0],
            #    "coverage_s": output.coverage[1],
            #    "reads_t": output.reads[0],
            #    "reads_s": output.reads[1],
            #    "ureads_t": output.ureads[0],
            #    "ureads_s": output.ureads[1],
            #    "hicureads_t": output.hicureads[0],
            #    "hicureads_s": output.hicureads[1]
            #}
            output_dict = {
                "percentage_s": output.percentage[0],
                "coverage_s": output.coverage[0],
                "reads_s": output.reads[0],
                "ureads_s": output.ureads[0],
                "hicureads_s": output.hicureads[0],
            }
            quantpi.kmcp_profile_merge_species(input.profile, threads, output_dict)


    rule profiling_kmcp_profile_merge:
        input:
            profile = expand(os.path.join(
                config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.metaphlan_format.{{profiling_mode}}.profile"),
                sample=SAMPLES_ID_LIST)
        output:
            profiles = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/metaphlan_format/{{profiling_mode}}/merged.abundance.profile.{level}.tsv"),
                level=["all", "strain", "species", "genus", "family", "order", "class", "phylum", "superkingdom"])
        log:
            os.path.join(config["output"]["profiling"],
                "logs/kmcp/profile_merge/kmcp_profile_merge_{profiling_mode}.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                "benchmark/kmcp/profile_merge/kmcp_profile_merge_{profiling_mode}.benchmark.txt")
        threads:
            config["params"]["profiling"]["threads"]
        resources:
            mem_mb=config["params"]["profiling"]["kmcp"]["mem_mb"]
        priority:
            24
        params:
            metaphlan_report_version = config["params"]["profiling"]["kmcp"]["profile"]["metaphlan_report_version"]
        run:
            quantpi.metaphlan_init(int(params.metaphlan_report_version))
            profile_list = quantpi.merge_metaphlan_tables(input.profile, threads)
            for i in range(0, len(profile_list)):
                profile_list[i].to_csv(output.profiles[i], sep='\t', index=False)


    rule profiling_kmcp_all:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/metaphlan_format/{profiling_mode}/merged.abundance.profile.{level}.tsv"),
                level=["all", "strain", "species", "genus", "family", "order", "class", "phylum", "superkingdom"],
                profiling_mode=list(KMCP_PROFILING_MODE_DO.keys())),
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/kmcp_format/{profiling_mode}/merged.{target}.profile.{level}.tsv"),
                target=["percentage", "coverage", "reads", "ureads", "hicureads"],
                #level=["strain", "species"],
                level=["species"],
                profiling_mode=list(KMCP_PROFILING_MODE_DO.keys()))

else:
    rule profiling_kmcp_all:
        input: