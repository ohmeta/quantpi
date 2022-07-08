KMCP_DB_NUMBER = len(KMCP_DBS)

if KMCP_DB_NUMBER > 0:
    rule profiling_kmcp_search:
        input:
            reads = profiling_input_with_short_reads,
            db_dir = lambda wildcards: config["params"]["profiling"]["kmcp"]["database"][wildcards.kmcp_db]
        output:
            os.path.join(config["output"]["profiling"],
                "search/kmcp/{sample}/{sample}.kmcp_search@{kmcp_db}.tsv.gz")
        conda:
            config["envs"]["kmcp"]
        log:
            os.path.join(config["output"]["profiling"],
                "logs/kmcp/search/{sample}.kmcp_search@{kmcp_db}.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                "benchmark/kmcp/search/{sample}.kmcp_search@{kmcp_db}.benchmark.txt")
        params:
            reads_mode = config["params"]["profiling"]["kmcp"]["search"]["reads_mode"],
            reads_layout = 1 if IS_PE else 0,
            min_query_len = config["params"]["profiling"]["kmcp"]["search"]["min_query_len"],
            min_query_cov = config["params"]["profiling"]["kmcp"]["search"]["min_query_cov"],
            external_opts = config["params"]["profiling"]["kmcp"]["search"]["external_opts"]
        priority:
            20
        threads:
            config["params"]["profiling"]["kmcp"]["search"]["threads"]
        shell:
            '''
            if [ {params.reads_layout} -eq 1 ] && [ "{params.reads_mode}" == "Paired-end" ]
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


    rule profiling_kmcp_search_merge:
        input:
            expand(os.path.join(config["output"]["profiling"],
                "search/kmcp/{{sample}}/{{sample}}.kmcp_search@{kmcp_db}.tsv.gz"),
                kmcp_db=KMCP_DBS)
        output:
            os.path.join(config["output"]["profiling"],
                "search/kmcp/{sample}/{sample}.kmcp_search@all.tsv.gz")
        conda:
            config["envs"]["kmcp"]
        log:
            os.path.join(config["output"]["profiling"],
                "logs/kmcp/search_merge/{sample}.kmcp_search_merge.log")
        priority:
            21
        benchmark:
            os.path.join(config["output"]["profiling"],
                "benchmark/kmcp/search_merge/{sample}.kmcp_search_merge.benchmark.txt")
        params:
            kmcp_db_number = KMCP_DB_NUMBER,
            outdir = os.path.join(config["output"]["profiling"], "search/kmcp/{sample}")
        threads:
            config["params"]["profiling"]["kmcp"]["profile"]["threads"]
        shell:
            '''
            rm -rf {output}

            if [ {params.kmcp_db_number} == 1 ]
            then
                pushd {params.outdir} && \
                ln -s $(basename {input[0]}) $(basename {output}) && \
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
        shell:
            '''
            taxidmap=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.taxidmap})

            kmcp profile \
            --threads {threads} \
            --taxid-map $taxidmap \
            --taxdump {input.taxdump} \
            --level species \
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


    rule profiling_kmcp_profile_merge:
        input:
           abuns = expand(os.path.join(
                config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.metaphlan_format.{{profiling_mode}}.profile"),
                          sample=SAMPLES_ID_LIST)
        output:
            profiles = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/kmcp.metaphlan_format.{{profiling_mode}}.merged.abundance.profile.{level}.tsv"),
                level=["all", "strain", "species", "genus", "family",
                       "order", "class", "phylum", "superkingdom"])
        log:
            os.path.join(config["output"]["profiling"],
                "logs/kmcp/profile_merge/kmcp_profile_merge_{profiling_mode}.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                "benchmark/kmcp/profile_merge/kmcp_profile_merge_{profiling_mode}.benchmark.txt")
        threads:
            config["params"]["profiling"]["threads"]
        priority:
            24
        params:
            metaphlan_report_version = config["params"]["profiling"]["kmcp"]["profile"]["metaphlan_report_version"]
        run:
            quantpi.metaphlan_init(int(params.metaphlan_report_version))
            profile_list = quantpi.merge_metaphlan_tables(input.abuns, threads)
            for i in range(0, len(profile_list)):
                profile_list[i].to_csv(output.profiles[i], sep='\t', index=False)


    rule profiling_kmcp_all:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/kmcp/kmcp.metaphlan_format.{profiling_mode}.merged.abundance.profile.{level}.tsv"),
                level=["all", "strain", "species", "genus", "family",
                       "order", "class", "phylum", "superkingdom"],
                profiling_mode=list(KMCP_PROFILING_MODE_DO.keys()))
 
else:
    rule profiling_kmcp_all:
        input: