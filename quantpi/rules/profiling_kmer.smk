if len(KMCP_DBS) > 0:
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
            config["params"]["profiling"]["threads"]
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
        benchmark:
            os.path.join(config["output"]["profiling"],
                "benchmark/kmcp/search_merge/{sample}.kmcp_search_merge.benchmark.txt")
        shell:
            '''
            kmcp merge \
            {input} \
            --out-file {output} \
            --log {log}
            '''


    KMCP_PROFILING_MODE = {
        "pathogen_detection": 0,
        "higher_recall": 1,
        "high_recall": 2,
        "default": 3,
        "high_precision": 4,
        "higher_precision": 5
        }


    rule profiling_kmcp_profile:
        input:
            search = os.path.join(config["output"]["profiling"],
                "search/kmcp/{sample}/{sample}.kmcp_search@all.tsv.gz"),
            taxid = KMCP_TAXIDMAP,
            taxdump = config["params"]["profiling"]["kmcp"]["database"]["ncbi_taxdump"]
        output:
            default_profile = os.path.join(config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.default_format.profile.{profiling_mode}.tsv.gz"),
            metaphlan_profile = os.path.join(config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.metaphlan_format.profile.{profiling_mode}.tsv.gz"),
            cami_profile = os.path.join(config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.CAMI_format.profile.{profiling_mode}.tsv.gz"),
            binning_result = os.path.join(config["output"]["profiling"],
                "profile/kmcp/{sample}/{sample}.kmcp.binning.{profiling_mode}.gz")
        conda:
            config["envs"]["kmcp"]
        log:
            os.path.join(config["output"]["profiling"],
                "logs/kmcp/profile/{sample}.kmcp_profile_{profiling_mode}.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                "benchmark/kmcp/profile/{sample}.kmcp_profile_{profiling_mode}.benchmark.txt")
        params:
            profiling_mode = lambda wildcards: KMCP_PROFILING_MODE[wildcards.profiling_mode],
            min_query_cov = config["params"]["profiling"]["kmcp"]["profile"]["min_query_cov"],
            min_chunks_reads = config["params"]["profiling"]["kmcp"]["profile"]["min_chunks_reads"],
            min_chunks_fraction = config["params"]["profiling"]["kmcp"]["profile"]["min_chunks_fraction"],
            max_chunks_depth_stdev = config["params"]["profiling"]["kmcp"]["profile"]["max_chunks_depth_stdev"],
            min_uniq_reads = config["params"]["profiling"]["kmcp"]["profile"]["min_uniq_reads"],
            min_hic_ureads = config["params"]["profiling"]["kmcp"]["profile"]["min_hic_ureads"],
            min_hic_ureads_qcov = config["params"]["profiling"]["kmcp"]["profile"]["min_hic_ureads_qcov"],
            min_hic_ureads_prop = config["params"]["profiling"]["kmcp"]["profile"]["min_hic_ureads_prop"],
            external_opts = config["params"]["profiling"]["kmcp"]["profile"]["external_opts"]
        shell:
            '''
            kmcp profile \
            {input.search} \
            --taxid-map {input.taxid} \
            --taxdump {input.taxdump} \
            --level species \
            --mode {params.profiling_mode} \
            --min-query-cov {params.min_query_cov} \
            --min-chunks-reads {params.min_chunks_reads} \
            --min-chunks-fraction {params.min_chunks_fraction} \
            --max-chunks-depth-stdev {params.max_chunks_depth_stdev} \
            --min-uniq-reads {params.min_uniq_reads} \
            --min-hic-ureads {params.min_hic_ureads} \
            --min-hic-ureads-qcov {params.min_hic_ureads_qcov} \
            --min-hic-ureads-prop {params.min_hic_ureads_prop} \
            {params.external_opts} \
            --out-prefix {output.default_profile} \
            --metaphlan-report {output.metaphlan_profile} \
            --cami-report {output.cami_profile} \
            --binning-result {output.binning_result} \
            --log {log}
            '''


    rule profiling_kmcp_all:
        input:
            expand([
                os.path.join(config["output"]["profiling"],
                    "profile/kmcp/{sample}/{sample}.kmcp.{profile_format}.profile.{profiling_mode}.tsv.gz"),
                os.path.join(config["output"]["profiling"],
                    "profile/kmcp/{sample}/{sample}.kmcp.binning.{profiling_mode}.gz")],
                sample=SAMPLES_ID_LIST,
                profile_format=["default_format", "metaphlan_format", "CAMI_format"],
                profiling_mode=list(KMCP_PROFILING_MODE.keys()))
 
else:
    rule profiling_kmcp_all:
        input: