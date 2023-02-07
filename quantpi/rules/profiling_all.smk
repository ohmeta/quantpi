rule profiling_all:
    input:
        rules.profiling_kraken2_all.input,
        rules.profiling_bracken_all.input,
        rules.profiling_metaphlan2_all.input,
        rules.profiling_metaphlan3_all.input,
        rules.profiling_metaphlan4_all.input,
        rules.profiling_strainphlan3_all.input,
        rules.profiling_strainphlan4_all.input,
        rules.profiling_kmcp_all.input,
        #rules.profiling_custom_bgi_soap_all.input,
        #rules.profiling_custom_bowtie2_all.input,
        #rules.profiling_custom_jgi_all.input,
        rules.profiling_genomecov_all.input,
        rules.profiling_genome_coverm_all.input,
        rules.profiling_humann2_all.input,
        rules.profiling_humann3_all.input,
        rules.profiling_humann4_all.input


localrules:
    profiling_kraken2_all,
    profiling_bracken_all,
    profiling_metaphlan2_all,
    profiling_metaphlan3_all,
    profiling_metaphlan4_all,
    profiling_strainphlan3_all,
    profiling_strainphlan4_all,
    profiling_kmcp_all,
    profiling_genomecov_all,
    profiling_genome_coverm_all,
    #profiling_custom_bgi_soap_all,
    #profiling_custom_bowtie2_all,
    #profiling_custom_jgi_all,
    profiling_humann2_all,
    profiling_humann3_all,
    profiling_humann4_all,
    profiling_all