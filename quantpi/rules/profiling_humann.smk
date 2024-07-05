# HUMAnN2 v2.8
## utility_mapping : full = http://huttenhower.sph.harvard.edu/humann2_data/full_mapping_1_1.tar.gz
## chocophlan : DEMO = http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/DEMO_chocophlan.v0.1.1.tar.gz
## chocophlan : full = http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan_plus_viral.v0.1.1.tar.gz
## uniref : DEMO_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref90_DEMO_diamond.tar.gz
## uniref : uniref90_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref90_annotated_1_1.tar.gz
## uniref : uniref50_ec_filtered_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_ec_filtered/uniref50_ec_filtered_1_1.tar.gz
## uniref : uniref50_GO_filtered_rapsearch2 = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref50_GO_filtered/uniref50_GO_filtered_rapsearch2.tar.gz
## uniref : uniref50_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref50_annotated_1_1.tar.gz
## uniref : uniref90_ec_filtered_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_ec_filtered/uniref90_ec_filtered_1_1.tar.gz

# HUMAnN v3.0
## chocophlan : full = http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v296_201901b.tar.gz
## chocophlan : DEMO = http://huttenhower.sph.harvard.edu/humann_data/chocophlan/DEMO_chocophlan.v296_201901b.tar.gz
## uniref : uniref50_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref50_annotated_v201901b_full.tar.gz
## uniref : uniref90_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz
## uniref : uniref50_ec_filtered_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_ec_filtered/uniref50_ec_filtered_201901b_subset.tar.gz
## uniref : uniref90_ec_filtered_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_ec_filtered/uniref90_ec_filtered_201901b_subset.tar.gz
## uniref : DEMO_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_DEMO_diamond_v201901b.tar.gz
## utility_mapping : full = http://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz

# HUMAnN v3.1, v3.5, v3.6, v3.7, v3.8, v3.9
## chocophlan : full = http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz
## chocophlan : DEMO = http://huttenhower.sph.harvard.edu/humann_data/chocophlan/DEMO_chocophlan.v201901_v31.tar.gz
## uniref : uniref50_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref50_annotated_v201901b_full.tar.gz
## uniref : uniref90_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz
## uniref : uniref50_ec_filtered_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_ec_filtered/uniref50_ec_filtered_201901b_subset.tar.gz
## uniref : uniref90_ec_filtered_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_ec_filtered/uniref90_ec_filtered_201901b_subset.tar.gz
## uniref : DEMO_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_DEMO_diamond_v201901b.tar.gz
## utility_mapping : full = http://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz


rule profiling_humann2_config:
    output:
        touch(os.path.join(config["output"]["profiling"], "config/humann2.config.done"))
    log:
        os.path.join(config["output"]["profiling"], "logs/humann2_config/humann2.config.log")
    conda:
        config["envs"]["biobakery2"]
    params:
        database_utility_mapping = config["params"]["profiling"]["humann"]["database_utility_mapping_v2"],
        database_nucleotide = config["params"]["profiling"]["humann"]["database_nucleotide_v2"],
        database_protein = config["params"]["profiling"]["humann"]["database_protein_v2"],
        threads = config["params"]["profiling"]["threads"]
    priority:
        20
    shell:
        '''
        humann2_config > {log}
        humann2_config --update database_folders utility_mapping {params.database_utility_mapping}
        humann2_config --update database_folders nucleotide {params.database_nucleotide}
        humann2_config --update database_folders protein {params.database_protein}
        humann2_config --update run_modes threads {params.threads}
        echo "####" >> {log}
        humann2_config >> {log}
        '''

localrules:
    profiling_humann2_config


rule profiling_humann2_build_chocophlan_pangenome_db:
    input:
        tag = os.path.join(config["output"]["profiling"], "config/humann2.config.done"),
        profile = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv")
    output:
        expand(os.path.join(
            config["output"]["profiling"],
            "database/humann2/{{sample}}/{{sample}}_bowtie2_index.{suffix}"),
            suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    conda:
        config["envs"]["biobakery2"]
    log:
        os.path.join(config["output"]["profiling"], "logs/humann2/{sample}.humann2.build_pandb.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/humann2/{sample}.bowtie2_index.benchmark.txt")
    params:
        basename = "{sample}",
        wrapper_dir = WRAPPER_DIR,
        db_dir = os.path.join(config["output"]["profiling"], "database/humann2/{sample}"),
        prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"]
    priority:
        20
    shell:
        '''
        python {params.wrapper_dir}/humann2_db_wrapper.py \
        --log {log} \
        --basename {params.basename} \
        --db_dir {params.db_dir} \
        --prescreen_threshold {params.prescreen_threshold} \
        --taxonomic_profile {input.profile}
        '''


rule profiling_humann2:
    input:
        tag = os.path.join(config["output"]["profiling"], "config/humann2.config.done"),
        reads = profiling_input_with_short_reads,
        index = expand(os.path.join(
            config["output"]["profiling"],
            "database/humann2/{{sample}}/{{sample}}_bowtie2_index.{suffix}"),
            suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    output:
        genefamilies = os.path.join(
            config["output"]["profiling"],
            "profile/humann2/{sample}/{sample}_genefamilies.tsv"),
        pathabundance = os.path.join(
            config["output"]["profiling"],
            "profile/humann2/{sample}/{sample}_pathabundance.tsv"),
        pathcoverage = os.path.join(
            config["output"]["profiling"],
            "profile/humann2/{sample}/{sample}_pathcoverage.tsv")
    log:
        os.path.join(config["output"]["profiling"], "logs/{sample}.humann2.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/humann2/{sample}.humann2.benchmark.txt")
    conda:
        config["envs"]["biobakery2"]
    params:
        basename = "{sample}",
        output_dir = os.path.join(config["output"]["profiling"], "profile/humann2/{sample}"),
        index = os.path.join(config["output"]["profiling"], "database/humann2/{sample}/{sample}_bowtie2_index"),
        prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"],
        identity_threshold = config["params"]["profiling"]["humann"]["identity_threshold"],
        translated_subject_coverage_threshold = config["params"]["profiling"]["humann"]["translated_subject_coverage_threshold"],
        translated_query_coverage_threshold = config["params"]["profiling"]["humann"]["translated_query_coverage_threshold"],
        external_opts = config["params"]["profiling"]["humann"]["external_opts_v2"]
    priority:
        20
    threads:
        config["params"]["profiling"]["threads"]
    shell:
        '''
        zcat {input.reads} | \
        bowtie2 \
        --threads {threads} \
        -x {params.index} \
        -U - 2>> {log} | \
        humann2 \
        --threads {threads} \
        --input - \
        --input-format sam \
        --prescreen-threshold {params.prescreen_threshold} \
        --identity-threshold {params.identity_threshold} \
        --translated-subject-coverage-threshold {params.translated_subject_coverage_threshold} \
        --translated-query-coverage-threshold {params.translated_query_coverage_threshold} \
        --output-basename {params.basename} \
        --output {params.output_dir} \
        --remove-temp-output \
        {params.external_otps} \
        --o-log {log}
        '''


rule profiling_humann2_postprocess:
    input:
        expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann2/{{sample}}/{{sample}}_{target}.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"])
    output:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann2/{{sample}}/{{sample}}_{target}_{norm}.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"],
            norm=config["params"]["profiling"]["humann"]["normalize_method"]),
        groupprofiles = expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann2/{{sample}}/{{sample}}_{group}_groupped.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann2/{sample}.humann2_postprocess.log")
    conda:
        config["envs"]["biobakery2"]
    params:
        wrapper_dir =WRAPPER_DIR,
        normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
        regroup_method = config["params"]["profiling"]["humann"]["regroup_method"],
        map_database =  config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        humann2_renorm_table \
        --input {input[0]} \
        --update-snames \
        --output {output.targets[0]} \
        --units {params.normalize_method} \
        > {log} 2>&1

        humann2_renorm_table \
        --input {input[1]} \
        --update-snames \
        --output {output.targets[1]} \
        --units {params.normalize_method} \
        >> {log} 2>&1

        humann2_renorm_table \
        --input {input[2]} \
        --update-snames \
        --output {output.targets[2]} \
        --units {params.normalize_method} \
        >> {log} 2>&1

        python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
        regroup_table \
        --input {input[0]} \
        --groups {params.map_database} \
        --function {params.regroup_method} \
        --output {output.groupprofiles} \
        >> {log} 2>&1
        '''


rule profiling_humann2_join:
    input:
        expand([
            os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}_{target}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}_{target}_{norm}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}_{group}_groupped.tsv")],
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm = config["params"]["profiling"]["humann"]["normalize_method"],
                group=config["params"]["profiling"]["humann"]["map_database"],
                sample=SAMPLES_ID_LIST)
    output:
        targets = expand(
            os.path.join(
                config["output"]["profiling"],
                "report/humann2/humann2_{target}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        targets_norm = expand(
            os.path.join(
                config["output"]["profiling"],
                "report/humann2/humann2_{target}_{norm}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"],
            norm=config["params"]["profiling"]["humann"]["normalize_method"]),
        groupprofile = expand(
            os.path.join(
                config["output"]["profiling"],
                "report/humann2/humann2.{group}.joined.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann2/humann2_join.log")
    conda:
        config["envs"]["biobakery2"]
    params:
        wrapper_dir =WRAPPER_DIR,
        input_dir = os.path.join(config["output"]["profiling"], "profile/humann2"),
        normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
        map_database = config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.targets} \
        --file_name genefamilies.tsv pathabundance.tsv pathcoverage.tsv \
        > {log} 2>&1

        python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.targets_norm} \
        --file_name \
        genefamilies_{params.normalize_method}.tsv \
        pathabundance_{params.normalize_method}.tsv \
        pathcoverage_{params.normalize_method}.tsv \
        > {log} 2>&1

        python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.groupprofile} \
        --file_name {params.map_database} \
        >> {log} 2>&1
        '''


rule profiling_humann2_split_stratified:
    input:
        targets = expand(
            os.path.join(
                config["output"]["profiling"],
                "report/humann2/humann2_{target}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        targets_norm = expand(
            os.path.join(
                config["output"]["profiling"],
                "report/humann2/humann2_{target}_{norm}_joined.tsv"),
            norm = config["params"]["profiling"]["humann"]["normalize_method"],
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        groupprofile = expand(
            os.path.join(
                config["output"]["profiling"],
                "report/humann2/humann2.{group}.joined.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    output:
        expand([
            os.path.join(
                config["output"]["profiling"],
                "report/humann2/humann2_{target}_joined_{suffix}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "report/humann2/humann2_{target}_{norm}_joined_{suffix}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "report/humann2/humann2.{group}.joined_{suffix}.tsv")],
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm = config["params"]["profiling"]["humann"]["normalize_method"],
                group=config["params"]["profiling"]["humann"]["map_database"],
                suffix=["stratified", "unstratified"])
    log:
        os.path.join(config["output"]["profiling"],
                        "logs/humann2/humann2_split_stratified.log")
    conda:
        config["envs"]["biobakery2"]
    params:
        wrapper_dir = WRAPPER_DIR,
        output_dir = os.path.join(config["output"]["profiling"], "report/humann2"),
        map_database = config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.targets} \
        --output {params.output_dir} \
        > {log} 2>&1

        python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.targets_norm} \
        --output {params.output_dir} \
        > {log} 2>&1

        python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.groupprofile} \
        --output {params.output_dir} \
        >> {log} 2>&1
        '''


if config["params"]["profiling"]["metaphlan"]["do_v2"] and \
   config["params"]["profiling"]["humann"]["do_v2"]:
    rule profiling_humann2_all:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann2/humann2_{target}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann2/humann2_{target}_joined_{suffix}.tsv"),
               os.path.join(
                    config["output"]["profiling"],
                    "report/humann2/humann2_{target}_{norm}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann2/humann2_{target}_{norm}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann2/humann2.{group}.joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann2/humann2.{group}.joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm = config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   suffix=["stratified", "unstratified"])

else:
    rule profiling_humann2_all:
        input:


rule profiling_humann35_config:
    output:
        touch(os.path.join(config["output"]["profiling"], "config/humann35.config.done"))
    log:
        os.path.join(config["output"]["profiling"], "logs/humann35/humann35.config.log")
    conda:
        config["envs"]["biobakery3"]
    params:
        database_utility_mapping = config["params"]["profiling"]["humann"]["database_utility_mapping_v31"],
        database_nucleotide = config["params"]["profiling"]["humann"]["database_nucleotide_v31"],
        database_protein = config["params"]["profiling"]["humann"]["database_protein_v31"],
        threads = config["params"]["profiling"]["threads"]
    priority:
        20
    shell:
        '''
        humann_config > {log}
        humann_config --update database_folders utility_mapping {params.database_utility_mapping}
        humann_config --update database_folders nucleotide {params.database_nucleotide}
        humann_config --update database_folders protein {params.database_protein}
        humann_config --update run_modes threads {params.threads}
        echo "####" >> {log}
        humann_config >> {log}
        '''


localrules:
    profiling_humann35_config


rule profiling_humann35:
    input:
        tag = os.path.join(config["output"]["profiling"], "config/humann35.config.done"),
        reads = profiling_input_with_short_reads,
        profile = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv")
    output:
        genefamilies = os.path.join(
            config["output"]["profiling"],
            "profile/humann35/{sample}/{sample}_genefamilies.tsv"),
        pathabundance = os.path.join(
            config["output"]["profiling"],
            "profile/humann35/{sample}/{sample}_pathabundance.tsv"),
        pathcoverage = os.path.join(
            config["output"]["profiling"],
            "profile/humann35/{sample}/{sample}_pathcoverage.tsv")
    log:
        os.path.join(config["output"]["profiling"], "logs/humann35/{sample}.humann35.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/humann35/{sample}.humann35.benchmark.txt")
    priority:
        20
    conda:
        config["envs"]["biobakery3"]
    params:
        wrapper_dir = WRAPPER_DIR,
        basename = "{sample}",
        output_dir = os.path.join(config["output"]["profiling"], "profile/humann35/{sample}"),
        prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"],
        nucleotide_identity_threshold = config["params"]["profiling"]["humann"]["nucleotide_identity_threshold"],
        translated_identity_threshold = config["params"]["profiling"]["humann"]["translated_identity_threshold"],
        translated_subject_coverage_threshold = config["params"]["profiling"]["humann"]["translated_subject_coverage_threshold"],
        translated_query_coverage_threshold = config["params"]["profiling"]["humann"]["translated_query_coverage_threshold"],
        nucleotide_subject_coverage_threshold = config["params"]["profiling"]["humann"]["nucleotide_subject_coverage_threshold"],
        nucleotide_query_coverage_threshold = config["params"]["profiling"]["humann"]["nucleotide_query_coverage_threshold"],
        external_opts = config["params"]["profiling"]["humann"]["external_opts_v35"]
    threads:
        config["params"]["profiling"]["threads"]
    shell:
        '''
        rm -rf {params.output_dir}
        mkdir -p {params.output_dir}

        python {params.wrapper_dir}/misc.py \
        --basename {params.basename} \
        --input-file {input.reads} \
        --output-dir {params.output_dir}

        humann \
        --threads {threads} \
        --input {params.output_dir}/{params.basename}.fq.gz \
        --input-format fastq.gz \
        --taxonomic-profile {input.profile} \
        --prescreen-threshold {params.prescreen_threshold} \
        --nucleotide-identity-threshold {params.nucleotide_identity_threshold} \
        --translated-identity-threshold {params.translated_identity_threshold} \
        --translated-subject-coverage-threshold {params.translated_subject_coverage_threshold} \
        --nucleotide-subject-coverage-threshold {params.nucleotide_subject_coverage_threshold} \
        --translated-query-coverage-threshold {params.translated_query_coverage_threshold} \
        --nucleotide-query-coverage-threshold {params.nucleotide_query_coverage_threshold} \
        --output-basename {params.basename} \
        --output {params.output_dir} \
        --remove-temp-output \
        {params.external_opts} \
        --o-log {log}

        rm -rf {params.output_dir}/{params.basename}.fq.gz
        '''


rule profiling_humann35_postprocess:
    input:
        expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann35/{{sample}}/{{sample}}_{target}.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"])
    output:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann35/{{sample}}/{{sample}}_{target}_{norm}.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"],
            norm=config["params"]["profiling"]["humann"]["normalize_method"]),
        groupprofiles = expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann35/{{sample}}/{{sample}}_{group}_groupped.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann35/{sample}.humann35_postprocess.log")
    conda:
        config["envs"]["biobakery3"]
    params:
        wrapper_dir =WRAPPER_DIR,
        normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
        regroup_method = config["params"]["profiling"]["humann"]["regroup_method"],
        map_database =  config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        humann_renorm_table \
        --input {input[0]} \
        --update-snames \
        --output {output.targets[0]} \
        --units {params.normalize_method} \
        > {log} 2>&1

        humann_renorm_table \
        --input {input[1]} \
        --update-snames \
        --output {output.targets[1]} \
        --units {params.normalize_method} \
        >> {log} 2>&1

        humann_renorm_table \
        --input {input[2]} \
        --update-snames \
        --output {output.targets[2]} \
        --units {params.normalize_method} \
        >> {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        regroup_table \
        --input {input[0]} \
        --groups {params.map_database} \
        --function {params.regroup_method} \
        --output {output.groupprofiles} \
        >> {log} 2>&1
        '''


rule profiling_humann35_join:
    input:
        expand([
            os.path.join(
                config["output"]["profiling"],
                "profile/humann35/{sample}/{sample}_{target}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/humann35/{sample}/{sample}_{target}_{norm}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/humann35/{sample}/{sample}_{group}_groupped.tsv")],
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann"]["normalize_method"],
                group=config["params"]["profiling"]["humann"]["map_database"],
                sample=SAMPLES_ID_LIST)
    output:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann35/humann35_{target}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        targets_norm = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann35/humann35_{target}_{norm}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"],
            norm=config["params"]["profiling"]["humann"]["normalize_method"]),
        groupprofile = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann35/humann35.{group}.joined.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann35/humann35_join.log")
    conda:
        config["envs"]["biobakery3"]
    params:
        wrapper_dir =WRAPPER_DIR,
        input_dir = os.path.join(config["output"]["profiling"], "profile/humann35"),
        normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
        map_database = config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.targets} \
        --file_name genefamilies.tsv pathabundance.tsv pathcoverage.tsv \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.targets_norm} \
        --file_name \
        genefamilies_{params.normalize_method}.tsv \
        pathabundance_{params.normalize_method}.tsv \
        pathcoverage_{params.normalize_method}.tsv \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.groupprofile} \
        --file_name {params.map_database} \
        >> {log} 2>&1
        '''


rule profiling_humann35_split_stratified:
    input:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann35/humann35_{target}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        targets_norm = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann35/humann35_{target}_{norm}_joined.tsv"),
            norm = config["params"]["profiling"]["humann"]["normalize_method"],
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        groupprofile = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann35/humann35.{group}.joined.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    output:
        expand([
            os.path.join(config["output"]["profiling"], "report/humann35/humann35_{target}_joined_{suffix}.tsv"),
            os.path.join(config["output"]["profiling"], "report/humann35/humann35_{target}_{norm}_joined_{suffix}.tsv"),
            os.path.join(config["output"]["profiling"], "report/humann35/humann35.{group}.joined_{suffix}.tsv")],
            target=["genefamilies", "pathabundance", "pathcoverage"],
            norm = config["params"]["profiling"]["humann"]["normalize_method"],
            group=config["params"]["profiling"]["humann"]["map_database"],
            suffix=["stratified", "unstratified"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann35/humann35_split_stratified.log")
    conda:
        config["envs"]["biobakery3"]
    params:
        wrapper_dir = WRAPPER_DIR,
        output_dir = os.path.join(config["output"]["profiling"], "report/humann35"),
        map_database = config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.targets} \
        --output {params.output_dir} \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.targets_norm} \
        --output {params.output_dir} \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.groupprofile} \
        --output {params.output_dir} \
        >> {log} 2>&1
        '''


if config["params"]["profiling"]["metaphlan"]["do_v3"] and \
   config["params"]["profiling"]["humann"]["do_v35"]:
    rule profiling_humann35_all:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann35/humann35_{target}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann35/humann35_{target}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann35/humann35_{target}_{norm}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann35/humann35_{target}_{norm}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann35/humann35.{group}.joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann35/humann35.{group}.joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm = config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   suffix=["stratified", "unstratified"])

else:
    rule profiling_humann35_all:
        input:


rule profiling_humann38_config:
    output:
        touch(os.path.join(config["output"]["profiling"], "config/humann38.config.done"))
    log:
        os.path.join(config["output"]["profiling"], "logs/humann38/humann38.config.log")
    conda:
        config["envs"]["biobakery40"]
    params:
        database_utility_mapping = config["params"]["profiling"]["humann"]["database_utility_mapping_v31"],
        database_nucleotide = config["params"]["profiling"]["humann"]["database_nucleotide_v31"],
        database_protein = config["params"]["profiling"]["humann"]["database_protein_v31"],
        threads = config["params"]["profiling"]["threads"]
    priority:
        20
    shell:
        '''
        humann_config > {log}
        humann_config --update database_folders utility_mapping {params.database_utility_mapping}
        humann_config --update database_folders nucleotide {params.database_nucleotide}
        humann_config --update database_folders protein {params.database_protein}
        humann_config --update run_modes threads {params.threads}
        echo "####" >> {log}
        humann_config >> {log}
        '''


localrules:
    profiling_humann38_config


rule profiling_humann38:
    input:
        tag = os.path.join(config["output"]["profiling"], "config/humann38.config.done"),
        reads = profiling_input_with_short_reads,
        profile = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan40/{sample}/{sample}.metaphlan40.abundance.profile.tsv")
    output:
        genefamilies = os.path.join(
            config["output"]["profiling"],
            "profile/humann38/{sample}/{sample}_genefamilies.tsv"),
        pathabundance = os.path.join(
            config["output"]["profiling"],
            "profile/humann38/{sample}/{sample}_pathabundance.tsv"),
        pathcoverage = os.path.join(
            config["output"]["profiling"],
            "profile/humann38/{sample}/{sample}_pathcoverage.tsv")
    log:
        os.path.join(config["output"]["profiling"], "logs/humann38/{sample}.humann38.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/humann38/{sample}.humann38.benchmark.txt")
    priority:
        20
    conda:
        config["envs"]["biobakery40"]
    params:
        wrapper_dir = WRAPPER_DIR,
        basename = "{sample}",
        output_dir = os.path.join(config["output"]["profiling"], "profile/humann38/{sample}"),
        prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"],
        nucleotide_identity_threshold = config["params"]["profiling"]["humann"]["nucleotide_identity_threshold"],
        translated_identity_threshold = config["params"]["profiling"]["humann"]["translated_identity_threshold"],
        translated_subject_coverage_threshold = config["params"]["profiling"]["humann"]["translated_subject_coverage_threshold"],
        translated_query_coverage_threshold = config["params"]["profiling"]["humann"]["translated_query_coverage_threshold"],
        nucleotide_subject_coverage_threshold = config["params"]["profiling"]["humann"]["nucleotide_subject_coverage_threshold"],
        nucleotide_query_coverage_threshold = config["params"]["profiling"]["humann"]["nucleotide_query_coverage_threshold"],
        external_opts = config["params"]["profiling"]["humann"]["external_opts_v38"]
    threads:
        config["params"]["profiling"]["threads"]
    shell:
        '''
        rm -rf {params.output_dir}
        mkdir -p {params.output_dir}

        python {params.wrapper_dir}/misc.py \
        --basename {params.basename} \
        --input-file {input.reads} \
        --output-dir {params.output_dir}

        humann \
        --threads {threads} \
        --input {params.output_dir}/{params.basename}.fq.gz \
        --input-format fastq.gz \
        --taxonomic-profile {input.profile} \
        --prescreen-threshold {params.prescreen_threshold} \
        --nucleotide-identity-threshold {params.nucleotide_identity_threshold} \
        --translated-identity-threshold {params.translated_identity_threshold} \
        --translated-subject-coverage-threshold {params.translated_subject_coverage_threshold} \
        --nucleotide-subject-coverage-threshold {params.nucleotide_subject_coverage_threshold} \
        --translated-query-coverage-threshold {params.translated_query_coverage_threshold} \
        --nucleotide-query-coverage-threshold {params.nucleotide_query_coverage_threshold} \
        --output-basename {params.basename} \
        --output {params.output_dir} \
        --remove-temp-output \
        {params.external_opts} \
        --o-log {log}

        rm -rf {params.output_dir}/{params.basename}.fq.gz
        '''


rule profiling_humann38_postprocess:
    input:
        expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann38/{{sample}}/{{sample}}_{target}.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"])
    output:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann38/{{sample}}/{{sample}}_{target}_{norm}.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"],
            norm=config["params"]["profiling"]["humann"]["normalize_method"]),
        groupprofiles = expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann38/{{sample}}/{{sample}}_{group}_groupped.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann38/{sample}.humann38_postprocess.log")
    conda:
        config["envs"]["biobakery40"]
    params:
        wrapper_dir =WRAPPER_DIR,
        normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
        regroup_method = config["params"]["profiling"]["humann"]["regroup_method"],
        map_database =  config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        humann_renorm_table \
        --input {input[0]} \
        --update-snames \
        --output {output.targets[0]} \
        --units {params.normalize_method} \
        > {log} 2>&1

        humann_renorm_table \
        --input {input[1]} \
        --update-snames \
        --output {output.targets[1]} \
        --units {params.normalize_method} \
        >> {log} 2>&1

        humann_renorm_table \
        --input {input[2]} \
        --update-snames \
        --output {output.targets[2]} \
        --units {params.normalize_method} \
        >> {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        regroup_table \
        --input {input[0]} \
        --groups {params.map_database} \
        --function {params.regroup_method} \
        --output {output.groupprofiles} \
        >> {log} 2>&1
        '''


rule profiling_humann38_join:
    input:
        expand([
            os.path.join(
                config["output"]["profiling"],
                "profile/humann38/{sample}/{sample}_{target}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/humann38/{sample}/{sample}_{target}_{norm}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/humann38/{sample}/{sample}_{group}_groupped.tsv")],
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann"]["normalize_method"],
                group=config["params"]["profiling"]["humann"]["map_database"],
                sample=SAMPLES_ID_LIST)
    output:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann38/humann38_{target}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        targets_norm = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann38/humann38_{target}_{norm}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"],
            norm=config["params"]["profiling"]["humann"]["normalize_method"]),
        groupprofile = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann38/humann38.{group}.joined.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann38/humann38_join.log")
    conda:
        config["envs"]["biobakery40"]
    params:
        wrapper_dir =WRAPPER_DIR,
        input_dir = os.path.join(config["output"]["profiling"], "profile/humann38"),
        normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
        map_database = config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.targets} \
        --file_name genefamilies.tsv pathabundance.tsv pathcoverage.tsv \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.targets_norm} \
        --file_name \
        genefamilies_{params.normalize_method}.tsv \
        pathabundance_{params.normalize_method}.tsv \
        pathcoverage_{params.normalize_method}.tsv \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.groupprofile} \
        --file_name {params.map_database} \
        >> {log} 2>&1
        '''


rule profiling_humann38_split_stratified:
    input:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann38/humann38_{target}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        targets_norm = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann38/humann38_{target}_{norm}_joined.tsv"),
            norm = config["params"]["profiling"]["humann"]["normalize_method"],
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        groupprofile = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann38/humann38.{group}.joined.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    output:
        expand([
            os.path.join(
                config["output"]["profiling"],
                "report/humann38/humann38_{target}_joined_{suffix}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "report/humann38/humann38_{target}_{norm}_joined_{suffix}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "report/humann38/humann38.{group}.joined_{suffix}.tsv")],
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm = config["params"]["profiling"]["humann"]["normalize_method"],
                group=config["params"]["profiling"]["humann"]["map_database"],
                suffix=["stratified", "unstratified"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann38/humann38_split_stratified.log")
    conda:
        config["envs"]["biobakery40"]
    params:
        wrapper_dir = WRAPPER_DIR,
        output_dir = os.path.join(config["output"]["profiling"], "report/humann38"),
        map_database = config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.targets} \
        --output {params.output_dir} \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.targets_norm} \
        --output {params.output_dir} \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.groupprofile} \
        --output {params.output_dir} \
        >> {log} 2>&1
        '''


if config["params"]["profiling"]["metaphlan"]["do_v40"] and \
   config["params"]["profiling"]["humann"]["do_v38"]:
        rule profiling_humann38_all:
            input:
                expand([
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann38/humann38_{target}_joined.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann38/humann38_{target}_joined_{suffix}.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann38/humann38_{target}_{norm}_joined.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann38/humann38_{target}_{norm}_joined_{suffix}.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann38/humann38.{group}.joined.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann38/humann38.{group}.joined_{suffix}.tsv")],
                    target=["genefamilies", "pathabundance", "pathcoverage"],
                    norm = config["params"]["profiling"]["humann"]["normalize_method"],
                    group=config["params"]["profiling"]["humann"]["map_database"],
                    suffix=["stratified", "unstratified"])

else:
    rule profiling_humann38_all:
        input:


rule profiling_humann39_config:
    output:
        touch(os.path.join(config["output"]["profiling"], "config/humann39.config.done"))
    log:
        os.path.join(config["output"]["profiling"], "logs/humann39/humann39.config.log")
    conda:
        config["envs"]["biobakery41"]
    params:
        database_utility_mapping = config["params"]["profiling"]["humann"]["database_utility_mapping_v31"],
        database_nucleotide = config["params"]["profiling"]["humann"]["database_nucleotide_v31"],
        database_protein = config["params"]["profiling"]["humann"]["database_protein_v31"],
        threads = config["params"]["profiling"]["threads"]
    priority:
        20
    shell:
        '''
        humann_config > {log}
        humann_config --update database_folders utility_mapping {params.database_utility_mapping}
        humann_config --update database_folders nucleotide {params.database_nucleotide}
        humann_config --update database_folders protein {params.database_protein}
        humann_config --update run_modes threads {params.threads}
        echo "####" >> {log}
        humann_config >> {log}
        '''


localrules:
    profiling_humann39_config


rule profiling_humann39:
    input:
        tag = os.path.join(config["output"]["profiling"], "config/humann39.config.done"),
        reads = profiling_input_with_short_reads,
        profile = os.path.join(
            config["output"]["profiling"],
            "profile/metaphlan41/{sample}/{sample}.metaphlan41.abundance.profile.tsv")
    output:
        genefamilies = os.path.join(
            config["output"]["profiling"],
            "profile/humann39/{sample}/{sample}_genefamilies.tsv"),
        pathabundance = os.path.join(
            config["output"]["profiling"],
            "profile/humann39/{sample}/{sample}_pathabundance.tsv"),
        pathcoverage = os.path.join(
            config["output"]["profiling"],
            "profile/humann39/{sample}/{sample}_pathcoverage.tsv")
    log:
        os.path.join(config["output"]["profiling"], "logs/humann39/{sample}.humann39.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/humann39/{sample}.humann39.benchmark.txt")
    priority:
        20
    conda:
        config["envs"]["biobakery41"]
    params:
        wrapper_dir = WRAPPER_DIR,
        basename = "{sample}",
        output_dir = os.path.join(config["output"]["profiling"], "profile/humann39/{sample}"),
        prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"],
        nucleotide_identity_threshold = config["params"]["profiling"]["humann"]["nucleotide_identity_threshold"],
        translated_identity_threshold = config["params"]["profiling"]["humann"]["translated_identity_threshold"],
        translated_subject_coverage_threshold = config["params"]["profiling"]["humann"]["translated_subject_coverage_threshold"],
        translated_query_coverage_threshold = config["params"]["profiling"]["humann"]["translated_query_coverage_threshold"],
        nucleotide_subject_coverage_threshold = config["params"]["profiling"]["humann"]["nucleotide_subject_coverage_threshold"],
        nucleotide_query_coverage_threshold = config["params"]["profiling"]["humann"]["nucleotide_query_coverage_threshold"],
        external_opts = config["params"]["profiling"]["humann"]["external_opts_v39"]
    threads:
        config["params"]["profiling"]["threads"]
    shell:
        '''
        rm -rf {params.output_dir}
        mkdir -p {params.output_dir}

        python {params.wrapper_dir}/misc.py \
        --basename {params.basename} \
        --input-file {input.reads} \
        --output-dir {params.output_dir}

        humann \
        --threads {threads} \
        --input {params.output_dir}/{params.basename}.fq.gz \
        --input-format fastq.gz \
        --taxonomic-profile {input.profile} \
        --prescreen-threshold {params.prescreen_threshold} \
        --nucleotide-identity-threshold {params.nucleotide_identity_threshold} \
        --translated-identity-threshold {params.translated_identity_threshold} \
        --translated-subject-coverage-threshold {params.translated_subject_coverage_threshold} \
        --nucleotide-subject-coverage-threshold {params.nucleotide_subject_coverage_threshold} \
        --translated-query-coverage-threshold {params.translated_query_coverage_threshold} \
        --nucleotide-query-coverage-threshold {params.nucleotide_query_coverage_threshold} \
        --output-basename {params.basename} \
        --output {params.output_dir} \
        --remove-temp-output \
        {params.external_opts} \
        --o-log {log}

        rm -rf {params.output_dir}/{params.basename}.fq.gz
        '''


rule profiling_humann39_postprocess:
    input:
        expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann39/{{sample}}/{{sample}}_{target}.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"])
    output:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann39/{{sample}}/{{sample}}_{target}_{norm}.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"],
            norm=config["params"]["profiling"]["humann"]["normalize_method"]),
        groupprofiles = expand(os.path.join(
            config["output"]["profiling"],
            "profile/humann39/{{sample}}/{{sample}}_{group}_groupped.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann39/{sample}.humann39_postprocess.log")
    conda:
        config["envs"]["biobakery41"]
    params:
        wrapper_dir =WRAPPER_DIR,
        normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
        regroup_method = config["params"]["profiling"]["humann"]["regroup_method"],
        map_database =  config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        humann_renorm_table \
        --input {input[0]} \
        --update-snames \
        --output {output.targets[0]} \
        --units {params.normalize_method} \
        > {log} 2>&1

        humann_renorm_table \
        --input {input[1]} \
        --update-snames \
        --output {output.targets[1]} \
        --units {params.normalize_method} \
        >> {log} 2>&1

        humann_renorm_table \
        --input {input[2]} \
        --update-snames \
        --output {output.targets[2]} \
        --units {params.normalize_method} \
        >> {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        regroup_table \
        --input {input[0]} \
        --groups {params.map_database} \
        --function {params.regroup_method} \
        --output {output.groupprofiles} \
        >> {log} 2>&1
        '''


rule profiling_humann39_join:
    input:
        expand([
            os.path.join(
                config["output"]["profiling"],
                "profile/humann39/{sample}/{sample}_{target}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/humann39/{sample}/{sample}_{target}_{norm}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/humann39/{sample}/{sample}_{group}_groupped.tsv")],
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann"]["normalize_method"],
                group=config["params"]["profiling"]["humann"]["map_database"],
                sample=SAMPLES_ID_LIST)
    output:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann39/humann39_{target}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        targets_norm = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann39/humann39_{target}_{norm}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"],
            norm=config["params"]["profiling"]["humann"]["normalize_method"]),
        groupprofile = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann39/humann39.{group}.joined.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann39/humann39_join.log")
    conda:
        config["envs"]["biobakery41"]
    params:
        wrapper_dir =WRAPPER_DIR,
        input_dir = os.path.join(config["output"]["profiling"], "profile/humann39"),
        normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
        map_database = config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.targets} \
        --file_name genefamilies.tsv pathabundance.tsv pathcoverage.tsv \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.targets_norm} \
        --file_name \
        genefamilies_{params.normalize_method}.tsv \
        pathabundance_{params.normalize_method}.tsv \
        pathcoverage_{params.normalize_method}.tsv \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        join_tables \
        --input {params.input_dir} \
        --output {output.groupprofile} \
        --file_name {params.map_database} \
        >> {log} 2>&1
        '''


rule profiling_humann39_split_stratified:
    input:
        targets = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann39/humann39_{target}_joined.tsv"),
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        targets_norm = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann39/humann39_{target}_{norm}_joined.tsv"),
            norm = config["params"]["profiling"]["humann"]["normalize_method"],
            target=["genefamilies", "pathabundance", "pathcoverage"]),
        groupprofile = expand(os.path.join(
            config["output"]["profiling"],
            "report/humann39/humann39.{group}.joined.tsv"),
            group=config["params"]["profiling"]["humann"]["map_database"])
    output:
        expand([
            os.path.join(
                config["output"]["profiling"],
                "report/humann39/humann39_{target}_joined_{suffix}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "report/humann39/humann39_{target}_{norm}_joined_{suffix}.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "report/humann39/humann39.{group}.joined_{suffix}.tsv")],
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm = config["params"]["profiling"]["humann"]["normalize_method"],
                group=config["params"]["profiling"]["humann"]["map_database"],
                suffix=["stratified", "unstratified"])
    log:
        os.path.join(config["output"]["profiling"], "logs/humann39/humann39_split_stratified.log")
    conda:
        config["envs"]["biobakery41"]
    params:
        wrapper_dir = WRAPPER_DIR,
        output_dir = os.path.join(config["output"]["profiling"], "report/humann39"),
        map_database = config["params"]["profiling"]["humann"]["map_database"]
    priority:
        20
    shell:
        '''
        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.targets} \
        --output {params.output_dir} \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.targets_norm} \
        --output {params.output_dir} \
        > {log} 2>&1

        python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
        split_stratified_table \
        --input {input.groupprofile} \
        --output {params.output_dir} \
        >> {log} 2>&1
        '''


if config["params"]["profiling"]["metaphlan"]["do_v41"] and \
   config["params"]["profiling"]["humann"]["do_v39"]:
        rule profiling_humann39_all:
            input:
                expand([
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann39/humann39_{target}_joined.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann39/humann39_{target}_joined_{suffix}.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann39/humann39_{target}_{norm}_joined.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann39/humann39_{target}_{norm}_joined_{suffix}.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann39/humann39.{group}.joined.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/humann39/humann39.{group}.joined_{suffix}.tsv")],
                    target=["genefamilies", "pathabundance", "pathcoverage"],
                    norm = config["params"]["profiling"]["humann"]["normalize_method"],
                    group=config["params"]["profiling"]["humann"]["map_database"],
                    suffix=["stratified", "unstratified"])

else:
    rule profiling_humann39_all:
        input:
