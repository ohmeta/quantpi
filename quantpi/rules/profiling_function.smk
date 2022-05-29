if config["params"]["profiling"]["metaphlan"]["do_v2"] and \
   config["params"]["profiling"]["humann"]["do_v2"]:
    if config["params"]["profiling"]["humann"]["update_config"]:
        rule profiling_humann2_config:
            output:
                touch(os.path.join(config["output"]["profiling"], ".humann2.config.done"))
            log:
                os.path.join(config["output"]["profiling"], "logs/humann2/humann2.config.log")
            conda:
                config["envs"]["humann2"]
            params:
                database_utility_mapping = config["params"]["profiling"]["humann"]["database_utility_mapping"],
                database_nucleotide = config["params"]["profiling"]["humann"]["database_nucleotide"],
                database_protein = config["params"]["profiling"]["humann"]["database_protein"],
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
    else:
        rule profiling_humann2_config:
            output:
                touch(os.path.join(config["output"]["profiling"], ".humann2.config.done"))
            shell:
                '''
                echo "hello"
                '''

    localrules:
        profiling_humann2_config


    rule profiling_humann2_build_chocophlan_pangenome_db:
        input:
            tag = os.path.join(config["output"]["profiling"], ".humann2.config.done"),
            profile = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv")
        output:
            expand(os.path.join(
                config["output"]["profiling"],
                "database/humann2/{{sample}}/{{sample}}_bowtie2_index.{suffix}"),
                   suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        conda:
            config["envs"]["humann2"]
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/humann2/{sample}.humann2.build_pandb.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/humann2/{sample}.bowtie2_index.benchmark.txt")
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
            tag = os.path.join(config["output"]["profiling"], ".humann2.config.done"),
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
            os.path.join(config["output"]["profiling"],
                         "benchmark/humann2/{sample}.humann2.benchmark.txt")
        conda:
            config["envs"]["humann2"]
        params:
            basename = "{sample}",
            index = os.path.join(config["output"]["profiling"],
                                 "database/humann2/{sample}/{sample}_bowtie2_index"),
            evalue = config["params"]["profiling"]["humann"]["evalue"],
            prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"],
            identity_threshold = config["params"]["profiling"]["humann"]["identity_threshold"],
            translated_subject_coverage_threshold = \
                config["params"]["profiling"]["humann"]["translated_subject_coverage_threshold"],
            translated_query_coverage_threshold = \
                config["params"]["profiling"]["humann"]["translated_query_coverage_threshold"],
            xipe = "on" if config["params"]["profiling"]["humann"]["xipe"] else "off",
            minpath = "on" if config["params"]["profiling"]["humann"]["minpath"] else "off",
            pick_frames = "on" if config["params"]["profiling"]["humann"]["pick_frames"] else "off",
            gap_fill = "on" if config["params"]["profiling"]["humann"]["gap_fill"] else "off",
            remove_temp_output = "--remove-temp-output" \
                if config["params"]["profiling"]["humann"]["remove_temp_output"] \
                   else "",
            memory_use = config["params"]["profiling"]["humann"]["memory_use"],
            output_dir = os.path.join(config["output"]["profiling"],
                                      "profile/humann2/{sample}")
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
            --evalue {params.evalue} \
            --prescreen-threshold {params.prescreen_threshold} \
            --identity-threshold {params.identity_threshold} \
            --translated-subject-coverage-threshold {params.translated_subject_coverage_threshold} \
            --translated-query-coverage-threshold {params.translated_query_coverage_threshold} \
            --xipe {params.xipe} \
            --minpath {params.minpath} \
            --pick-frames {params.pick_frames} \
            --gap-fill {params.gap_fill} \
            --memory-use {params.memory_use} \
            --output-basename {params.basename} \
            --output {params.output_dir} \
            {params.remove_temp_output} \
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
            os.path.join(config["output"]["profiling"],
                         "logs/humann2/{sample}.humann2_postprocess.log")
        conda:
            config["envs"]["humann2"]
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
            os.path.join(config["output"]["profiling"],
                         "logs/humann2/humann2_join.log")
        conda:
            config["envs"]["humann2"]
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
            config["envs"]["humann2"]
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
                   suffix=["stratified", "unstratified"]),

            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_humann2_all:
        input:


if config["params"]["profiling"]["metaphlan"]["do_v3"] and \
   config["params"]["profiling"]["humann"]["do_v3"]:
    if config["params"]["profiling"]["humann"]["update_config"]:
        rule profiling_humann3_config:
            output:
                touch(os.path.join(config["output"]["profiling"], ".humann3.config.done"))
            log:
                os.path.join(config["output"]["profiling"], "logs/humann3/humann3.config.log")
            conda:
                config["envs"]["humann3"]
            params:
                database_utility_mapping = config["params"]["profiling"]["humann"]["database_utility_mapping_v3"],
                database_nucleotide = config["params"]["profiling"]["humann"]["database_nucleotide_v3"],
                database_protein = config["params"]["profiling"]["humann"]["database_protein_v3"],
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
    else:
        rule profiling_humann3_config:
            output:
                touch(os.path.join(config["output"]["profiling"], ".humann3.config.done"))
            shell:
                '''
                echo "hello"
                '''

    localrules:
        profiling_humann3_config


    #rule profiling_humann3_build_chocophlan_pangenome_db:
    #    input:
    #        profile = os.path.join(
    #            config["output"]["profiling"],
    #            "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv")
    #    output:
    #        expand(os.path.join(
    #            config["output"]["profiling"],
    #            "database/humann3/{{sample}}/{{sample}}_bowtie2_index.{suffix}"),
    #               suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    #    log:
    #        os.path.join(config["output"]["profiling"],
    #                     "logs/humann3/{sample}.humann3.build_pandb.log")
    #    params:
    #        basename = "{sample}",
    #        wrapper_dir = WRAPPER_DIR,
    #        db_dir = os.path.join(config["output"]["profiling"], "database/humann3/{sample}"),
    #        prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"]
    #    shell:
    #        '''
    #        python {params.wrapper_dir}/humann3_db_wrapper.py \
    #        --log {log} \
    #        --basename {params.basename} \
    #        --db_dir {params.db_dir} \
    #        --prescreen_threshold {params.prescreen_threshold} \
    #        --taxonomic_profile {input.profile}
    #        '''


    rule profiling_humann3:
        input:
            tag = os.path.join(config["output"]["profiling"], ".humann3.config.done"),
            reads = profiling_input_with_short_reads,
            profile = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv")
        output:
            genefamilies = os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{sample}/{sample}_genefamilies.tsv"),
            pathabundance = os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{sample}/{sample}_pathabundance.tsv"),
            pathcoverage = os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{sample}/{sample}_pathcoverage.tsv")
        log:
            os.path.join(config["output"]["profiling"], "logs/humann3/{sample}.humann3.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/humann3/{sample}.humann3.benchmark.txt")
        priority:
            20
        conda:
            config["envs"]["humann3"]
        params:
            basename = "{sample}",
            wrapper_dir = WRAPPER_DIR,
            index = os.path.join(config["output"]["profiling"],
                                 "database/humann3/{sample}/{sample}_bowtie2_index"),
            evalue = config["params"]["profiling"]["humann"]["evalue"],
            prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"],
            nucleotide_identity_threshold = \
                config["params"]["profiling"]["humann"]["nucleotide_identity_threshold"],
            translated_identity_threshold = \
                config["params"]["profiling"]["humann"]["translated_identity_threshold"],
            translated_subject_coverage_threshold = \
                config["params"]["profiling"]["humann"]["translated_subject_coverage_threshold"],
            translated_query_coverage_threshold = \
                config["params"]["profiling"]["humann"]["translated_query_coverage_threshold"],
            nucleotide_subject_coverage_threshold = \
                config["params"]["profiling"]["humann"]["nucleotide_subject_coverage_threshold"],
            nucleotide_query_coverage_threshold = \
                config["params"]["profiling"]["humann"]["nucleotide_query_coverage_threshold"],
            minpath = "on" if config["params"]["profiling"]["humann"]["minpath"] else "off",
            gap_fill = "on" if config["params"]["profiling"]["humann"]["gap_fill"] else "off",
            xipe = "on" if config["params"]["profiling"]["humann"]["xipe"] else "off",
            pathways = config["params"]["profiling"]["humann"]["pathways"],
            remove_temp_output = "--remove-temp-output" \
                if config["params"]["profiling"]["humann"]["remove_temp_output"] \
                   else "",
            memory_use = config["params"]["profiling"]["humann"]["memory_use"],
            output_dir = os.path.join(config["output"]["profiling"],
                                      "profile/humann3/{sample}")
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            mkdir -p {params.output_dir}

            python {params.wrapper_dir}/misc.py \
            --basename {params.basename} \
            --input-file {input.reads} \
            --output-dir {params.output_dir}

            rm -rf {params.output_dir}/{params.basename}_humann_temp_*

            humann \
            --threads {threads} \
            --input {params.output_dir}/{params.basename}.fq.gz \
            --input-format fastq.gz \
            --taxonomic-profile {input.profile} \
            --evalue {params.evalue} \
            --prescreen-threshold {params.prescreen_threshold} \
            --nucleotide-identity-threshold {params.nucleotide_identity_threshold} \
            --translated-identity-threshold {params.translated_identity_threshold} \
            --translated-subject-coverage-threshold {params.translated_subject_coverage_threshold} \
            --nucleotide-subject-coverage-threshold {params.nucleotide_subject_coverage_threshold} \
            --translated-query-coverage-threshold {params.translated_query_coverage_threshold} \
            --nucleotide-query-coverage-threshold {params.nucleotide_query_coverage_threshold} \
            --gap-fill {params.gap_fill} \
            --minpath {params.minpath} \
            --xipe {params.xipe} \
            --pathways {params.pathways} \
            --memory-use {params.memory_use} \
            --output-basename {params.basename} \
            --output {params.output_dir} \
            {params.remove_temp_output} \
            --o-log {log}

            rm -rf {params.output_dir}/{params.basename}.fq.gz
            '''


    rule profiling_humann3_postprocess:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{{sample}}/{{sample}}_{target}.tsv"),
                   target=["genefamilies", "pathabundance", "pathcoverage"])
        output:
            targets = expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{{sample}}/{{sample}}_{target}_{norm}.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann"]["normalize_method"]),
            groupprofiles = expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{{sample}}/{{sample}}_{group}_groupped.tsv"),
                group=config["params"]["profiling"]["humann"]["map_database"])
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/humann3/{sample}.humann3_postprocess.log")
        conda:
            config["envs"]["humann3"]
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


    rule profiling_humann3_join:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3/{sample}/{sample}_{target}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3/{sample}/{sample}_{target}_{norm}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3/{sample}/{sample}_{group}_groupped.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm=config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   sample=SAMPLES_ID_LIST)
        output:
            targets = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3_{target}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            targets_norm = expand(
               os.path.join(
                   config["output"]["profiling"],
                   "report/humann3/humann3_{target}_{norm}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann"]["normalize_method"]),
            groupprofile = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3.{group}.joined.tsv"),
                group=config["params"]["profiling"]["humann"]["map_database"])
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/humann3/humann3_join.log")
        conda:
            config["envs"]["humann3"]
        params:
            wrapper_dir =WRAPPER_DIR,
            input_dir = os.path.join(config["output"]["profiling"], "profile/humann3"),
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


    rule profiling_humann3_split_stratified:
        input:
            targets = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3_{target}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            targets_norm = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3_{target}_{norm}_joined.tsv"),
                norm = config["params"]["profiling"]["humann"]["normalize_method"],
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            groupprofile = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3.{group}.joined.tsv"),
                group=config["params"]["profiling"]["humann"]["map_database"])
        output:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3_{target}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3_{target}_{norm}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3.{group}.joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm = config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   suffix=["stratified", "unstratified"])
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/humann3/humann3_split_stratified.log")
        conda:
            config["envs"]["humann3"]
        params:
            wrapper_dir = WRAPPER_DIR,
            output_dir = os.path.join(config["output"]["profiling"], "report/humann3"),
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


    rule profiling_humann3_all:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3_{target}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3_{target}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3_{target}_{norm}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3_{target}_{norm}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3.{group}.joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/humann3/humann3.{group}.joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm = config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   suffix=["stratified", "unstratified"]),

            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_humann3_all:
        input:


rule profiling_all:
    input:
        rules.profiling_kraken2_all.input,
        rules.profiling_bracken_all.input,
        rules.profiling_metaphlan2_all.input,
        rules.profiling_metaphlan3_all.input,
        rules.profiling_kmcp_all.input,
        rules.profiling_bgi_soap_all.input,
        rules.profiling_bowtie2_all.input,
        rules.profiling_jgi_all.input,
        rules.profiling_humann2_all.input,
        rules.profiling_humann3_all.input


localrules:
    profiling_kraken2_all,
    profiling_bracken_all,
    profiling_metaphlan2_all,
    profiling_metaphlan3_all,
    profiling_kmcp_all,
    profiling_bgi_soap_all,
    profiling_bowtie2_all,
    profiling_jgi_all,
    profiling_humann2_all,
    profiling_humann3_all,
    profiling_all