def profiling_input_with_short_reads(wildcards, have_single=False):
    if RMHOST_DO:
        return get_reads(wildcards, "rmhost", have_single)
    if TRIMMING_DO:
        return get_reads(wildcards, "trimming", have_single)
    else:
        return get_reads(wildcards, "raw", have_single)


if config["params"]["profiling"]["kraken2"]["do"]:
    rule profiling_kraken2:
        input:
            reads = profiling_input_with_short_reads
        output:
            report = os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report"),
            report_mpa_reads_count = os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report.mpa.reads_count"),
            report_mpa_percentages = os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report.mpa.percentages")
        log:
            os.path.join(
                config["output"]["profiling"], "logs/kraken2/{sample}.kraken2.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"], "benchmark/kraken2/{sample}.kraken2.benchmark.txt")
        params:
            save_table = config["params"]["profiling"]["kraken2"]["save_table"],
            paired = "--paired" if IS_PE else "",
            database = config["params"]["profiling"]["kraken2"]["database"],
            quick = "--quick" \
                if config["params"]["profiling"]["kraken2"]["quick"] \
                else "",
            memory_mapping = "--memory-mapping" \
                if config["params"]["profiling"]["kraken2"]["memory_mapping"] \
                else "",
            use_names = "--use-names" \
                if config["params"]["profiling"]["kraken2"]["use_names"] \
                else "",
            use_mpa_style = "--use-mpa-style" \
                if config["params"]["profiling"]["kraken2"]["use_mpa_style"] \
                else "",
            report_zero_counts = "--report-zero-counts" \
                if config["params"]["profiling"]["kraken2"]["report_zero_counts"] \
                else "",
            confidence = config["params"]["profiling"]["kraken2"]["confidence"],
            min_base_quality = config["params"]["profiling"]["kraken2"]["min_base_quality"],
            min_hit_groups = config["params"]["profiling"]["kraken2"]["min_hit_groups"],
            unclassified_out = "--unclassified-out %s" % \
                os.path.join(
                    config["output"]["profiling"],
                    "profile/kraken2/{sample}/{sample}.kraken2.unclassified%s.fq" \
                    % "#" if IS_PE else "") \
                    if config["params"]["profiling"]["kraken2"]["unclassified_out"] \
                    else "",
            classified_out = "--classified-out %s" % \
                os.path.join(
                    config["output"]["profiling"],
                    "profile/kraken2/{sample}/{sample}.kraken2.classified%s.fq" \
                    % "#" if IS_PE else "") \
                    if config["params"]["profiling"]["kraken2"]["classified_out"] \
                    else "",
            table = "--output %s" % \
                os.path.join(
                    config["output"]["profiling"],
                    "profile/kraken2/{sample}/{sample}.kraken2.table") \
                if config["params"]["profiling"]["kraken2"]["save_table"] \
                    else "--output /dev/null"
        threads:
            config["params"]["profiling"]["threads"]
        resources:
            mem_mb=config["params"]["profiling"]["kraken2"]["mem_mb"]
        conda:
            config["envs"]["kraken2"]
        shell:
            '''
            kraken2 \
            {params.quick} \
            {params.memory_mapping} \
            {params.use_mpa_style} \
            {params.use_names} \
            {params.report_zero_counts} \
            --threads {threads} \
            --db {params.database} \
            --confidence {params.confidence} \
            --minimum-base-quality {params.min_base_quality} \
            --minimum-hit-groups {params.min_hit_groups} \
            {params.unclassified_out} \
            {params.classified_out} \
            {params.table} \
            --report {output.report} \
            --gzip-compressed \
            {params.paired} \
            {input.reads} \
            2> {log}

            kreport2mpa.py \
            --report-file {output.report} \
            --display-header \
            --no-intermediate-ranks \
            --read_count \
            --output {output.report_mpa_reads_count}

            kreport2mpa.py \
            --report {output.report} \
            --no-intermediate-ranks \
            --percentages \
            --output {output.report_mpa_percentages}

            if [ "{params.save_table}" == "True" ];
            then
                tablef=`awk '{{print $2}}' {params.table}`
                pigz $tablef
            fi
            '''


    rule profiling_kraken2_combine_kreport:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/kraken2/{sample}/{sample}.kraken2.report"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(
                config["output"]["profiling"],
                "report/kraken2/kraken2_report.all.tsv")
        params:
            samples_name = " ".join(list(SAMPLES_ID_LIST))
        conda:
            config["envs"]["kraken2"]
        log:
            os.path.join(config["output"]["profiling"], "logs/krakentools/combine_kreports.log")
        threads:
            1
        resources:
            mem_mb=config["params"]["profiling"]["kraken2"]["mem_mb"]
        shell:
            '''
            combine_kreports.py \
            --report-file {input} \
            --sample-names {params.samples_name} \
            --display-headers \
            --output {output} \
            > {log} 2>&1
            '''


    rule profiling_kraken2_combine_kreport_mpa:
        input:
            report_mpa_reads_count = expand(os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report.mpa.reads_count"),
                sample=SAMPLES_ID_LIST),
            report_mpa_percentages = expand(os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report.mpa.percentages"),
                sample=SAMPLES_ID_LIST)
        output:
            report_mpa_reads_count = os.path.join(
                config["output"]["profiling"],
                "report/kraken2/kraken2_report.mpa.reads_count.tsv"),
            report_mpa_percentages = os.path.join(
                config["output"]["profiling"],
                "report/kraken2/kraken2_report.mpa.percentages.tsv")
        conda:
            config["envs"]["kraken2"]
        log:
            os.path.join(
                config["output"]["profiling"], "logs/krakentools/combine_kreports_mpa.log")
        threads:
            1
        resources:
            mem_mb=config["params"]["profiling"]["kraken2"]["mem_mb"]
        shell:
            '''
            echo "process 1:" > {log} 2>&1

            combine_mpa.py \
            --input {input.report_mpa_reads_count} \
            --output {output.report_mpa_reads_count} \
            >> {log} 2>&1

            echo "process 1:" >> {log} 2>&1

            combine_mpa.py \
            --input {input.report_mpa_percentages} \
            --output {output.report_mpa_percentages} \
            >> {log} 2>&1
            '''


    if config["params"]["profiling"]["kraken2"]["krona"]["do"]:
        rule profiling_kraken2_krona_report:
            input:
                taxonomy = config["params"]["profiling"]["kraken2"]["taxonomy"],
                report = expand(
                    os.path.join(
                        config["output"]["profiling"],
                        "profile/kraken2/{sample}/{sample}.kraken2.report"),
                    sample=SAMPLES_ID_LIST)
            output:
                os.path.join(
                    config["output"]["profiling"],
                    "report/kraken2/kraken2_krona.all.html")
            conda:
                config["envs"]["kraken2"]
            log:
                os.path.join(config["output"]["profiling"], "logs/krona/krona_report.log")
            threads:
                1
            resources:
                mem_mb=config["params"]["profiling"]["kraken2"]["mem_mb"]
            shell:
                '''
                taxtab={input.taxonomy}/taxonomy.tab
                if [ ! -e $taxtab ];
                then
                    kronadir=$(dirname $(realpath $(which ktUpdateTaxonomy.sh)))
                    extaxpl=$kronadir/scripts/extractTaxonomy.pl

                    echo "No $taxtab found, run $extaxpl" >{log} 2>&1
                    perl $extaxpl {input.taxonomy} >>{log} 2>&1
                else
                    echo "Found $taxtab." >{log} 2>&1
                fi

                if [ -e $taxtab ];
                then
                    echo "Running ktImportTaxonomy." >>{log} 2>&1
                    ktImportTaxonomy \
                    -t 5 \
                    -m 3 \
                    -tax {input.taxonomy} \
                    {input.report} \
                    -o {output} \
                    >>{log} 2>&1
                else
                    echo "No $taxtab found, run ktImportTaxonomy failed." >> {log} 2>&1
                fi
                '''


        rule profiling_kraken2_all:
            input:
                expand([
                    os.path.join(
                        config["output"]["profiling"],
                        "profile/kraken2/{sample}/{sample}.kraken2.report{suffix}"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/kraken2/kraken2_report.{report}.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/kraken2/kraken2_krona.all.html")],

                        suffix=["", ".mpa.reads_count", ".mpa.percentages"],
                        report=["all", "mpa.reads_count", "mpa.percentages"],
                        sample=SAMPLES_ID_LIST)

    else:
        rule profiling_kraken2_all:
            input:
                expand([
                    os.path.join(
                        config["output"]["profiling"],
                        "profile/kraken2/{sample}/{sample}.kraken2.report{suffix}"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/kraken2/kraken2_report.{report}.tsv")],
                        suffix=["", ".mpa.reads_count", ".mpa.percentages"],
                        report=["all", "mpa.reads_count", "mpa.percentages"],
                        sample=SAMPLES_ID_LIST)


else:
    rule profiling_kraken2_all:
        input:


if config["params"]["profiling"]["kraken2"]["do"] and config["params"]["profiling"]["kraken2"]["bracken"]["do"]:
    rule profiling_kraken2_bracken:
        input:
            os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report")
        output:
            profile = os.path.join(
                config["output"]["profiling"],
                "profile/kraken2_bracken/{sample}/{sample}.bracken.{level}.profile"),
            report = os.path.join(
                config["output"]["profiling"],
                "profile/kraken2_bracken/{sample}/{sample}.bracken.{level}.report")
        log:
            os.path.join(
                config["output"]["profiling"], "logs/kraken2_bracken/{sample}.bracken.{level}.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"], "benchmark/kraken2_bracken/{sample}.bracken.{level}.benchmark.txt")
        params:
            database = config["params"]["profiling"]["kraken2"]["database"],
            reads_len = config["params"]["profiling"]["kraken2"]["bracken"]["reads_len"],
            level = "{level}"
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        resources:
            mem_mb=config["params"]["profiling"]["kraken2"]["mem_mb"]
        conda:
            config["envs"]["kraken2"]
        shell:
            '''
            set +e

            bracken \
            -d {params.database} \
            -i {input} \
            -o {output.profile} \
            -w {output.report} \
            -r {params.reads_len} \
            -l {params.level} \
            -t {threads} \
            > {log} 2>&1


            exitcode=$?
            echo "Exit code is: $exitcode" >> {log}

            if [ $exitcode -eq 1 ];
            then
                grep -oEi "Error: no reads found" {log}
                grepcode=$?
                if [ $grepcode -eq 0 ];
                then
                    echo "Touch {output.report}" >> {log}
                    echo "Touch {output.profile}" >> {log}
                    touch {output.report} >> {log} 2>&1
                    touch {output.profile} >> {log} 2>&1
                    exit 0
                else
                    echo "Runing failed, check kraken2 report please." >> {log} 2>&1
                    exit $exitcode
                fi
            else
                exit $exitcode
            fi
            '''


    rule profiling_kraken2_bracken_merge:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/kraken2_bracken/{sample}/{sample}.bracken.{{level}}.profile"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(
                config["output"]["profiling"],
                "report/kraken2_bracken/bracken.merged.abundance.profile.{level}.tsv")
        priority:
            20
        log:
            os.path.join(
                config["output"]["profiling"], "logs/kraken2_bracken/bracken.merged.{level}.log")
        params:
            samples_id_list = ",".join(SAMPLES_ID_LIST)
        conda:
            config["envs"]["kraken2"]
        threads:
            1
        resources:
            mem_mb=config["params"]["profiling"]["kraken2"]["mem_mb"]
        shell:
            '''
            combine_bracken_outputs.py \
            --files {input} \
            --names {params.samples_id_list} \
            --output {output} \
            > {log} 2>&1
            '''


    rule profiling_kraken2_bracken_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "report/kraken2_bracken/bracken.merged.abundance.profile.{level}.tsv"),
                level=config["params"]["profiling"]["kraken2"]["bracken"]["level"])

else:
    rule profiling_kraken2_bracken_all:
        input:
