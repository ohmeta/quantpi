if config["params"]["profiling"]["krakenuniq"]["do"]:
    rule profiling_krakenuniq_preload_database:
        input:
            database = config["params"]["profiling"]["krakenuniq"]["database"]
        output:
            done = os.path.join(config["output"]["profiling"], "config/krakenuniq/db_preload_done")
        log:
            os.path.join(
                config["output"]["profiling"], "logs/krakenuniq/krakenuniq_preload_database.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"], "benchmark/krakenuniq/krakenuniq_preload_database.benchmark.txt")
        conda:
            config["envs"]["krakenuniq"]
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            krakenuniq --db {input.database} --preload --threads {threads} >{log} 2>&1
            touch {output.done}
            '''


    rule profiling_krakenuniq:
        input:
            done = os.path.join(config["output"]["profiling"], "config/krakenuniq/db_preload_done"),
            reads = profiling_input_with_short_reads
        output:
            report = os.path.join(
                config["output"]["profiling"],
                "profile/krakenuniq/{sample}/{sample}.krakenuniq.report"),
            report_mpa_reads_count = os.path.join(
                config["output"]["profiling"],
                "profile/krakenuniq/{sample}/{sample}.krakenuniq.report.mpa.reads_count"),
            report_mpa_percentages = os.path.join(
                config["output"]["profiling"],
                "profile/krakenuniq/{sample}/{sample}.krakenuniq.report.mpa.percentages")
        log:
            os.path.join(
                config["output"]["profiling"], "logs/krakenuniq/{sample}.krakenuniq.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"], "benchmark/krakenuniq/{sample}.krakenuniq.benchmark.txt")
        params:
            save_table = config["params"]["profiling"]["krakenuniq"]["save_table"],
            paired = "--paired" if IS_PE else "",
            check_names = "--check-names" if IS_PE else "",
            database = config["params"]["profiling"]["krakenuniq"]["database"],
            hll_precision = config["params"]["profiling"]["krakenuniq"]["hll_precision"],
            exact = "--exact" if config["params"]["profiling"]["krakenuniq"]["exact"] else "",
            quick = "--quick" \
                if config["params"]["profiling"]["krakenuniq"]["quick"] \
                else "",
            unclassified_out = "--unclassified-out %s" % \
                os.path.join(
                    config["output"]["profiling"],
                    "profile/krakenuniq/{sample}/{sample}.krakenuniq.unclassified%s.fq" \
                    % "#" if IS_PE else "") \
                    if config["params"]["profiling"]["krakenuniq"]["unclassified_out"] \
                    else "",
            classified_out = "--classified-out %s" % \
                os.path.join(
                    config["output"]["profiling"],
                    "profile/krakenuniq/{sample}/{sample}.krakenuniq.classified%s.fq" \
                    % "#" if IS_PE else "") \
                    if config["params"]["profiling"]["krakenuniq"]["classified_out"] \
                    else "",
            table = "--output %s" % \
                os.path.join(
                    config["output"]["profiling"],
                    "profile/krakenuniq/{sample}/{sample}.krakenuniq.table") \
                if config["params"]["profiling"]["krakenuniq"]["save_table"] \
                    else "--output /dev/null"
        threads:
            config["params"]["profiling"]["threads"]
        conda:
            config["envs"]["krakenuniq"]
        shell:
            '''
            krakenuniq \
            --report-file {output.report} \
            --threads {threads} \
            --db {params.database} \
            {params.quick} \
            --hll-precision {params.hll_precision} \
            {params.exact} \
            {params.unclassified_out} \
            {params.classified_out} \
            {params.table} \
            {params.paired} \
            {params.check_names} \
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


    rule profiling_krakenuniq_combine_kreport:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/krakenuniq/{sample}/{sample}.krakenuniq.report"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(
                config["output"]["profiling"],
                "report/krakenuniq/krakenuniq_report.all.tsv")
        params:
            samples_name = " ".join(list(SAMPLES_ID_LIST))
        conda:
            config["envs"]["krakenuniq"]
        log:
            os.path.join(config["output"]["profiling"],
                            "logs/krakentools/combine_kreports.log")
        shell:
            '''
            combine_kreports.py \
            --report-file {input} \
            --sample-names {params.samples_name} \
            --display-headers \
            --output {output} \
            > {log} 2>&1
            '''


    rule profiling_krakenuniq_combine_kreport_mpa:
        input:
            report_mpa_reads_count = expand(os.path.join(
                config["output"]["profiling"],
                "profile/krakenuniq/{sample}/{sample}.krakenuniq.report.mpa.reads_count"),
                sample=SAMPLES_ID_LIST),
            report_mpa_percentages = expand(os.path.join(
                config["output"]["profiling"],
                "profile/krakenuniq/{sample}/{sample}.krakenuniq.report.mpa.percentages"),
                sample=SAMPLES_ID_LIST)
        output:
            report_mpa_reads_count = os.path.join(
                config["output"]["profiling"],
                "report/krakenuniq/krakenuniq_report.mpa.reads_count.tsv"),
            report_mpa_percentages = os.path.join(
                config["output"]["profiling"],
                "report/krakenuniq/krakenuniq_report.mpa.percentages.tsv")
        conda:
            config["envs"]["krakenuniq"]
        log:
            os.path.join(
                config["output"]["profiling"], "logs/krakentools/combine_kreports_mpa.log")
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


    if config["params"]["profiling"]["krakenuniq"]["krona"]["do"]:
        rule profiling_krakenuniq_krona_report:
            input:
                taxonomy = config["params"]["profiling"]["krakenuniq"]["taxonomy"],
                report = expand(
                    os.path.join(
                        config["output"]["profiling"],
                        "profile/krakenuniq/{sample}/{sample}.krakenuniq.report"),
                    sample=SAMPLES_ID_LIST)
            output:
                os.path.join(
                    config["output"]["profiling"],
                    "report/krakenuniq/krakenuniq_krona.all.html")
            conda:
                config["envs"]["krakenuniq"]
            log:
                os.path.join(config["output"]["profiling"],
                            "logs/krona/krona_report.log")
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


        rule profiling_krakenuniq_all:
            input:
                expand([
                    os.path.join(
                        config["output"]["profiling"],
                        "profile/krakenuniq/{sample}/{sample}.krakenuniq.report{suffix}"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/krakenuniq/krakenuniq_report.{report}.tsv"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/krakenuniq/krakenuniq_krona.all.html")],

                        suffix=["", ".mpa.reads_count", ".mpa.percentages"],
                        report=["all", "mpa.reads_count", "mpa.percentages"],
                        sample=SAMPLES_ID_LIST),
                rules.qcreport_all.input

    else:
        rule profiling_krakenuniq_all:
            input:
                expand([
                    os.path.join(
                        config["output"]["profiling"],
                        "profile/krakenuniq/{sample}/{sample}.krakenuniq.report{suffix}"),
                    os.path.join(
                        config["output"]["profiling"],
                        "report/krakenuniq/krakenuniq_report.{report}.tsv")],
                        suffix=["", ".mpa.reads_count", ".mpa.percentages"],
                        report=["all", "mpa.reads_count", "mpa.percentages"],
                        sample=SAMPLES_ID_LIST),
                rules.qcreport_all.input


else:
    rule profiling_krakenuniq_all:
        input:


if config["params"]["profiling"]["krakenuniq"]["do"] and config["params"]["profiling"]["krakenuniq"]["bracken"]["do"]:
    rule profiling_krakenuniq_bracken:
        input:
            os.path.join(
                config["output"]["profiling"],
                "profile/krakenuniq/{sample}/{sample}.krakenuniq.report")
        output:
            profile = os.path.join(
                config["output"]["profiling"],
                "profile/krakenuniq_bracken/{sample}/{sample}.bracken.{level}.profile"),
            report = os.path.join(
                config["output"]["profiling"],
                "profile/krakenuniq_bracken/{sample}/{sample}.bracken.{level}.report")
        log:
            os.path.join(
                config["output"]["profiling"], "logs/krakenuniq_bracken/{sample}.bracken.{level}.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"], "benchmark/krakenuniq_bracken/{sample}.bracken.{level}.benchmark.txt")
        params:
            database = config["params"]["profiling"]["krakenuniq"]["database"],
            reads_len = config["params"]["profiling"]["krakenuniq"]["bracken"]["reads_len"],
            level = "{level}"
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        conda:
            config["envs"]["krakenuniq"]
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
                    echo "Runing failed, check krakenuniq report please." >> {log} 2>&1
                    exit $exitcode
                fi
            else
                exit $exitcode
            fi
            '''


    rule profiling_krakenuniq_bracken_merge:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/krakenuniq_bracken/{sample}/{sample}.bracken.{{level}}.profile"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(
                config["output"]["profiling"],
                "report/krakenuniq_bracken/bracken.merged.abundance.profile.{level}.tsv")
        priority:
            20
        log:
            os.path.join(
                config["output"]["profiling"], "logs/krakenuniq_bracken/bracken.merged.{level}.log")
        params:
            samples_id_list = ",".join(SAMPLES_ID_LIST)
        conda:
            config["envs"]["krakenuniq"]
        shell:
            '''
            combine_bracken_outputs.py \
            --files {input} \
            --names {params.samples_id_list} \
            --output {output} \
            > {log} 2>&1
            '''


    rule profiling_krakenuniq_bracken_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "report/krakenuniq_bracken/bracken.merged.abundance.profile.{level}.tsv"),
                level=config["params"]["profiling"]["krakenuniq"]["bracken"]["level"]),

            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_krakenuniq_bracken_all:
        input:
