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
            report = protected(os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report")),
            report_mpa_reads_count = protected(os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report.mpa.reads_count")),
            report_mpa_percentages = protected(os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report.mpa.percentages"))
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/kraken2/{sample}.kraken2.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/kraken2/{sample}.kraken2.benchmark.txt")
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
                    else "",
        threads:
            config["params"]["profiling"]["threads"]
        run:
            shell(
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
                ''')

            shell(
                '''
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
                ''')

            if params.save_table:
                shell('''pigz %s''' % params.table.split(" ")[-1])


    rule profiling_kraken2_krona_report:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/kraken2/{sample}/{sample}.kraken2.report"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(
                config["output"]["profiling"],
                "report/kraken2/kraken2_krona.all.html")
        shell:
            '''
            ktImportTaxonomy -q 2 -t 3 {input} -o {output}
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
        shell:
            '''
            combine_kreports.py \
            --report-file {input} \
            --sample-names {params.samples_name} \
            --display-headers \
            --output {output}
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
        shell:
            '''
            combine_mpa.py \
            --input {input.report_mpa_reads_count} \
            --output {output.report_mpa_reads_count}

            combine_mpa.py \
            --input {input.report_mpa_percentages} \
            --output {output.report_mpa_percentages}
            '''


    rule profiling_kraken2_all:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/kraken2/{sample}/{sample}.kraken2.report{suffix}"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/kraken2/kraken2_krona.all.html"),
                os.path.join(
                    config["output"]["profiling"],
                    "report/kraken2/kraken2_report.{report}.tsv")
                    ],
                    suffix=["", ".mpa.reads_count", ".mpa.percentages"],
                    report=["all", "mpa.reads_count", "mpa.percentages"],
                    sample=SAMPLES_ID_LIST),

            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_short_reads_kraken2_all:
        input:


if config["params"]["profiling"]["kraken2"]["do"] and \
   config["params"]["profiling"]["bracken"]["do"]:
    rule profiling_bracken:
        input:
            os.path.join(
                config["output"]["profiling"],
                "profile/kraken2/{sample}/{sample}.kraken2.report")
        output:
            profile = protected(os.path.join(
                config["output"]["profiling"],
                "profile/bracken/{sample}/{sample}.bracken.{level}.profile")),
            report = protected(os.path.join(
                config["output"]["profiling"],
                "profile/bracken/{sample}/{sample}.bracken.{level}.report"))
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/bracken/{sample}.bracken.{level}.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/bracken/{sample}.bracken.{level}.benchmark.txt")
        params:
            database = config["params"]["profiling"]["kraken2"]["database"],
            reads_len = config["params"]["profiling"]["bracken"]["reads_len"],
            level = "{level}"
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            bracken \
            -d {params.database} \
            -i {input} \
            -o {output.profile} \
            -w {output.report} \
            -r {params.reads_len} \
            -l {params.level} \
            -t {threads} \
            > {log} 2>&1
            '''


    rule profiling_bracken_merge:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/bracken/{sample}/{sample}.bracken.{{level}}.profile"),
                 sample=SAMPLES_ID_LIST)
        output:
            os.path.join(
                config["output"]["profiling"],
                "report/bracken/bracken.merged.abundance.profile.{level}.tsv")
        priority:
            20
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/bracken/bracken.merged.{level}.log")
        run:
            shell(
                '''
                combine_bracken_outputs.py \
                --files {input} \
                --names %s \
                --output {output} \
                > {log} 2>&1
                ''' % ",".join(SAMPLES_ID_LIST))


    rule profiling_bracken_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "report/bracken/bracken.merged.abundance.profile.{level}.tsv"),
                   level=config["params"]["profiling"]["bracken"]["level"]),
           
            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_bracken_all:
        input: