if config["params"]["trimming"]["sickle"]["do"]:
    rule trimming_sickle:
        input:
            reads = lambda wildcards: get_reads(wildcards, "raw", False)
        output:
            done = touch(os.path.join(config["output"]["trimming"], "short_reads/{sample}/{sample}.sickle.done")),
            reads = expand(
                os.path.join(
                    config["output"]["trimming"],
                    "short_reads/{{sample}}/{{sample}}.trimming{read}.fq.gz"),
                read=[".1", ".2", ".single"] if IS_PE else "") \
                if config["params"]["trimming"]["save_reads"] else \
                temp(expand(os.path.join(
                    config["output"]["trimming"],
                    "short_reads/{{sample}}/{{sample}}.trimming{read}.fq.gz"),
                    read=[".1", ".2", ".single"] if IS_PE else ""))
        log:
            os.path.join(config["output"]["trimming"], "logs/{sample}.sickle.log")
        benchmark:
            os.path.join(config["output"]["trimming"], "benchmark/sickle/{sample}.sickle.benchmark.txt")
        conda:
            config["envs"]["trimming"]
        params:
            output_prefix = os.path.join(config["output"]["trimming"], "short_reads/{sample}/{sample}"),
            quality_type = config["params"]["trimming"]["sickle"]["quality_type"],
            quality_cutoff = config["params"]["trimming"]["sickle"]["quality_cutoff"],
            length_cutoff = config["params"]["trimming"]["sickle"]["length_cutoff"],
            pe = "pe" if IS_PE else "se"
        resources:
            mem_mb=config["params"]["trimming"]["sickle"]["mem_mb"]
        priority:
            10
        shell:
            '''
            READSIN=({input.reads})
            READSOUT=({output.reads})

            if [ "{params.pe}" == "pe" ]; then
                sickle pe \
                --pe-file1 ${{READSIN[0]}} \
                --pe-file2 ${{READSIN[1]}} \
                --output-pe1 ${{READSOUT[0]}} \
                --output-pe2 ${{READSOUT[1]}} \
                --output-single ${{READSOUT[2]}} \
                --qual-type {params.quality_type} \
                --qual-threshold {params.quality_cutoff} \
                --length-threshold {params.length_cutoff} \
                --gzip-output 2> {log}
            else
                sickle se \
                --fastq-file ${{READSIN[0]}} \
                --output-file ${{READSOUT[0]}} \
                --qual-type {params.quality_type} \
                --qual-threshold {params.quality_cutoff} \
                --length-threshold {params.length_cutoff} \
                --gzip-output 2> {log}
            fi
            '''


    rule trimming_sickle_all:
        input:
            expand(
                os.path.join(config["output"]["trimming"], "short_reads/{sample}/{sample}.sickle.done"),
                sample=SAMPLES_ID_LIST)

else:
    rule trimming_sickle_all:
        input:


ADAPTER_OPERATION = ""
adapter_sequence = config["params"]["trimming"]["fastp"]["adapter_sequence"]
adapter_sequence_r2 = config["params"]["trimming"]["fastp"]["adapter_sequence_r2"]

if config["params"]["trimming"]["fastp"]["disable_adapter_trimming"]:
    adapter_operation = "--disable_adapter_trimming"
else:
    if IS_PE:
        if config["params"]["trimming"]["fastp"]["detect_adapter_for_pe"]:
            ADAPTER_OPERATION = "--detect_adapter_for_pe"
        else:
            ADAPTER_OPERATION = f'''--adapter_sequence {adapter_sequence} --adapter_sequence_r2 {adapter_sequence_r2}'''
    else:
        if config["params"]["trimming"]["fastp"]["detect_adapter_for_se"]:
            ADAPTER_OPERATION = ""
        else:
            ADAPTER_OPERATION = f'''--adapter_sequence {adapter_sequence}'''


if config["params"]["trimming"]["fastp"]["do"]:
    rule trimming_fastp:
        input:
            reads = lambda wildcards: get_reads(wildcards, "raw")
        output:
            html = os.path.join(config["output"]["trimming"],
                                "short_reads/{sample}/{sample}.fastp.html"),
            json = os.path.join(config["output"]["trimming"],
                                "short_reads/{sample}/{sample}.fastp.json"),
            reads = expand(
                os.path.join(
                    config["output"]["trimming"],
                    "short_reads/{{sample}}/{{sample}}.trimming{read}.fq.gz"),
                read=[".1", ".2"] if IS_PE else "") \
                if config["params"]["trimming"]["save_reads"] else \
                temp(expand(os.path.join(
                    config["output"]["trimming"],
                    "short_reads/{{sample}}/{{sample}}.trimming{read}.fq.gz"),
                    read=[".1", ".2"] if IS_PE else ""))
        params:
            output_prefix = os.path.join(config["output"]["trimming"], "short_reads/{sample}/{sample}"),
            compression = config["params"]["trimming"]["fastp"]["compression"],
            cut_front_window_size = config["params"]["trimming"]["fastp"]["cut_front_window_size"],
            cut_front_mean_quality = config["params"]["trimming"]["fastp"]["cut_front_mean_quality"],
            cut_tail_window_size = config["params"]["trimming"]["fastp"]["cut_tail_window_size"],
            cut_tail_mean_quality = config["params"]["trimming"]["fastp"]["cut_tail_mean_quality"],
            cut_right_window_size = config["params"]["trimming"]["fastp"]["cut_right_window_size"],
            cut_right_mean_quality = config["params"]["trimming"]["fastp"]["cut_right_mean_quality"],
            length_required = config["params"]["trimming"]["fastp"]["length_required"],
            n_base_limit = config["params"]["trimming"]["fastp"]["n_base_limit"],
            use_slide_window = "yes" if config["params"]["trimming"]["fastp"]["use_slide_window"] else "no",
            pe = "pe" if IS_PE else "se",
            dedup = f'''--dedup --dup_calc_accuracy {config["params"]["trimming"]["fastp"]["dup_calc_accuracy"]}''' \
                if config["params"]["trimming"]["fastp"]["dedup"] else ""
        log:
            os.path.join(config["output"]["trimming"], "logs/{sample}.fastp.log")
        benchmark:
            os.path.join(config["output"]["trimming"], "benchmark/fastp/{sample}.fastp.benchmark.txt")
        conda:
            config["envs"]["trimming"]
        threads:
            config["params"]["trimming"]["threads"]
        resources:
            mem_mb=config["params"]["trimming"]["fastp"]["mem_mb"]
        priority:
            10
        shell:
            '''
            READSIN=({input.reads})
            READSOUT=({output.reads})

            if [ "{params.pe}" == "pe" ]; then
                if [ "{params.use_slide_window}" == "yes" ]; then
                    fastp \
                    --in1 ${{READSIN[0]}} \
                    --in2 ${{READSIN[1]}} \
                    --out1 ${{READSOUT[0]}} \
                    --out2 ${{READSOUT[1]}} \
                    --compression {params.compression} \
                    {ADAPTER_OPERATION} \
                    {params.dedup} \
                    --cut_front \
                    --cut_right \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_right_window_size {params.cut_right_window_size} \
                    --cut_right_mean_quality {params.cut_right_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html {output.html} \
                    --json {output.json} 2> {log}
                else
                    fastp \
                    --in1 ${{READSIN[0]}} \
                    --in2 ${{READSIN[1]}} \
                    --out1 ${{READSOUT[0]}} \
                    --out2 ${{READSOUT[1]}} \
                    --compression {params.compression} \
                    {ADAPTER_OPERATION} \
                    {params.dedup} \
                    --cut_front \
                    --cut_tail \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_tail_window_size {params.cut_tail_window_size} \
                    --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html {output.html} \
                    --json {output.json} 2> {log}
                fi
            else
                if [ "{params.use_slide_window}" == "yes" ]; then
                    fastp \
                    --in1 ${{READSIN[0]}} \
                    --out1 ${{READSOUT[0]}} \
                    --compression {params.compression} \
                    {ADAPTER_OPERATION} \
                    {params.dedup} \
                    --cut_front \
                    --cut_right \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_right_window_size {params.cut_right_window_size} \
                    --cut_right_mean_quality {params.cut_right_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html {output.html} \
                    --json {output.json} 2> {log}
                else
                    fastp \
                    --in1 ${{READSIN[0]}} \
                    --out1 ${{READSOUT[0]}} \
                    --compression {params.compression} \
                    {ADAPTER_OPERATION} \
                    {params.dedup} \
                    --cut_front \
                    --cut_tail \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_tail_window_size {params.cut_tail_window_size} \
                    --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html {output.html} \
                    --json {output.json} 2> {log}
                fi
            fi
            '''


    rule trimming_fastp_multiqc:
        input:
            expand(
                os.path.join(config["output"]["trimming"], "short_reads/{sample}/{sample}.fastp.json"),
                sample=SAMPLES_ID_LIST)
        output:
            html = os.path.join(
                config["output"]["trimming"],
                "report/fastp_multiqc_report.html")
        log:
            os.path.join(config["output"]["trimming"], "logs/multiqc.fastp.log")
        params:
            outdir = os.path.join(config["output"]["trimming"], "report")
        conda:
            config["envs"]["multiqc"]
        threads:
            1
        resources:
            mem_mb=config["params"]["qcreport"]["mem_mb"]
        priority:
            10
        shell:
            '''
            multiqc \
            --outdir {params.outdir} \
            --title fastp \
            --module fastp \
            {input} \
            2> {log}
            '''


    rule trimming_fastp_all:
        input:
            expand([
                os.path.join(config["output"]["trimming"], "short_reads/{sample}/{sample}.fastp.html"),
                os.path.join(config["output"]["trimming"], "short_reads/{sample}/{sample}.fastp.json"),
                os.path.join(config["output"]["trimming"], "report/fastp_multiqc_report.html")
            ], sample=SAMPLES_ID_LIST)

else:
    rule trimming_fastp_all:
        input:


if config["params"]["trimming"]["trimmomatic"]["do"]:
    rule trimming_trimmomatic:
        input:
            reads = lambda wildcards: get_reads(wildcards, "raw")
        output:
            summary = os.path.join(config["output"]["trimming"], "short_reads/{sample}/{sample}.trimmomatic.summary.txt"),
            log = os.path.join(config["output"]["trimming"], "short_reads/{sample}/{sample}.trimmomatic.log"),
            reads = expand(
                os.path.join(
                    config["output"]["trimming"],
                    "short_reads/{{sample}}/{{sample}}.trimming{read}.fq.gz"),
                read=[".1", ".2"] if IS_PE else "") \
                if config["params"]["trimming"]["save_reads"] else \
                temp(expand(os.path.join(
                    config["output"]["trimming"],
                    "short_reads/{{sample}}/{{sample}}.trimming{read}.fq.gz"),
                    read=[".1", ".2"] if IS_PE else ""))
        params:
            trimmomatic_options = config["params"]["trimming"]["trimmomatic"]["trimmomatic_options"],
            phred = config["params"]["trimming"]["trimmomatic"]["phred"],
            save_unpaired = "true" if config["params"]["trimming"]["trimmomatic"]["save_unpaired"] else "false",
            pe = "pe" if IS_PE else "se"
        log:
            os.path.join(config["output"]["trimming"], "logs/{sample}.trimmomatic.log")
        benchmark:
            os.path.join(
                config["output"]["trimming"],
                "benchmark/trimmomatic/{sample}.trimmomatic.benchmark.txt")
        conda:
            config["envs"]["trimming"]
        threads:
            config["params"]["trimming"]["threads"]
        resources:
            mem_mb=config["params"]["trimming"]["trimmomatic"]["mem_mb"]
        priority:
            10
        shell:
            '''
            READSIN=({input.reads})
            READSOUT=({output.reads})

            if [ "{params.pe}" == "pe" ]; then
                trimmomatic \
                PE \
                {params.phred} \
                -threads {threads} \
                -summary {output.summary} \
                ${{READSIN[0]}} ${{READSIN[1]}} \
                ${{READSOUT[0]}} ${{READSOUT[0]}}.unpaired.gz \
                ${{READSOUT[1]}} ${{READSOUT[1]}}.unpaired.gz \
                {params.trimmomatic_options} \
                >{output.log} 2>&1

                if [ "{params.save_unpaired}" == "false" ]; then
                    rm -rf ${{READSOUT[0]}}.unpaired.gz
                    rm -rf ${{READSOUT[1]}}.unpaired.gz
                fi
            else
                trimmomatic \
                SE \
                {params.phred} \
                -threads {threads} \
                -summary {output.summary} \
                ${{READSIN[0]}} \
                ${{READSOUT[0]}} \
                {params.trimmomatic_options} \
                >{output.log} 2>&1
            fi

            cp {output.log} {log}
            '''


    rule trimming_trimmomatic_multiqc:
        input:
            expand(
                os.path.join(
                    config["output"]["trimming"],
                    "short_reads/{sample}/{sample}.trimmomatic.log"),
                    sample=SAMPLES_ID_LIST)
        output:
            html = os.path.join(config["output"]["trimming"], "report/trimmomatic_multiqc_report.html")
        log:
            os.path.join(config["output"]["trimming"], "logs/multiqc.trimmomatic.log")
        params:
            outdir = os.path.join(config["output"]["trimming"], "report")
        conda:
            config["envs"]["multiqc"]
        threads:
            1
        resources:
            mem_mb=config["params"]["qcreport"]["mem_mb"]
        priority:
            10
        shell:
            '''
            multiqc \
            --outdir {params.outdir} \
            --title trimmomatic \
            --module trimmomatic \
            {input} \
            2> {log}
            '''


    rule trimming_trimmomatic_all:
        input:
            expand([
                os.path.join(
                    config["output"]["trimming"],
                    "short_reads/{sample}/{sample}.trimmomatic.summary.txt"),
                os.path.join(
                    config["output"]["trimming"],
                    "short_reads/{sample}/{sample}.trimmomatic.log"),
                os.path.join(
                    config["output"]["trimming"],
                    "report/trimmomatic_multiqc_report.html")
                ],
                sample=SAMPLES_ID_LIST)

else:
    rule trimming_trimmomatic_all:
        input:


if TRIMMING_DO and config["params"]["qcreport"]["do"]:
    rule trimming_report:
        input:
            lambda wildcards: get_reads(wildcards, "trimming")
        output:
            temp(os.path.join(
                config["output"]["trimming"],
                "report/stats/{sample}_trimming_stats.tsv.raw"))
        log:
            os.path.join(
                config["output"]["trimming"],
                "logs/{sample}.seqkit.log")
        params:
            fq_encoding = config["params"]["fq_encoding"]
        conda:
            config["envs"]["report"]
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        resources:
            mem_mb=config["params"]["qcreport"]["mem_mb"]
        priority:
            10
        shell:
            '''
            seqkit stats \
            --all \
            --basename \
            --tabular \
            --fq-encoding {params.fq_encoding} \
            --out-file {output} \
            --threads {threads} \
            {input} 2> {log}
            '''


    rule trimming_report_refine:
        input:
            os.path.join(
                config["output"]["trimming"],
                "report/stats/{sample}_trimming_stats.tsv.raw")
        output:
            os.path.join(
                config["output"]["trimming"],
                "report/stats/{sample}_trimming_stats.tsv")
        params:
            sample_id = "{sample}"
        threads:
            1
        resources:
            mem_mb=config["params"]["qcreport"]["mem_mb"]
        priority:
            10
        run:
            if IS_PE:
                quantpi.change(
                    str(input), str(output), params.sample_id, "trimming",
                    "pe", ["fq1", "fq2"])
            else:
                quantpi.change(
                    str(input), str(output), params.sample_id, "trimming",
                    "se", ["fq1"])


    rule trimming_report_merge:
        input:
            expand(
                os.path.join(
                    config["output"]["trimming"],
                    "report/stats/{sample}_trimming_stats.tsv"),
                    sample=SAMPLES_ID_LIST)
        output:
            os.path.join(config["output"]["qcreport"], "trimming_stats.tsv")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        resources:
            mem_mb=config["params"]["qcreport"]["mem_mb"] * 5
        run:
            quantpi.merge(input, quantpi.parse, threads, output=output[0])


    rule trimming_report_all:
        input:
            os.path.join(config["output"]["qcreport"], "trimming_stats.tsv")

else:
    rule trimming_report_all:
        input:


rule trimming_all:
    input:
        rules.trimming_sickle_all.input,
        rules.trimming_fastp_all.input,
        rules.trimming_trimmomatic_all.input,
        rules.trimming_report_all.input


localrules:
    trimming_fastp_all,
    trimming_sickle_all,
    trimming_trimmomatic_all,
    trimming_report_all,
    trimming_all