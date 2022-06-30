rule profiling_alignment_bowtie2:
    input:
        reads = profiling_input_with_short_reads
    output:
        bam = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.sorted.bam"),
        flagstat = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.flagstats")
    log:
        os.path.join(config["output"]["profiling"], "logs/bowtie2_samtools/{sample}.bowtie2.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/bowtie2_samtools/{sample}.bowtie2.txt")
    threads:
        config["params"]["profiling"]["threads"]
    params:
        bowtie2_db_prefix = config["params"]["profiling"]["genomecov"]["bowtie2_db_prefix"],
        reads_layout = 1 if IS_PE else 0,
        tmp_bam = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.bam.tmp")
    shell:
        '''
        rm -rf {output.bam}
        rm -rf {output.flagstat}
        rm -rf {params.tmp_bam}*

        if [ {params.reads_layout} -eq 1 ]
        then
            bowtie2 \
            --no-unal -p {threads} -x {params.bowtie2_db_prefix} \
            -1 {input.reads[0]} \
            -2 {input.reads[1]} \
            2> {log} | \
            tee >(samtools flagstat \
                  -@{threads} - \
                  > {output.flagstat}) | \
            samtools sort -T {params.tmp_bam} -@1 -m 8G -O BAM -o {output.bam} >{log} 2>&1
        else
            bowtie2 \
            --no-unal -p {threads} -x {params.bowtie2_db_prefix} \
            -U {input.reads} \
            2> {log} | \
            tee >(samtools flagstat \
                  -@{threads} - \
                  > {output.flagstat}) | \
            samtools sort -T {params.tmp_bam} -@1 -m 8G -O BAM -o {output.bam} >{log} 2>&1
        fi
        '''


rule profiling_alignment_bam_postprocess:
    input:
        bam = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.sorted.bam")
    output:
        bam = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.sorted.uniq.bam")
    threads:
        config["params"]["profiling"]["threads"]
    log:
        os.path.join(config["output"]["profiling"], "logs/sambamba/{sample}.sambamba.log")
    benchmark:
        os.path.join(config["output"]["profiling"], "benchmark/sambamba/{sample}.sambamba.txt")
    shell:
        '''
        sambamba view \
        -h -t {threads} -f bam \
        -F "[XS] == null and not unmapped  and not duplicate" \
        {input.bam} -o {output.bam} >{log} 2>&1
        '''


if config["params"]["profiling"]["genomecov"]["do"]:
    rule profiling_genomecov_gen_bed:
        input:
            bam = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.sorted.uniq.bam")
        output:
            bed = os.path.join(config["output"]["profiling"], "profile/genomecov/{sample}/{sample}.coverage.bed")
        log:
            os.path.join(config["output"]["profiling"], "logs/genomecov/{sample}.genomecov.log")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/genomecov/{sample}.genomecov.txt")
        threads:
            1
        shell:
            '''
            bedtools genomecov -ibam {input.bam} > {output.bed} 2> {log}
            '''


    rule profiling_genomecov_gen_cov:
        input:
            bowtie2_db_fasta = config["params"]["profiling"]["genomecov"]["bowtie2_db_fasta"],
            bed = os.path.join(config["output"]["profiling"], "profile/genomecov/{sample}/{sample}.coverage.bed")
        output:
            coverage = os.path.join(config["output"]["profiling"], "profile/genomecov/{sample}/{sample}.coverage.tsv")
        params:
            gen_contig_cov_script = config["params"]["profiling"]["genomecov"]["gen_contig_cov_script"]
        log:
            os.path.join(config["output"]["profiling"], "logs/gen_contig_cov/{sample}.gen_contig_cov.log")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/gen_contig_cov/{sample}.gen_contig_cov.txt")
        threads:
            1
        shell:
            '''
            python {params.gen_contig_cov_script} \
            --isbedfiles \
            {input.bowtie2_db_fasta} \
            {input.bed} \
            > {output.coverage} 2> {log}
            '''


    rule profiling_genomecov_gen_cov_merge:
        input:
            expand(os.path.join(config["output"]["profiling"],
                                "profile/genomecov/{sample}/{sample}.coverage.tsv"),
                                sample=SAMPLES_ID_LIST)
        output:
            report_cov = os.path.join(config["output"]["profiling"], "report/genomecov/contigs_coverage.cov.tsv"), 
            report_per = os.path.join(config["output"]["profiling"], "report/genomecov/contigs_coverage.per.tsv")
        threads:
            config["params"]["profiling"]["threads"]
        run:
            metapi.genomecov_merge(input, threads, output_cov=output.report_cov, output_per=output.report_per)


    rule profiling_genomecov_all:
        input:
            expand(os.path.join(config["output"]["profiling"],
                                "report/genomecov/contigs_coverage.{suffix}"),
                                suffix=["cov.tsv", "per.tsv"])

else:
    rule profiling_genomecov_all:
        input: