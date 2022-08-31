if IS_PE:
    rule profiling_alignment_bowtie2:
        input:
            reads = profiling_input_with_short_reads
        output:
            bam = temp(os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.sorted.bam")),
            flagstat = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.flagstats")
        log:
            os.path.join(config["output"]["profiling"], "logs/bowtie2_samtools/{sample}.bowtie2.log")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/bowtie2_samtools/{sample}.bowtie2.txt")
        threads:
            config["params"]["profiling"]["threads"]
        params:
            bowtie2_db_prefix = config["params"]["profiling"]["genomecov"]["bowtie2_db_prefix"],
            tmp_bam = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.bam.tmp")
        conda:
            config["envs"]["align"]
        shell:
            '''
            rm -rf {output.bam}
            rm -rf {output.flagstat}
            rm -rf {params.tmp_bam}*

            bowtie2 \
            --no-unal -p {threads} -x {params.bowtie2_db_prefix} \
            -1 {input.reads[0]} \
            -2 {input.reads[1]} \
            2> {log} | \
            tee >(samtools flagstat \
                -@{threads} - \
                > {output.flagstat}) | \
            samtools sort -T {params.tmp_bam} -@1 -m 8G -O BAM -o {output.bam} >{log} 2>&1
            '''

else:
    rule profiling_alignment_bowtie2:
        input:
            reads = profiling_input_with_short_reads
        output:
            bam = temp(os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.sorted.bam")),
            flagstat = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.flagstats")
        log:
            os.path.join(config["output"]["profiling"], "logs/bowtie2_samtools/{sample}.bowtie2.log")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/bowtie2_samtools/{sample}.bowtie2.txt")
        threads:
            config["params"]["profiling"]["threads"]
        params:
            bowtie2_db_prefix = config["params"]["profiling"]["genomecov"]["bowtie2_db_prefix"],
            tmp_bam = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.bam.tmp")
        conda:
            config["envs"]["align"]
        shell:
            '''
            rm -rf {output.bam}
            rm -rf {output.flagstat}
            rm -rf {params.tmp_bam}*

            bowtie2 \
            --no-unal -p {threads} -x {params.bowtie2_db_prefix} \
            -U {input.reads} \
            2> {log} | \
            tee >(samtools flagstat \
                 -@{threads} - \
                 > {output.flagstat}) | \
            samtools sort -T {params.tmp_bam} -@1 -m 8G -O BAM -o {output.bam} >{log} 2>&1
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
    conda:
        config["envs"]["align"]
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
            bed = os.path.join(config["output"]["profiling"], "profile/genomecov/{sample}/{sample}.coverage.bed.gz")
        log:
            os.path.join(config["output"]["profiling"], "logs/genomecov/{sample}.genomecov.log")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/genomecov/{sample}.genomecov.txt")
        params:
            bed = os.path.join(config["output"]["profiling"], "profile/genomecov/{sample}/{sample}.coverage.bed")
        threads:
            1
        conda:
            config["envs"]["align"]
        shell:
            '''
            bedtools genomecov -ibam {input.bam} > {params.bed} 2> {log}
            gzip {params.bed} 2>> {log}
            '''


    rule profiling_genomecov_gen_cov:
        input:
            bowtie2_db_fasta = config["params"]["profiling"]["genomecov"]["bowtie2_db_fasta"],
            bed = os.path.join(config["output"]["profiling"], "profile/genomecov/{sample}/{sample}.coverage.bed.gz")
        output:
            coverage = os.path.join(config["output"]["profiling"], "profile/genomecov/{sample}/{sample}.coverage.tsv.gz")
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
            | gzip -c > {output.coverage} 2> {log}
            '''


    rule profiling_genomecov_gen_cov_merge:
        input:
            expand(os.path.join(config["output"]["profiling"],
                                "profile/genomecov/{sample}/{sample}.coverage.tsv.gz"),
                                sample=SAMPLES_ID_LIST)
        output:
            report_cov = os.path.join(config["output"]["profiling"], "report/genomecov/contigs_coverage.cov.tsv.gz"), 
            report_per = os.path.join(config["output"]["profiling"], "report/genomecov/contigs_coverage.per.tsv.gz")
        threads:
            config["params"]["profiling"]["threads"]
        run:
            quantpi.genomecov_merge(input, threads, output_cov=output.report_cov, output_per=output.report_per)


    rule profiling_genomecov_all:
        input:
            expand(os.path.join(config["output"]["profiling"],
                                "report/genomecov/contigs_coverage.{suffix}.tsv.gz"),
                                suffix=["cov", "per"])

else:
    rule profiling_genomecov_all:
        input: