if config["params"]["profiling"]["coverm"]["do"]:
    rule profiling_genome_coverm:
        input:
            genome_dir = config["params"]["profiling"]["coverm"]["genome_dir"],
            bam = os.path.join(config["output"]["profiling"], "align/bowtie2/{sample}/{sample}.sorted.uniq.bam"),
        output:
            table = os.path.join(config["output"]["profiling"], "profile/coverm/{sample}/{sample}.coverm.tsv.gz")
        log:
            os.path.join(config["output"]["profiling"], "logs/coverm/{sample}.coverm.log")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/coverm/{sample}.coverm.txt")
        params:
            methods = " ".join(config["params"]["profiling"]["coverm"]["methods"]),
            genome_fasta_extension = config["params"]["profiling"]["coverm"]["genome_fasta_extension"],
            table = os.path.join(config["output"]["profiling"], "profile/coverm/{sample}/{sample}.coverm.tsv"),
            outdir = os.path.join(config["output"]["profiling"], "profile/coverm/{sample}"),
            min_covered_fraction = config["params"]["profiling"]["coverm"]["min_covered_fraction"],
            contig_end_exclusion = config["params"]["profiling"]["coverm"]["contig_end_exclusion"],
            trim_min = config["params"]["profiling"]["coverm"]["trim_min"],
            trim_max = config["params"]["profiling"]["coverm"]["trim_max"]
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            rm -rf {params.outdir}
            mkdir -p {params.outdir}

            coverm genome \
            --threads {threads} \
            --bam-files {input.bam} \
            --genome-fasta-directory {input.genome_dir} \
            --genome-fasta-extension {params.genome_fasta_extension} \
            --methods {params.methods} \
            --min-covered-fraction {params.min_covered_fraction} \
            --contig-end-exclusion {params.contig_end_exclusion} \
            --trim-min {params.trim_min} \
            --trim-max {params.trim_max} \
            --output-file {params.table} \
            > {log} 2>&1

            gzip {params.table} 2>> {log}
            '''


    rule profiling_genome_coverm_merge:
        input:
            expand(os.path.join(config["output"]["profiling"], "profile/coverm/{sample}/{sample}.coverm.tsv.gz"),
                   sample=SAMPLES_ID_LIST)
        output:
            expand(os.path.join(config["output"]["profiling"], "report/coverm/coverm.{method}.tsv.gz"),
                   method=config["params"]["profiling"]["coverm"]["methods"])
        threads:
            config["params"]["profiling"]["threads"]
        params:
            methods = config["params"]["profiling"]["coverm"]["methods"],
        run:
            quantpi.coverm_merge(input, params.methods, threads, output)


    rule profiling_genome_coverm_all:
        input:
           expand(os.path.join(config["output"]["profiling"], "report/coverm/coverm.{method}.tsv.gz"),
                  method=config["params"]["profiling"]["coverm"]["methods"])
  
 
else:
    rule profiling_genome_coverm_all:
        input: