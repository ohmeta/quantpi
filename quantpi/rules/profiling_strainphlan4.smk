## reference
## https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4


if config["params"]["profiling"]["strainphlan"]["do_v4"]:
    rule profiling_strainphlan4_sample2markers:
        input:
            database_pkl = expand(os.path.join(
                config["params"]["profiling"]["metaphlan"]["bowtie2db"], "{index}.pkl"),
                index = config["params"]["profiling"]["metaphlan"]["index_v4"]),
            sam = os.path.join(config["output"]["profiling"],
                               "profile/metaphlan4/{sample}/{sample}.sam.bz2"),
            aln = os.path.join(config["output"]["profiling"],
                               "profile/metaphlan4/{sample}/{sample}.bowtie2.bz2")
        output:
            os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan4/{sample}/consensus_markers/{sample}.pkl")
        log:
            os.path.join(
                config["output"]["profiling"],
                "logs/strainphlan4_sample2markers/{sample}.strainphlan4_sample2markers.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"],
                "benchmark/strainphlan4_sample2markers/{sample}.strainphlan4_sample2markers.benchmark.txt")
        conda:
            config["envs"]["biobakery4"]
        params:
            outdir = os.path.join(config["output"]["profiling"], "profile/strainphlan4/{sample}/consensus_markers")
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            mkdir -p {params.outdir}

            sample2markers.py \
            --database {input.database_pkl} \
            --input {input.sam} \
            --output_dir {params.outdir}/ \
            --nprocs {threads} \
            > {log} 2>&1
            '''


    STRAINPHLAN_CLADES_V4 = \
        pd.read_csv(config["params"]["profiling"]["strainphlan"]["clades_tsv_v4"], sep="\t")\
          .set_index("clades")
    STRAINPHLAN_CLADES_LIST_V4 = STRAINPHLAN_CLADES_V4.index.unique()


    rule profiling_strainphlan4_extract_markers:
        input:
            clades_tsv = config["params"]["profiling"]["strainphlan"]["clades_tsv"],
            database_pkl = expand(os.path.join(
                config["params"]["profiling"]["metaphlan"]["bowtie2db"], "{index}.pkl"),
                index = config["params"]["profiling"]["metaphlan"]["index_v4"])
        output:
            clade_marker = os.path.join(
                config["output"]["profiling"],
                "databases/strainphlan4/clade_markers/{clade}.fna")
        log:
            os.path.join(
                config["output"]["profiling"], 
                "logs/strainphlan4_extract_markers/{clade}.strainphlan4_extract_markers.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"],
                "benchmark/strainphlan4_extract_markers/{clade}.strainphlan4_extract_markers.benchmark.txt")
        conda:
            config["envs"]["biobakery4"]
        params:
            clade = "{clade}",
            outdir = os.path.join(config["output"]["profiling"], "databases/strainphlan4/clade_markers")
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            mkdir -p {params.outdir}

            extract_markers.py \
            --database {input.database_pkl} \
            --clade {params.clade} \
            --output_dir {params.outdir}/ \
            > {log} 2>&1
            '''


    rule profiling_strainphlan4:
        input:
            clade_marker = os.path.join(
                config["output"]["profiling"],
                "databases/strainphlan4/clade_markers/{clade}.fna")
            consensus_markers = expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan4/{sample}/consensus_markers/{sample}.pkl"),
                sample=SAMPLES_ID_LIST),
            reference_genome = lambda wildcards: STRAINPHLAN_CLADES_V3.loc[wildcards.clade, "fna_path"]
        output:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan4/{{sample}}/{{clade}}/RAxML_{prefix}.{{clade}}.StrainPhlAn4.tre"),
                prefix=["bestTree", "info", "log", "parsimonyTree", "result"]),
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan4/{{sample}}/{{clade}}/{{clade}}{suffix}"),
                suffix=["_mutation_rates", ".info", ".mutation", ".polymorphic",
                        ".StrainPhlAn4_concatenated.aln"])
        log:
            os.path.join(
                config["output"]["profiling"],
                "logs/strainphlan4/{sample}.{clade}.strainphlan4.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"],
                "benchmark/strainphlan4/{sample}.{clade}.strainphlan4.benchmark.txt")
        conda:
            config["envs"]["biobakery4"]
        params:
            clade = "{clade}",
            outdir = os.path.join(config["output"]["profiling"], "profile/strainphlan4/{sample}/{clade}"),
            opts = config["params"]["profiling"]["strainphlan"]["external_opts_v4"]
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            mkdir -p {params.outdir}

            strainphlan \
            --samples {input.consensus_markers} \
            --clade_markers {input.clade_marker} \
            --references {input.reference_genome} \
            --output_dir {params.outdir}/ \
            --nprocs {threads} \
            --clade {params.clade} \
            --mutation_rates \
            {params.opts} \
            >{log} 2>&1
            '''


    rule profiling_strainphlan4_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan4/{sample}/{clade}/RAxML_{prefix}.{clade}.StrainPhlAn4.tre"),
                prefix=["bestTree", "info", "log", "parsimonyTree", "result"],
                sample=SAMPLES_ID_LIST,
                clade=STRAINPHLAN_CLADES_LIST_V4),
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan4/{sample}/{clade}/{clade}{suffix}"),
                suffix=["_mutation_rates", ".info", ".mutation", ".polymorphic",
                        ".StrainPhlAn4_concatenated.aln"],
                sample=SAMPLES_ID_LIST,
                clade=STRAINPHLAN_CLADES_LIST_V4)

else:
    rule profiling_strainphlan4_all:
        input: