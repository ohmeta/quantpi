## reference
## https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-3


if config["params"]["profiling"]["strainphlan"]["do_v3"]:
    rule profiling_strainphlan3_sample2markers:
        input:
            sam = os.path.join(config["output"]["profiling"],
                               "profile/metaphlan3/{sample}/{sample}.sam.bz2"),
            aln = os.path.join(config["output"]["profiling"],
                               "profile/metaphlan3/{sample}/{sample}.bowtie2.bz2")
        output:
            os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/{sample}/consensus_markers/{sample}.pkl")
        log:
            os.path.join(
                config["output"]["profiling"],
                "logs/strainphlan3_sample2markers/{sample}.strainphlan3_sample2markers.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"],
                "benchmark/strainphlan3_sample2markers/{sample}.strainphlan3_sample2markers.benchmark.txt")
        conda:
            config["envs"]["biobakery3"]
        params:
            outdir = os.path.join(config["output"]["profiling"], "profile/strainphlan3/{sample}/consensus_markers")
        priority:
            20
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            mkdir -p {params.outdir}

            sample2markers.py \
            --input {input.sam} \
            --output_dir {params.outdir}/ \
            --nprocs {threads} \
            > {log} 2>&1
            '''


    STRAINPHLAN_CLADES_V3 = \
        pd.read_csv(config["params"]["profiling"]["strainphlan"]["clades_tsv_v3"], sep="\t")\
          .set_index("clades")
    STRAINPHLAN_CLADES_LIST_V3 = STRAINPHLAN_CLADES_V3.index.unique()


    rule profiling_strainphlan3_extract_markers:
        input:
            clades_tsv = config["params"]["profiling"]["strainphlan"]["clades_tsv"],
            database_pkl = expand(os.path.join(
                config["params"]["profiling"]["metaphlan"]["bowtie2db"], "{index}.pkl"),
                index = config["params"]["profiling"]["metaphlan"]["index_v3"])
        output:
            clade_marker = os.path.join(
                config["output"]["profiling"],
                "databases/strainphlan3/clade_markers/{clade}.fna")
        log:
            os.path.join(
                config["output"]["profiling"], 
                "logs/strainphlan3_extract_markers/{clade}.strainphlan3_extract_markers.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"],
                "benchmark/strainphlan3_extract_markers/{clade}.strainphlan3_extract_markers.benchmark.txt")
        conda:
            config["envs"]["biobakery3"]
        params:
            clade = "{clade}",
            outdir = os.path.join(config["output"]["profiling"], "databases/strainphlan3/clade_markers")
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


    rule profiling_strainphlan3:
        input:
            clade_marker = os.path.join(
                config["output"]["profiling"],
                "databases/strainphlan3/clade_markers/{clade}.fna")
            consensus_markers = expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/{sample}/consensus_markers/{sample}.pkl"),
                sample=SAMPLES_ID_LIST),
            reference_genome = lambda wildcards: STRAINPHLAN_CLADES_V3.loc[wildcards.clade, "fna_path"]
        output:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/{{sample}}/{{clade}}/RAxML_{prefix}.{{clade}}.StrainPhlAn3.tre"),
                prefix=["bestTree", "info", "log", "parsimonyTree", "result"]),
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/{{sample}}/{{clade}}/{{clade}}{suffix}"),
                suffix=["_mutation_rates", ".info", ".mutation", ".polymorphic",
                        ".StrainPhlAn3_concatenated.aln"])
        log:
            os.path.join(
                config["output"]["profiling"],
                "logs/strainphlan3/{sample}.{clade}.strainphlan3.log")
        benchmark:
            os.path.join(
                config["output"]["profiling"],
                "benchmark/strainphlan3/{sample}.{clade}.strainphlan3.benchmark.txt")
        conda:
            config["envs"]["biobakery3"]
        params:
            clade = "{clade}",
            outdir = os.path.join(config["output"]["profiling"], "profile/strainphlan3/{sample}/{clade}"),
            opts = config["params"]["profiling"]["strainphlan"]["external_opts_v3"]
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


    rule profiling_strainphlan3_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/{sample}/{clade}/RAxML_{prefix}.{clade}.StrainPhlAn3.tre"),
                prefix=["bestTree", "info", "log", "parsimonyTree", "result"],
                sample=SAMPLES_ID_LIST,
                clade=STRAINPHLAN_CLADES_LIST_V3),
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/strainphlan3/{sample}/{clade}/{clade}{suffix}"),
                suffix=["_mutation_rates", ".info", ".mutation", ".polymorphic",
                        ".StrainPhlAn3_concatenated.aln"],
                sample=SAMPLES_ID_LIST,
                clade=STRAINPHLAN_CLADES_LIST_V3)

else:
    rule profiling_strainphlan3_all:
        input: