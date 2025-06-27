# TODO: Include bowtie2 and/or STAR alignment.
rule alignment:
    input:
        fastqs_trimmed = expand("results/trimming/{sample}_trimmed.fastq", sample=SAMPLES),
    output:
        count_table="results/alignment/{project}.counts.txt",
        count_summary="results/alignment/{project}.countsummary.txt",
    conda:
        config["conda_envs"]["mageck"],
    params:
        library_file=config["library_file"],
        norm_method=config["parameters"]["alignment"]["norm_method"],
        project=config["project"],
        sample_labels= ",".join(SAMPLES),
        extra=config["parameters"]["alignment"]["extra"],
    log:
        "logs/alignment/alignment.log",
    benchmark:
        "benchmarks/alignment/alignment.bmk"
    shell:
        """
        mageck count --list-seq {params.library_file} \
            --fastq {input.fastqs_trimmed} \
            --count-table {output.count_table} \
            --norm-method {params.norm_method} \
            --sample-label {params.sample_labels} \
            --output-prefix {params.project} \
            {params.extra}
        """