# TODO: Include bowtie2 and/or STAR alignment.
rule mageck_count_norm:
    input:
        fastqs_trimmed=expand("results/trimming/{sample}_trimmed.fastq", sample=SAMPLES),
    output:
        count_table="results/mageck_count_norm/{norm_state}.counts.txt",
        count_summary="results/mageck_count_norm/{norm_state}.countsummary.txt",
    conda:
        config["conda_envs"]["mageck"]
    params:
        library_file=config["library_file"],
        norm_method=config["parameters"]["mageck_count_norm"]["norm_method"],
        out_prefix=config["parameters"]["mageck_count_norm"]["out_prefix"],
        sample_labels=",".join(SAMPLES),
        extra=config["parameters"]["mageck_count_norm"]["extra"],
    log:
        "logs/mageck_count_norm/mageck_count_norm.log",
    benchmark:
        "benchmarks/mageck_count_norm/mageck_count_norm.bmk"
    shell:
        """
        mageck count --list-seq {params.library_file} \
            --fastq {input.fastqs_trimmed} \
            --count-table {output.count_table} \
            --norm-method {params.norm_method} \
            --sample-label {params.sample_labels} \
            --output-prefix {params.out_prefix} \
            {params.extra} 2> {log}
        """