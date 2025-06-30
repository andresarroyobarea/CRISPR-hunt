# TODO: Include bowtie2 and/or STAR alignment.
rule mageck_count_raw:
    input:
        fastqs_trimmed=expand("results/trimming/{sample}_trimmed.fastq", sample=SAMPLES),
    output:
        count_table="results/mageck_count_raw/raw.counts.txt",
        count_summary="results/mageck_count_raw/raw.countsummary.txt",
    conda:
        config["conda_envs"]["mageck"]
    params:
        library_file=config["library_file"],
        sample_labels=",".join(SAMPLES),
        extra=config["parameters"]["mageck_count_raw"]["extra"],
    log:
        "logs/mageck_count_raw/mageck_count_raw.log",
    benchmark:
        "benchmarks/mageck_count_raw/mageck_count_raw.bmk"
    shell:
        """
        mageck count --list-seq {params.library_file} \
            --fastq {input.fastqs_trimmed} \
            --count-table {output.count_table} \
            --norm-method "none" \
            --sample-label {params.sample_labels} \
            --output-prefix raw \
            {params.extra} 2> {log}
        """