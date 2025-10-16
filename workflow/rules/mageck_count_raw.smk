rule mageck_count_raw:
    input:
        fastqs_trimmed=expand("results/trimming/{sample}_trimmed.fastq.gz", sample=SAMPLES),
    output:
        count_table="results/mageck_count_raw/{project}.count.txt",
        count_summary="results/mageck_count_raw/{project}.countsummary.txt",
    conda:
        config["conda_envs"]["mageck"]
    params:
        library_file=config["library_file"],
        label=lambda wc: f"results/mageck_count_raw/{wc.project}",
        sample_labels=",".join(SAMPLES),
        #sample_labels=lambda wc: ",".join(config["project_samples"][wc.project]),
        extra=config["parameters"]["mageck_count_raw"]["extra"],
    log:
        "logs/mageck_count_raw/{project}_mageck_count_raw.log",
    benchmark:
        "benchmarks/mageck_count_raw/{project}_mageck_count_raw.bmk"
    shell:
        """
        mageck count --list-seq "{params.library_file}" \
            --fastq {input.fastqs_trimmed} \
            --norm-method none \
            --sample-label "{params.sample_labels}" \
            --output-prefix "{params.label}" \
            {params.extra} 2> "{log}"
        """
