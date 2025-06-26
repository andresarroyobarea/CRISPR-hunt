# TODO: Include bowtie2 and/or STAR alignment.
rule alignment:
    input:
        expand("results/trimming/{sample}_trimmed.fastq", sample=SAMPLES)
    output:
        count_table="results/alignment/{project}.counts.txt",
        count_summary="results/alignment/{project}.countsummary.txt",
    conda:
        config["conda_envs"]["mageck"]
    params:
        library_file=config["library_file"],
        # nombre proyecto
        # muestras control
        # muetra tratamiento o dia x
        # NORM method
        # Extra
    log:
        "logs/alignment/{sample}.log",
    benchmark:
        "logs/alignment/{sample}.bmk"
    shell:
        """
        mageck count -l {params.library_file} \
            -n {params.project} # CONFIRMAR
            --fastq
            -k
            --norm-method 
        """