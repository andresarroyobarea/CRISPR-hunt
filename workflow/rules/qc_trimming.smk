rule qc_trimming:
    input:
        fastq_trimmed="results/trimming/{sample}_trimmed.fastq",
    output:
        html="results/qc_trimming/{sample}/{sample}_trimmed_fastqc.html",
        zip="results/qc_trimming/{sample}/{sample}_trimmed_fastqc.zip",
    threads: get_resource("qc_trimming", "threads")
    resources:
        mem_mb=get_resource("qc_trimming", "mem_mb"),
        runtime=get_resource("qc_trimming", "runtime"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.html),
    log:
        "logs/qc_trimming/{sample}.log",
    benchmark:
        "benchmarks/qc_trimming/{sample}.bmk"
    wrapper:
        config["wrapper"]["qc"]
