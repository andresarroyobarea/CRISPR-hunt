rule qc_raw:
    input:
        fastq="data/{sample}.fastq.gz",
    output:
        html="results/qc_raw/{sample}/{sample}_fastqc.html",
        zip="results/qc_raw/{sample}/{sample}_fastqc.zip",
    threads: get_resource("qc_raw", "threads")
    resources:
        mem_mb=get_resource("qc_raw", "mem_mb"),
        runtime=get_resource("qc_raw", "runtime"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.html),
    log:
        fastqc="logs/qc_raw/{sample}.log",
    benchmark:
        "benchmarks/qc_raw/{sample}.bmk"
    wrapper:
        config["wrapper"]["qc"]
