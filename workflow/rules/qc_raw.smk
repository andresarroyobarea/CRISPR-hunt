rule qc_raw:
    input:
        fastq="data/{sample}.fastq.gz",
    output:
        html="results/qc_raw/{sample}/{sample}_fastqc.html",
        zip="results/qc_raw/{sample}/{sample}_fastqc.zip",
    conda:
        config["conda_envs"]["qc"]
    threads: get_resource("qc_raw", "threads")
    resources:
        mem_mb=get_resource("qc_raw", "mem_mb"),
        runtime=get_resource("qc_raw", "runtime"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.html),
    log:
        fastqc="log/qc_raw/{sample}_fastqc.log",
    shell:
        "mkdir -p {params.outdir} &&"
        "fastqc --outdir {params.outdir} --threads {threads} {input.fastq} 2> {log.fastqc} "
