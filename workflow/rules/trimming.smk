# TODO: cutadapt paired-end
rule trimming:
    input:
        fastq="data/{sample}.fastq.gz",
    output:
        fastq_trimmed="results/trimming/{sample}_trimmed.fastq",
        trimming_qc="results/trimming/{sample}_trimming.qc.txt",
    threads: get_resource("trimming", "threads")
    resources:
        mem_mb=get_resource("trimming", "mem_mb"),
        runtime=get_resource("trimming", "runtime"),
    params:
        extra=config["parameters"]["trimming"]["extra"],
        adapters=config["parameters"]["trimming"]["adapters"],
    log:
        "logs/trimming/{sample}.log",
    benchmark:
        "benchmarks/trimming/{sample}.bmk"
    wrapper:
        config["wrapper"]["trimming"]
