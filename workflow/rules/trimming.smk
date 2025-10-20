# TODO: cutadapt paired-end
rule trimming:
    input:
        "data/{sample}.fastq.gz",
    output:
        fastq="results/trimming/{sample}_trimmed.fastq.gz",
        qc="results/trimming/{sample}_trimmed.qc.txt",
    threads: get_resource(config, "trimming", "threads")
    resources:
        mem_mb=get_resource(config, "trimming", "mem_mb"),
        runtime=get_resource(config, "trimming", "runtime"),
    params:
        adapters=config["parameters"]["trimming"]["adapters"],
        extra=config["parameters"]["trimming"]["extra"],
    log:
        "logs/trimming/{sample}.log",
    benchmark:
        "benchmarks/trimming/{sample}.bmk"
    wrapper:
        config["wrapper"]["trimming"]
