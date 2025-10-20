rule qc_raw:
    input:
        fastq="data/{sample}.fastq.gz",
    output:
        html="results/qc_raw/{sample}/{sample}_fastqc.html",
        zip="results/qc_raw/{sample}/{sample}_fastqc.zip",
    resources:
        mem_mb=get_resource(config, "qc", "mem_mb"),
        runtime=get_resource(config, "qc", "runtime"),
    params:
        lambda wc: "-t {}".format(get_resource(config, "qc", "threads")),
        outdir=lambda wildcards, output: os.path.dirname(output.html),
    log:
        "logs/qc_raw/{sample}.log",
    benchmark:
        "benchmarks/qc_raw/{sample}.bmk"
    wrapper:
        config["wrapper"]["qc"]


rule qc_trimming:
    input:
        fastq_trimmed="results/trimming/{sample}_trimmed.fastq.gz",
    output:
        html="results/qc_trimming/{sample}/{sample}_trimmed_fastqc.html",
        zip="results/qc_trimming/{sample}/{sample}_trimmed_fastqc.zip",
    resources:
        mem_mb=get_resource(config, "qc", "mem_mb"),
        runtime=get_resource(config, "qc", "runtime"),
    params:
        lambda wc: "-t {}".format(get_resource(config, "qc", "threads")),
        outdir=lambda wildcards, output: os.path.dirname(output.html),
    log:
        "logs/qc_trimming/{sample}.log",
    benchmark:
        "benchmarks/qc_trimming/{sample}.bmk"
    wrapper:
        config["wrapper"]["qc"]


rule multiqc_raw:
    input:
        fastqc=expand("results/qc_raw/{sample}/{sample}_fastqc.html", sample=SAMPLES),
    output:
        multiqc_report="results/qc_raw/multiqc_report.html",
    params:
        extra="--verbose",
        outdir=lambda wc, output: os.path.dirname(output.multiqc_report),
    log:
        "logs/qc_raw/multiqc_raw.log",
    benchmark:
        "benchmarks/qc_raw/multiqc_raw.bmk"
    wrapper:
        config["wrapper"]["multiqc"]


rule multiqc_trimmed:
    input:
        fastqc_trimmed=expand("results/qc_trimming/{sample}/{sample}_trimmed_fastqc.html", sample=SAMPLES),
        # TODO: Check if this log file is ok for cutadapt.
        cutadapt_stats=expand(
            "logs/trimming/{sample}.log",
            sample=SAMPLES
        ),
    output:
        multiqc_report="results/qc_trimming/multiqc_report.html",
    params:
        extra="--verbose",
        outdir=lambda wc, output: os.path.dirname(output.multiqc_report),
    log:
        "logs/qc_trimming/multiqc_trimming.log",
    benchmark:
        "benchmarks/qc_trimming/multiqc_trimming.bmk"
    wrapper:
        config["wrapper"]["multiqc"]
