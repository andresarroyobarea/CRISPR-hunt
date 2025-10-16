rule qc_raw: 
    input: 
        fastq="data/{sample}.fastq.gz", 
    output: 
        html="results/qc_raw/{sample}/{sample}_fastqc.html", 
        zip="results/qc_raw/{sample}/{sample}_fastqc.zip", 
    resources: 
        mem_mb=get_resource("qc_raw", "mem_mb"), 
        runtime=get_resource("qc_raw", "runtime"), 
    params: 
        outdir=lambda wildcards, output: os.path.dirname(output.html), lambda wc: "-t {}".format(get_resource("qc_raw", "threads")), 
    log: 
        "logs/qc_raw/{sample}.log", 
    benchmark: 
        "benchmarks/qc_raw/{sample}.bmk" 
    wrapper: 
        config["wrapper"]["qc"]

rule qc_trimming: 
    input: 
        fastq_trimmed="results/trimming/{sample}_trimmed.fastq", 
    output: 
        html="results/qc_trimming/{sample}/{sample}_trimmed_fastqc.html", 
        zip="results/qc_trimming/{sample}/{sample}_trimmed_fastqc.zip", 
    resources: 
        mem_mb=get_resource("qc_trimming", "mem_mb"), 
        runtime=get_resource("qc_trimming", "runtime"), 
    params: 
        outdir=lambda wildcards, output: os.path.dirname(output.html), lambda wc: "-t {}".format(get_resource("qc_trimming", "threads")), 
    log: 
        "logs/qc_trimming/{sample}.log", 
    benchmark: 
        "benchmarks/qc_trimming/{sample}.bmk" 
    wrapper: 
        config["wrapper"]["qc"]