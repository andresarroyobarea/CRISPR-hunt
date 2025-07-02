rule cnv_correction:
    input: 
        count_table_raw="results/mageck_count_raw/raw.counts.txt",
    output:
        count_table_cnv_norm="results/cnv_correction/norm_cnv_correction.counts.txt",
    conda:
        config["conda_envs"]["crisprcleanR"],
    params:
        sgrna_library=config["parameters"]["cnv_correction"]["sgrna_library"],
        outdir=outdir=lambda wildcards, output: os.path.dirname(output.count_table),
        min_reads=config["parameters"]["cnv_correction"]["min_reads"],
        ncontrols=len(ctrl_samples),
        norm_method=config["parameters"]["cnv_correction"]["norm_method"],
        min_genes=config["parameters"]["cnv_correction"]["min_genes"],
        project=config["project"],
        extra=config["parameters"]["cnv_correction"]["extra"],
    log:
        "logs/cnv_correction/cnv_correction.log",
    benchmark:
        "benchmarks/cnv_correction/cnv_correction.bmk",
    shell: 
        "Rscript "