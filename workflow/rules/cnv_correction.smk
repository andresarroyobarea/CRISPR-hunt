rule cnv_correction:
    input: 
        count_table_raw="results/mageck_count_raw/raw.counts.txt",
    output:
        count_table = "results/cnv_correction/norm_cnv_correction.counts.txt"
    conda:
        config["conda_envs"]["crisprcleanR"],
    params:
        sgrna_library=config["parameters"]["cnv_correction"]["sgrna_library"],
        extra=config["parameters"]["cnv_correction"]["extra"],
    log:
        "logs/cnv_correction/cnv_correction.log",
    benchmark:
        "benchmarks/cnv_correction/cnv_correction.bmk",
    shell: 
        "Rscript "