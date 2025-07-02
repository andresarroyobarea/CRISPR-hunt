rule cnv_correction:
    input: 
        count_table_raw="results/mageck_count_raw/raw.counts.txt",
    output:
        count_table_cnv_norm="results/cnv_correction/norm_cnv_correction.counts.txt",
    conda:
        config["conda_envs"]["crisprcleanR"],
    params:
        library_type=config["parameters"]["cnv_correction"]["library_type"],
        sgrna_library=config["parameters"]["cnv_correction"]["sgrna_library"],
        norm_method=config["parameters"]["cnv_correction"]["norm_method"],
        min_reads=config["parameters"]["cnv_correction"]["min_reads"],
        min_genes=config["parameters"]["cnv_correction"]["min_genes"],
        ncontrols=len(ctrl_samples),
        project=config["project"],
        outdir=outdir=lambda wildcards, output: os.path.dirname(output.count_table),
        extra=config["parameters"]["cnv_correction"]["extra"],
    log:
        "logs/cnv_correction/cnv_correction.log",
    benchmark:
        "benchmarks/cnv_correction/cnv_correction.bmk",
    shell: 
        """
        Rscript scripts/run_crisprcleanR.R \
            --input {input.count_table_raw} \
            --output {output.count_table_cnv_norm} \
            --lib-type {params.library_type} \
            --sgrna-library {params.sgrna_library} \
            --norm-method {params.norm_method} \
            --min-reads {params.min_reads} \
            --min-genes {params.min_genes} \
            --n-control-samples {params.ncontrols} \
            --label {params.project} \
            --outdir {params.outdir} \
            {params.extra} 2> {log}
        """