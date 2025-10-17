rule cnv_correction: 
    input: 
        count_table_raw="results/filter_sgrna_counts/{project}_processed.count.txt", 
    output: 
        count_table_cnv_norm="results/cnv_correction/{project}_cnvcorr.count.txt", 
        lfc_cnv_norm="results/cnv_correction/{project}_cnvcorr.foldchange",
    conda: 
        config["conda_envs"]["crisprcleanR"] 
    params: 
        library_type=config["parameters"]["cnv_correction"]["library_type"],
        sgrna_library=config["parameters"]["cnv_correction"]["sgrna_library"], 
        norm_method=config["parameters"]["cnv_correction"]["norm_method"], 
        min_reads=config["parameters"]["cnv_correction"]["min_reads"],
        min_genes=config["parameters"]["cnv_correction"]["min_genes"], 
        treat_samples=",".join(treat_samples), 
        ctrl_samples=",".join(ctrl_samples), 
        exp_design=config["parameters"]["cnv_correction"]["exp_design"], 
        project=config["project"], outdir=lambda wildcards, output: os.path.dirname(output.count_table_cnv_norm), 
        extra=config["parameters"]["cnv_correction"]["extra"], 
    log: 
        "logs/cnv_correction/{project}_cnvcorr.log", 
    benchmark: 
        "benchmarks/cnv_correction/{project}_cnvcorr.bmk" 
    shell: 
        """ 
        Rscript workflow/scripts/run_crisprcleanR.R \
            --input {input.count_table_raw} \
            --out-count {output.count_table_cnv_norm} \
            --out-lfc {output.lfc_cnv_norm} \
            --lib-type {params.library_type} \
            --sgrna-library {params.sgrna_library} \
            --norm-method {params.norm_method} \
            --exp-design {params.exp_design} \
            --min-reads {params.min_reads} \
            --min-genes {params.min_genes} \
            --control-samples {params.ctrl_samples} \
            --treat-samples {params.treat_samples} \
            --label {params.project} \
            --outdir {params.outdir} \
            {params.extra} &> {log} 
        """