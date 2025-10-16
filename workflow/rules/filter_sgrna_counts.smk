rule filter_sgrna_counts: 
    input: 
        raw_counts="results/mageck_count_raw/{project}.count.txt", 
    output: 
        processed_count="results/filter_sgrna_counts/{project}_processed.count.txt",
        processed_count_NCT="results/filter_sgrna_counts/{project}_processed_NCT.count.txt", 
        filtering_report="results/filter_sgrna_counts/{project}_filter_report.txt", 
        sgrna_removed = "results/filter_sgrna_counts/{project}_removed_sgRNAs.txt", 
    conda: 
        config["conda_envs"]["sgrna_prepr"], 
    params: 
        project=config["project"], 
        ctrl_samples=",".join(ctrl_samples), 
        treat_samples=",".join(treat_samples),
        extra=config["parameters"]["filter_sgrna_counts"]["extra"], 
    log: 
        "logs/filter_sgrna_counts/{project}_{mageck_norms}_filter_sgrna_counts.log", 
    benchmark: 
        "benchmarks/filter_sgrna_counts/{project}_{mageck_norms}_filter_sgrna_counts.bmk" 
    shell: 
        """ 
        Rscript workflow/scripts/sgRNA_preprocessing.R \
            --input {input.raw_counts} \
            --out-count-filt {output.processed_count} \
            --out-count-filt-NCT {output.processed_count_NCT} \
            --out-filter-report {output.filtering_report} \
            --out-sgrna-removed {output.sgrna_removed} \
            --control-samples {params.ctrl_samples} \
            --treat-samples {params.treat_samples} \
            --project {wildcards.project} \ {params.extra} &> {log} 
        """